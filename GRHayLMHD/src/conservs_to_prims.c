/* We evolve forward in time a set of functions called the
 *  "conservative variables", and any time the conserv's
 *  are updated, we must solve for the primitive variables
 *  (rho, pressure, velocities) using a Newton-Raphson
 *  technique, before reconstructing & evaluating the RHSs
 *  of the MHD equations again.
 *
 * This file contains the driver routine for this Newton-
 *  Raphson solver. Truncation errors in conservative
 *  variables can lead to no physical solutions in
 *  primitive variables. We correct for these errors here
 *  through a number of tricks described in the appendices
 *  of http://arxiv.org/pdf/1112.0568.pdf.
 *
 * This is a wrapper for the 2d solver of Noble et al. See
 *  harm_utoprim_2d.c for references and copyright notice
 *  for that solver. This wrapper was primarily written by
 *  Zachariah Etienne & Yuk Tung Liu, in 2011-2013.
 *
 * For optimal compatibility, this wrapper is licensed under
 *  the GPL v2 or any later version.
 *
 * Note that this code assumes a simple gamma law for the
 *  moment, though it would be easy to extend to a piecewise
 *  polytrope. */

#include "GRHayLMHD.h"
#include "Symmetry.h"

void GRHayLMHD_conserv_to_prims(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLMHD_conserv_to_prims;
  DECLARE_CCTK_PARAMETERS;

  const int imax = cctk_lsh[0];
  const int jmax = cctk_lsh[1];
  const int kmax = cctk_lsh[2];

  if(perturb_every_con2prim > 1e-30) {
    srand(random_seed); // Use srand() as rand() is thread-safe.
#pragma omp parallel for
    for(int k=0; k<kmax; k++) {
      for(int j=0; j<jmax; j++) {
        for(int i=0; i<imax; i++) {
          const int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
          rho_star[index] *= one_plus_pert(perturb_every_con2prim);
          tau[index]      *= one_plus_pert(perturb_every_con2prim);
          Stildex[index]  *= one_plus_pert(perturb_every_con2prim);
          Stildey[index]  *= one_plus_pert(perturb_every_con2prim);
          Stildez[index]  *= one_plus_pert(perturb_every_con2prim);
        }
      }
    }
  }

  if(CCTK_EQUALS(Symmetry,"equatorial")) {
    // SET SYMMETRY GHOSTZONES ON ALL CONSERVATIVE VARIABLES!
    int ierr=0;
    ierr+=CartSymGN(cctkGH,"GRHayLMHD::grmhd_conservatives");
    // FIXME: UGLY. Filling metric ghostzones is needed for, e.g., Cowling runs.
    ierr+=CartSymGN(cctkGH,"lapse::lapse_vars");
    ierr+=CartSymGN(cctkGH,"bssn::BSSN_vars");
    ierr+=CartSymGN(cctkGH,"bssn::BSSN_AH");
    ierr+=CartSymGN(cctkGH,"shift::shift_vars");
    if(ierr!=0) CCTK_VERROR("Error with setting equatorial symmetries in con2prim.");
  }

  //Start the timer, so we can benchmark the primitives solver during evolution.
  //  Slower solver -> harder to find roots -> things may be going crazy!
  //FIXME: Replace this timing benchmark with something more meaningful, like the avg # of Newton-Raphson iterations per gridpoint!
  /*
    struct timeval start, end;
    long mtime, seconds, useconds;
    gettimeofday(&start, NULL);
  */

  const double poison = 0.0/0.0;

  // Diagnostic variables.
  int failures=0;
  int vel_limited_ptcount=0;
  int rho_star_fix_applied=0;
  int pointcount=0;
  int failures_inhoriz=0;
  int pointcount_inhoriz=0;
  int backup0=0;
  int backup1=0;
  int backup2=0;
  double error_int_numer=0;
  double error_int_denom=0;
  int n_iter=0;
  double dummy1, dummy2, dummy3;

#pragma omp parallel for reduction(+:failures,vel_limited_ptcount,rho_star_fix_applied,pointcount,failures_inhoriz,pointcount_inhoriz,backup0,backup1,backup2,error_int_numer,error_int_denom,n_iter) schedule(static)
  for(int k=0; k<kmax; k++) {
    for(int j=0; j<jmax; j++) {
      for(int i=0; i<imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        ghl_con2prim_diagnostics diagnostics;
        ghl_initialize_diagnostics(&diagnostics);

        // Read in ADM metric quantities from gridfunctions and
        // set auxiliary and ADM metric quantities
        ghl_metric_quantities ADM_metric;
        ghl_enforce_detgtij_and_initialize_ADM_metric(
              alp[index],
              betax[index], betay[index], betaz[index],
              gxx[index], gxy[index], gxz[index],
              gyy[index], gyz[index], gzz[index],
              &ADM_metric);

        ghl_ADM_aux_quantities metric_aux;
        ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

        // Read in primitive variables from gridfunctions
        ghl_primitive_quantities prims;
        ghl_initialize_primitives(
              rho_b[index], pressure[index], eps[index],
              vx[index], vy[index], vz[index],
              Bx_center[index], By_center[index], Bz_center[index],
              poison, poison, poison, &prims);

        // Read in conservative variables from gridfunctions
        ghl_conservative_quantities cons, cons_orig;
        ghl_initialize_conservatives(
              rho_star[index], tau[index],
              Stildex[index], Stildey[index], Stildez[index],
              poison, poison, &cons);

        // Here we save the original values of conservative variables in cons_orig for debugging purposes.
        cons_orig = cons;

        //FIXME: might slow down the code. Was formerly a CCTK_WARN
        if(isnan(cons.rho*cons.tau*cons.SD[0]*cons.SD[1]*cons.SD[2]*prims.BU[0]*prims.BU[1]*prims.BU[2])) {
          CCTK_VERROR("NaN found at start of C2P kernel:\n"
                      "  index %d %d %d, rho_* = %e, ~tau = %e, ~S_i = %e %e %e, Bi = %e %e %e\n"
                      "  lapse = %e, shift = %e %e %e, gij = %e %e %e %e %e %e, Psi6 = %e\n",
                      i, j, k, cons.rho, cons.tau, cons.SD[0], cons.SD[1], cons.SD[2], prims.BU[0], prims.BU[1], prims.BU[2],
                      ADM_metric.lapse, ADM_metric.betaU[0], ADM_metric.betaU[1], ADM_metric.betaU[2],
                      ADM_metric.gammaDD[0][0], ADM_metric.gammaDD[0][1], ADM_metric.gammaDD[0][2],
                      ADM_metric.gammaDD[1][1], ADM_metric.gammaDD[1][2], ADM_metric.gammaDD[2][2], ADM_metric.sqrt_detgamma);
        }

        /************* Main conservative-to-primitive logic ************/
        if(cons.rho>0.0) {
          // Apply the tau floor
          if( ghl_eos->eos_type == ghl_eos_hybrid )
            ghl_apply_conservative_limits(
                  ghl_params, ghl_eos, &ADM_metric,
                  &prims, &cons, &diagnostics);

          // declare some variables for the C2P routine.
          ghl_conservative_quantities cons_undens;

          // Set the conserved variables required by the con2prim routine
          ghl_undensitize_conservatives(ADM_metric.sqrt_detgamma, &cons, &cons_undens);

          /************* Conservative-to-primitive recovery ************/
          const int check = ghl_con2prim_multi_method(
                ghl_params, ghl_eos, &ADM_metric, &metric_aux,
                &cons_undens, &prims, &diagnostics);

          if(check==0) {
            //Check for NAN!
            if( isnan(prims.rho*prims.press*prims.eps*prims.vU[0]*prims.vU[1]*prims.vU[2]) ) {
              CCTK_VERROR("***********************************************************\n"
                          "NAN found after Con2Prim with routine index %d!\n"
                          "Input variables:\n"
                          "lapse, shift = %e %e %e %e\n"
                          "gij = %e %e %e %e %e %e\n"
                          "rho_*, ~tau, ~S_{i}: %e %e %e %e %e\n"
                          "Undensitized conserved variables:\n"
                          "D, tau, S_{i}: %e %e %e %e %e\n"
                          "Output primitive variables:\n"
                          "rho, P: %e %e\n"
                          "v: %e %e %e\n"
                          "***********************************************************",
                          diagnostics.which_routine, ADM_metric.lapse, ADM_metric.betaU[0], ADM_metric.betaU[1], ADM_metric.betaU[2],
                          ADM_metric.gammaDD[0][0], ADM_metric.gammaDD[0][1], ADM_metric.gammaDD[0][2],
                          ADM_metric.gammaDD[1][1], ADM_metric.gammaDD[1][2], ADM_metric.gammaDD[2][2],
                          cons.rho, cons.tau, cons.SD[0], cons.SD[1], cons.SD[2],
                          cons_undens.rho, cons_undens.tau, cons_undens.SD[0], cons_undens.SD[1], cons_undens.SD[2],
                          prims.rho, prims.press,
                          prims.vU[0], prims.vU[1], prims.vU[2]);
            }
          } else {
            //--------------------------------------------------
            //----------- Primitive recovery failed ------------
            //--------------------------------------------------
            failure_checker[index] += 100;
            ghl_set_prims_to_constant_atm(ghl_eos, &prims);

            failures++;
            if(ADM_metric.sqrt_detgamma > ghl_params->psi6threshold) {
              failures_inhoriz++;
              pointcount_inhoriz++;
            }
            CCTK_VINFO("Con2Prim failed! Resetting to atmosphere...\n");
            CCTK_VINFO("rho_* = %e, ~tau = %e, ~S_i = %e %e %e, Bi = %e %e %e,\n"
                       "lapse = %e, shift = %e %e %e, gij = %e %e %e %e %e %e, Psi6 = %e",
                       cons_orig.rho, cons_orig.tau, cons_orig.SD[0], cons_orig.SD[1], cons_orig.SD[2], prims.BU[0], prims.BU[1], prims.BU[2],
                       ADM_metric.lapse, ADM_metric.betaU[0], ADM_metric.betaU[1], ADM_metric.betaU[2],
                       ADM_metric.gammaDD[0][0], ADM_metric.gammaDD[0][1], ADM_metric.gammaDD[0][2],
                       ADM_metric.gammaDD[1][1], ADM_metric.gammaDD[1][2], ADM_metric.gammaDD[2][2], ADM_metric.sqrt_detgamma);
          }
        } else {
          failure_checker[index] += 1;
          ghl_set_prims_to_constant_atm(ghl_eos, &prims);
          rho_star_fix_applied++;
        } // if rho_star>0
        /***************************************************************/


        //--------------------------------------------------
        //---------- Primitive recovery succeeded ----------
        //--------------------------------------------------
        // Enforce limits on primitive variables and recompute conservatives.
        diagnostics.speed_limited += ghl_enforce_primitive_limits_and_compute_u0(
              ghl_params, ghl_eos, &ADM_metric, &prims);
        ghl_compute_conservs(
              &ADM_metric, &metric_aux, &prims, &cons);

        //Now we compute the difference between original & new conservatives, for diagnostic purposes:
        error_int_numer += fabs(cons.tau - cons_orig.tau) + fabs(cons.rho - cons_orig.rho) + fabs(cons.SD[0] - cons_orig.SD[0])
                           + fabs(cons.SD[1] - cons_orig.SD[1]) + fabs(cons.SD[2] - cons_orig.SD[2]);
        error_int_denom += cons_orig.tau + cons_orig.rho + fabs(cons_orig.SD[0]) + fabs(cons_orig.SD[1]) + fabs(cons_orig.SD[2]);

        ghl_return_primitives(
              &prims,
              &rho_b[index], &pressure[index], &eps[index],
              &vx[index], &vy[index], &vz[index],
              &Bx_center[index], &By_center[index], &Bz_center[index],
              &dummy1, &dummy2, &dummy3);
        u0[index] = prims.u0;

        ghl_return_conservatives(
              &cons,
              &rho_star[index], &tau[index],
              &Stildex[index], &Stildey[index], &Stildez[index],
              &dummy1, &dummy2);

        pointcount++;
        if(diagnostics.speed_limited) {
          failure_checker[index] += 10;
          vel_limited_ptcount++;
        }
        backup0 += diagnostics.backup[0];
        backup1 += diagnostics.backup[1];
        backup2 += diagnostics.backup[2];
        n_iter += diagnostics.n_iter;
        failure_checker[index] += 1000*diagnostics.backup[0]
                                + 10000*diagnostics.tau_fix
                                + 100000*diagnostics.Stilde_fix;
      }
    }
  }

  /*
    Failure checker decoder:
       1: atmosphere reset when rho_star < 0
      10: Limiting velocity u~ after C2P/Font Fix or v in ghl_enforce_primitive_limits_and_compute_u0
     100: Both C2P and Font Fix failed
      1k: backups used
     10k: tau~ was reset in ghl_apply_conservative_limits
    100k: S~ was reset in ghl_apply_conservative_limits
  */
  if(CCTK_Equals(verbose, "essential") || CCTK_Equals(verbose, "essential+iteration output")) {
    CCTK_VINFO("C2P: Lev: %d NumPts= %d | Fixes: Font= %d VL= %d rho*= %d | Failures: %d InHoriz= %d / %d | Error: %.3e, ErrDenom: %.3e | %.2f iters/gridpt",
               (int)GetRefinementLevel(cctkGH), pointcount,
               backup0, vel_limited_ptcount, rho_star_fix_applied,
               failures, failures_inhoriz, pointcount_inhoriz,
               error_int_numer/error_int_denom, error_int_denom,
               (double)n_iter/( (double)(cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]) ));
  }
}
