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

#include "GRHayLHDX.h"

extern "C" void GRHayLHDX_conservs_to_prims(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_conservs_to_prims;
  DECLARE_CCTK_PARAMETERS;

  // Diagnostic variables.
  int failures = 0;
  int vel_limited_ptcount = 0;
  int rho_star_fix_applied = 0;
  int pointcount = 0;
  int failures_inhoriz = 0;
  int pointcount_inhoriz = 0;
  int backup0 = 0;
  int backup1 = 0;
  int backup2 = 0;
  double error_int_numer = 0;
  double error_int_denom = 0;
  int n_iter = 0;

  const Loop::GF3D2layout layout(cctkGH, {1,1,1});

//I don't think I can do this sort of diagnostic summing in carpetx
//#pragma omp parallel for reduction(+: \
//      pointcount, backup0, backup1, backup2, vel_limited_ptcount, rho_star_fix_applied, failures, failures_inhoriz, pointcount_inhoriz, n_iter, \
//      error_rho_numer, error_tau_numer, error_Sx_numer, error_Sy_numer, error_Sz_numer, error_entropy_numer, error_Ye_numer,                    \
//      error_rho_denom, error_tau_denom, error_Sx_denom, error_Sy_denom, error_Sz_denom, error_entropy_denom, error_Ye_denom) schedule(static)
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index index(layout, p.I);

    double local_failure_checker = 0;

    ghl_con2prim_diagnostics diagnostics;
    ghl_initialize_diagnostics(&diagnostics);

    // Read in ADM metric quantities from gridfunctions and
    // set auxiliary and ADM metric quantities
    ghl_metric_quantities ADM_metric;
    ghl_initialize_metric(
          ccc_lapse(index),
          ccc_betax(index), ccc_betay(index), ccc_betaz(index),
          ccc_gxx(index), ccc_gxy(index), ccc_gxz(index),
          ccc_gyy(index), ccc_gyz(index), ccc_gzz(index),
          &ADM_metric);

    ghl_ADM_aux_quantities metric_aux;
    ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

    // Read in primitive variables from gridfunctions
    ghl_primitive_quantities prims;
    ghl_initialize_primitives(
          rho(index), press(index), eps(index),
          vx(index), vy(index), vz(index),
          0.0, 0.0, 0.0,
          entropy[index], Y_e[index], temperature[index],
          &prims);

    // Read in conservative variables from gridfunctions
    ghl_conservative_quantities cons, cons_orig;
    ghl_initialize_conservatives(
          rho_star(index), tau(index),
          Stildex(index), Stildey(index), Stildez(index),
          ent_star[index], Ye_star[index], &cons);

    // Here we save the original values of conservative variables in cons_orig for debugging purposes.
    cons_orig = cons;

    //FIXME: might slow down the code. Was formerly a CCTK_WARN
    if(isnan(cons.rho*cons.tau*cons.SD[0]*cons.SD[1]*cons.SD[2])) {
      CCTK_VERROR("NaN found at start of C2P kernel:\n"
                  "  index %d %d %d, rho_* = %e, ~tau = %e, ~S_i = %e %e %e\n"
                  "  lapse = %e, shift = %e %e %e, gij = %e %e %e %e %e %e, Psi6 = %e\n",
                  i, j, k, cons.rho, cons.tau, cons.SD[0], cons.SD[1], cons.SD[2],
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
                      "rho_*, ~tau, ~S_{i}, ~DS, ~DY_e: %e, %e, %e, %e, %e, %e, %e\n"
                      "Undensitized conserved variables:\n"
                      "D, tau, S_{i}: %e %e %e %e %e\n"
                      "Output primitive variables:\n"
                      "rho, P: %e %e\n"
                      "v: %e %e %e\n"
                      "***********************************************************",
                      diagnostics.which_routine, ADM_metric.lapse, ADM_metric.betaU[0], ADM_metric.betaU[1], ADM_metric.betaU[2],
                      ADM_metric.gammaDD[0][0], ADM_metric.gammaDD[0][1], ADM_metric.gammaDD[0][2],
                      ADM_metric.gammaDD[1][1], ADM_metric.gammaDD[1][2], ADM_metric.gammaDD[2][2],
                      cons.rho, cons.tau, cons.SD[0], cons.SD[1], cons.SD[2], cons.entropy, cons.Y_e,
                      cons_undens.rho, cons_undens.tau, cons_undens.SD[0], cons_undens.SD[1], cons_undens.SD[2],
                      prims.rho, prims.press,
                      prims.vU[0], prims.vU[1], prims.vU[2]);
            }
      } else {
        //--------------------------------------------------
        //----------- Primitive recovery failed ------------
        //--------------------------------------------------
        local_failure_checker += 100;
        ghl_set_prims_to_constant_atm(ghl_eos, &prims);

        failures++;
        if(ADM_metric.sqrt_detgamma > ghl_params->psi6threshold) {
          failures_inhoriz++;
          pointcount_inhoriz++;
        }
        CCTK_VINFO("Con2Prim failed! Resetting to atmosphere...\n");
        CCTK_VINFO("rho_* = %e, ~tau = %e, ~S_i = %e %e %e\n"
                   "~DS = %e, ~DY_e = %e\n"
                   "lapse = %e, shift = %e %e %e, gij = %e %e %e %e %e %e, Psi6 = %e",
                   cons_orig.rho, cons_orig.tau, cons_orig.SD[0], cons_orig.SD[1], cons_orig.SD[2],
                   cons.entropy, cons.Y_e,
                   ADM_metric.lapse, ADM_metric.betaU[0], ADM_metric.betaU[1], ADM_metric.betaU[2],
                   ADM_metric.gammaDD[0][0], ADM_metric.gammaDD[0][1], ADM_metric.gammaDD[0][2],
                   ADM_metric.gammaDD[1][1], ADM_metric.gammaDD[1][2], ADM_metric.gammaDD[2][2], ADM_metric.sqrt_detgamma);
      }
    } else {
      local_failure_checker += 1;
      ghl_set_prims_to_constant_atm(ghl_eos, &prims);
      //rho_star_fix_applied++;
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

    double dummy1, dummy2, dummy3;
    ghl_return_primitives(
          &prims,
          &rho(index), &press(index), &eps(index),
          &vx(index), &vy(index), &vz(index),
          &dummy1, &dummy2, &dummy3,
          &entropy(index), &Y_e(index), &temperature(index));
    u0(index) = prims.u0;

    ghl_return_conservatives(
          &cons,
          &rho_star(index), &tau(index),
          &Stildex(index), &Stildey(index), &Stildez(index),
          &ent_star[index], &Ye_star[index]);

    /* The following code involves the diagnostics requiring reduction
    //Now we compute the difference between original & new conservatives, for diagnostic purposes:
    error_rho_numer += fabs(cons.rho - cons_orig.rho);
    error_tau_numer += fabs(cons.tau - cons_orig.tau);
    error_Sx_numer  += fabs(cons.SD[0] - cons_orig.SD[0]);
    error_Sy_numer  += fabs(cons.SD[1] - cons_orig.SD[1]);
    error_Sz_numer  += fabs(cons.SD[2] - cons_orig.SD[2]);
    error_entropy_numer  += fabs(cons.entropy - cons_orig.entropy);
    error_Ye_numer  += fabs(cons.Y_e - cons_orig.Y_e);
    error_rho_denom += cons_orig.tau;
    error_tau_denom += cons_orig.rho;
    error_Sx_denom  += fabs(cons_orig.SD[0]);
    error_Sy_denom  += fabs(cons_orig.SD[1]);
    error_Sz_denom  += fabs(cons_orig.SD[2]);
    error_entropy_denom  += cons_orig.entropy;
    error_Ye_denom  += cons_orig.Y_e;

    pointcount++;
    if(diagnostics.speed_limited) {
      local_failure_checker += 10;
      vel_limited_ptcount++;
    }
    backup0 += diagnostics.backup[0];
    backup1 += diagnostics.backup[1];
    backup2 += diagnostics.backup[2];
    n_iter += diagnostics.n_iter;
    failure_checker[index] = local_failure_checker
                           + 1000*diagnostics.backup[0]
                           + 10000*diagnostics.tau_fix
                           + 100000*diagnostics.Stilde_fix;
    */
  }); // ccc loop everywhere

  /*
    Failure checker decoder:
       1: atmosphere reset when rho_star < 0
      10: Limiting velocity u~ after C2P/Font Fix or v in ghl_enforce_primitive_limits_and_compute_u0
     100: Both C2P and Font Fix failed
      1k: backups used
     10k: tau~ was reset in ghl_apply_conservative_limits
    100k: S~ was reset in ghl_apply_conservative_limits
  */
  /* Again, diagnostics require reductions
  if(CCTK_Equals(verbose, "yes")) {
    if(ghl_eos->eos_type == ghl_eos_hybrid) {
      if(ghl_params->evolve_entropy) {
        CCTK_VINFO("C2P: Iter. # %d, Lev: %d NumPts= %d | Backups: %d %d %d | Fixes: VL= %d rho*= %d | Failures: %d InHoriz= %d / %d | %.2f iters/gridpt\n"
                   "   Error, Sum: rho %.3e, %.3e | tau %.3e, %.3e | entropy %.3e, %.3e\n"
                   "               Sx %.3e, %.3e | Sy %.3e, %.3e | Sz %.3e, %.3e\n",
                   cctk_iteration, (int)GetRefinementLevel(cctkGH), pointcount,
                   backup0, backup1, backup2,
                   vel_limited_ptcount, rho_star_fix_applied,
                   failures, failures_inhoriz, pointcount_inhoriz,
                   (double)n_iter/( (double)(cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]) ),
                   rho_error, error_rho_denom,
                   tau_error, error_tau_denom,
                   entropy_error, error_entropy_denom,
                   Sx_error, error_Sx_denom,
                   Sy_error, error_Sy_denom,
                   Sz_error, error_Sz_denom);
      } else {
        CCTK_VINFO("C2P: Iter. # %d, Lev: %d NumPts= %d | Backups: %d %d %d | Fixes: VL= %d rho*= %d | Failures: %d InHoriz= %d / %d | %.2f iters/gridpt\n"
                   "   Error, Sum: rho %.3e, %.3e | tau %.3e, %.3e | Sx %.3e, %.3e | Sy %.3e, %.3e | Sz %.3e, %.3e\n",
                   cctk_iteration, (int)GetRefinementLevel(cctkGH), pointcount,
                   backup0, backup1, backup2,
                   vel_limited_ptcount, rho_star_fix_applied,
                   failures, failures_inhoriz, pointcount_inhoriz,
                   (double)n_iter/( (double)(cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]) ),
                   rho_error, error_rho_denom,
                   tau_error, error_tau_denom,
                   Sx_error, error_Sx_denom,
                   Sy_error, error_Sy_denom,
                   Sz_error, error_Sz_denom);
      }
    } else if( ghl_eos->eos_type == ghl_eos_tabulated ) {
      if(ghl_params->evolve_entropy) {
        CCTK_VINFO("C2P: Iter. # %d, Lev: %d NumPts= %d | Backups: %d %d %d | Fixes: VL= %d rho*= %d | Failures: %d InHoriz= %d / %d | %.2f iters/gridpt\n"
                   "   Error, Sum: rho %.3e, %.3e | tau %.3e, %.3e | entropy %.3e, %.3e | Y_e %.3e, %.3e\n"
                   "               Sx %.3e, %.3e | Sy %.3e, %.3e | Sz %.3e, %.3e\n",
                   cctk_iteration, (int)GetRefinementLevel(cctkGH), pointcount,
                   backup0, backup1, backup2,
                   vel_limited_ptcount, rho_star_fix_applied,
                   failures, failures_inhoriz, pointcount_inhoriz,
                   (double)n_iter/( (double)(cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]) ),
                   rho_error, error_rho_denom,
                   tau_error, error_tau_denom,
                   entropy_error, error_entropy_denom,
                   Ye_error, error_Ye_denom,
                   Sx_error, error_Sx_denom,
                   Sy_error, error_Sy_denom,
                   Sz_error, error_Sz_denom);
      } else {
        CCTK_VINFO("C2P: Iter. # %d, Lev: %d NumPts= %d | Backups: %d %d %d | Fixes: VL= %d rho*= %d | Failures: %d InHoriz= %d / %d | %.2f iters/gridpt\n"
                   "   Error, Sum: rho %.3e, %.3e | tau %.3e, %.3e | Y_e %.3e, %.3e\n"
                   "               Sx %.3e, %.3e | Sy %.3e, %.3e | Sz %.3e, %.3e\n",
                   cctk_iteration, (int)GetRefinementLevel(cctkGH), pointcount,
                   backup0, backup1, backup2,
                   vel_limited_ptcount, rho_star_fix_applied,
                   failures, failures_inhoriz, pointcount_inhoriz,
                   (double)n_iter/( (double)(cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]) ),
                   rho_error, error_rho_denom,
                   tau_error, error_tau_denom,
                   Ye_error, error_Ye_denom,
                   Sx_error, error_Sx_denom,
                   Sy_error, error_Sy_denom,
                   Sz_error, error_Sz_denom);
      }
    }
  }
  */
}
