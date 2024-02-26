#include "GRHayLHD.h"
#include "Symmetry.h"

void GRHayLHD_hybrid_conservs_to_prims(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLHD_hybrid_conservs_to_prims;
  DECLARE_CCTK_PARAMETERS;

  const int imax = cctk_lsh[0];
  const int jmax = cctk_lsh[1];
  const int kmax = cctk_lsh[2];

  if(CCTK_EQUALS(Symmetry, "equatorial")) {
    // SET SYMMETRY GHOSTZONES ON ALL CONSERVATIVE VARIABLES!
    int ierr = 0;
    ierr += CartSymGN(cctkGH, "GRHayLHD::grmhd_velocities");
    // FIXME: UGLY. Filling metric ghostzones is needed for, e.g., Cowling runs.
    ierr += CartSymGN(cctkGH, "lapse::lapse_vars");
    ierr += CartSymGN(cctkGH, "bssn::BSSN_vars");
    ierr += CartSymGN(cctkGH, "bssn::BSSN_AH");
    ierr += CartSymGN(cctkGH, "shift::shift_vars");
    if(ierr != 0)
      CCTK_VERROR("Error with setting equatorial symmetries in con2prim.");
  }

  // Diagnostic variables.
  int failures = 0;
  int vel_limited_ptcount = 0;
  int rho_star_fix_applied = 0;
  int navg_total = 0;
  int failures_inhoriz = 0;
  int pointcount_inhoriz = 0;
  int backup0 = 0;
  int backup1 = 0;
  int backup2 = 0;
  int n_iter = 0;
  int pointcount_avg = 0;
  int pointcount_Font = 0;

#pragma omp parallel for reduction(+:                              \
      backup0, backup1, backup2, vel_limited_ptcount,              \
      rho_star_fix_applied, failures, failures_inhoriz,            \
      pointcount_inhoriz, n_iter, pointcount_avg, pointcount_Font) \
      schedule(static)
  for(int k=0; k<kmax; k++) {
    for(int j=0; j<jmax; j++) {
      for(int i=0; i<imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);

        CCTK_REAL local_failure_checker = 0;

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
        // The code has only ever been tested using the default GRHayL guess,
        // so using the previous timelevel as an initial guess would need to
        // be implemented here.
        ghl_primitive_quantities prims;
        prims.BU[0] = prims.BU[1] = prims.BU[2] = 0.0;

        // Read in conservative variables from gridfunctions
        ghl_conservative_quantities cons, cons_undens;
        cons.rho   = rho_star[index];
        cons.tau   = tau[index];
        cons.SD[0] = Stildex[index];
        cons.SD[1] = Stildey[index];
        cons.SD[2] = Stildez[index];


        int check;

        /************* Main conservative-to-primitive logic ************/
        if (cons.rho > 0.0 &&
            isfinite(cons.rho*cons.tau*cons.SD[0]*cons.SD[1]*cons.SD[2])) {

          ghl_apply_conservative_limits(
              ghl_params, ghl_eos, &ADM_metric,
              &prims, &cons, &diagnostics);

          ghl_undensitize_conservatives(ADM_metric.sqrt_detgamma, &cons, &cons_undens);

          /************* Conservative-to-primitive recovery ************/
          check = ghl_con2prim_multi_method(
                ghl_params, ghl_eos, &ADM_metric, &metric_aux,
                &cons_undens, &prims, &diagnostics);
        } else {
          local_failure_checker += 1;
          check = 1;
        }

        if(check || isnan(prims.rho*prims.press*prims.eps*prims.vU[0]*prims.vU[1]*prims.vU[2])) {
          pointcount_avg++;
          check = 1;

          // Reload cons for consistency (hybrid applies limits to cons)
          cons.rho     = rho_star[index];
          cons.tau     = tau[index];
          cons.SD[0]   = Stildex[index];
          cons.SD[1]   = Stildey[index];
          cons.SD[2]   = Stildez[index];

          ghl_conservative_quantities cons_neigh_avg, cons_avg;
          cons_neigh_avg.rho     = 0.0;
          cons_neigh_avg.tau     = 0.0;
          cons_neigh_avg.SD[0]   = 0.0;
          cons_neigh_avg.SD[1]   = 0.0;
          cons_neigh_avg.SD[2]   = 0.0;

          const int iavg_min = MAX(0, i-1);
          const int javg_min = MAX(0, j-1);
          const int kavg_min = MAX(0, k-1);
          const int iavg_max = MIN(imax, i+2);
          const int javg_max = MIN(jmax, j+2);
          const int kavg_max = MIN(kmax, k+2);

          int n_avg = 0;
          // We compute the average of neighboring points once and
          // reuse it for the various weighted averages.
          for(int kavg=kavg_min; kavg<kavg_max; kavg++) {
            for(int javg=javg_min; javg<javg_max; javg++) {
              for(int iavg=iavg_min; iavg<iavg_max; iavg++) {
                // Skip this point
                const int inavg = CCTK_GFINDEX3D(cctkGH,iavg,javg,kavg);
                if((index==inavg))
                  continue;
                cons_neigh_avg.rho   += rho_star[inavg];
                cons_neigh_avg.tau   += tau[inavg];
                cons_neigh_avg.SD[0] += Stildex[inavg];
                cons_neigh_avg.SD[1] += Stildey[inavg];
                cons_neigh_avg.SD[2] += Stildez[inavg];
                n_avg++;
              }
            }
          }

          int avg_weight = 1;
          while(check && avg_weight < 5) {
            check = 0;
            // last point doesn't add central point and has 1 less point
            // being averaged.
            n_avg += (avg_weight!=4);

            const CCTK_REAL wfac = avg_weight/4.0;
            const CCTK_REAL cfac = 1.0 - wfac;
            cons_avg.rho   = wfac*cons_neigh_avg.rho     + cfac*cons.rho;
            cons_avg.tau   = wfac*cons_neigh_avg.tau     + cfac*cons.tau;
            cons_avg.SD[0] = wfac*cons_neigh_avg.SD[0]   + cfac*cons.SD[0];
            cons_avg.SD[1] = wfac*cons_neigh_avg.SD[1]   + cfac*cons.SD[1];
            cons_avg.SD[2] = wfac*cons_neigh_avg.SD[2]   + cfac*cons.SD[2];

            cons_avg.rho   /= n_avg;
            cons_avg.tau   /= n_avg;
            cons_avg.SD[0] /= n_avg;
            cons_avg.SD[1] /= n_avg;
            cons_avg.SD[2] /= n_avg;

            ghl_undensitize_conservatives(ADM_metric.sqrt_detgamma, &cons_avg, &cons_undens);

            /************* Conservative-to-primitive recovery ************/
            check = ghl_con2prim_multi_method(
                  ghl_params, ghl_eos, &ADM_metric, &metric_aux,
                  &cons_undens, &prims, &diagnostics);

            avg_weight++;
            if(isnan(prims.rho*prims.press*prims.eps*prims.vU[0]*prims.vU[1]*prims.vU[2]) )
              check = 1;
          }
          if(check) {
            // We are still failing after exhausting the averaging options.
            // Next, we try Font1D.
            pointcount_Font++;

            ghl_apply_conservative_limits(
                ghl_params, ghl_eos, &ADM_metric,
                &prims, &cons, &diagnostics);

            ghl_undensitize_conservatives(ADM_metric.sqrt_detgamma, &cons, &cons_undens);

            check = ghl_hybrid_Font1D(
                  ghl_params, ghl_eos, &ADM_metric, &metric_aux,
                  &cons_undens, &prims, &diagnostics);

            if(isnan(prims.rho*prims.press*prims.eps*prims.vU[0]*prims.vU[1]*prims.vU[2]*prims.entropy*prims.Y_e*prims.temperature) )
              check = 1;

            if(check) {
              // We are still failing after exhausting the averaging options.
              // We'll surrender and resort to atmospheric reset...

              failure_checker[index] += 100;

              ghl_set_prims_to_constant_atm(ghl_eos, &prims);

              failures++;
              if(ADM_metric.sqrt_detgamma > ghl_params->psi6threshold) {
                failures_inhoriz++;
                pointcount_inhoriz++;
              }
              CCTK_VINFO("***********************************************************\n"
                         "Con2Prim and averaging backups failed! Resetting to atmosphere...\n"
                         "position = %e %e %e\n"
                         "lapse, shift = %e, %e, %e, %e\n"
                         "gij = %e, %e, %e, %e, %e, %e\n"
                         "rho_*, ~tau, ~S_{i}: %e, %e, %e, %e, %e\n"
                       "***********************************************************",
                         x[index], y[index], z[index],
                         ADM_metric.lapse, ADM_metric.betaU[0], ADM_metric.betaU[1], ADM_metric.betaU[2],
                         ADM_metric.gammaDD[0][0], ADM_metric.gammaDD[0][1], ADM_metric.gammaDD[0][2],
                         ADM_metric.gammaDD[1][1], ADM_metric.gammaDD[1][2], ADM_metric.gammaDD[2][2],
                         cons.rho, cons.tau, cons.SD[0], cons.SD[1], cons.SD[2]);
            } // atmospheric backup
          } // Font1D backup
        } // if c2p failed
        /***************************************************************/

        //--------------------------------------------------
        //---------- Primitive recovery succeeded ----------
        //--------------------------------------------------
        // Enforce limits on primitive variables and recompute conservatives.
        diagnostics.speed_limited += ghl_enforce_primitive_limits_and_compute_u0(
              ghl_params, ghl_eos, &ADM_metric, &prims);

        rho[index]   = prims.rho;
        press[index] = prims.press;
        eps[index]   = prims.eps;
        u0[index]    = prims.u0;
        vx[index]    = prims.vU[0];
        vy[index]    = prims.vU[1];
        vz[index]    = prims.vU[2];

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
      }
    }
  }

  CCTK_REAL error_rho_numer = 0;
  CCTK_REAL error_tau_numer = 0;
  CCTK_REAL error_Sx_numer = 0;
  CCTK_REAL error_Sy_numer = 0;
  CCTK_REAL error_Sz_numer = 0;
  CCTK_REAL error_Ye_numer = 0;

  CCTK_REAL error_rho_denom = 0;
  CCTK_REAL error_tau_denom = 0;
  CCTK_REAL error_Sx_denom = 0;
  CCTK_REAL error_Sy_denom = 0;
  CCTK_REAL error_Sz_denom = 0;

// We split the prim2con into a separate loop so the averaging
// loop is deterministic and doesn't have a race condition
#pragma omp parallel for reduction(+:                                   \
      error_rho_numer, error_tau_numer, error_Sx_numer, error_Sy_numer, \
      error_Sz_numer, error_rho_denom, error_tau_denom, error_Sx_denom, \
      error_Sy_denom, error_Sz_denom) \
      schedule(static)
  for(int k=0; k<kmax; k++) {
    for(int j=0; j<jmax; j++) {
      for(int i=0; i<imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);

        ghl_metric_quantities ADM_metric;
        ghl_enforce_detgtij_and_initialize_ADM_metric(
              alp[index],
              betax[index], betay[index], betaz[index],
              gxx[index], gxy[index], gxz[index],
              gyy[index], gyz[index], gzz[index],
              &ADM_metric);

        ghl_ADM_aux_quantities metric_aux;
        ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

        ghl_primitive_quantities prims;
        prims.BU[0] = prims.BU[1] = prims.BU[2] = 0.0;
        prims.rho     = rho[index];
        prims.press   = press[index];
        prims.eps     = eps[index];
        prims.u0      = u0[index];
        prims.vU[0]   = vx[index];
        prims.vU[1]   = vy[index];
        prims.vU[2]   = vz[index];

        ghl_conservative_quantities cons, cons_orig;
        cons_orig.rho     = rho_star[index];
        cons_orig.tau     = tau[index];
        cons_orig.SD[0]   = Stildex[index];
        cons_orig.SD[1]   = Stildey[index];
        cons_orig.SD[2]   = Stildez[index];

        ghl_compute_conservs(&ADM_metric, &metric_aux, &prims, &cons);

        rho_star[index] = cons.rho;
        tau[index]      = cons.tau;
        Stildex[index]  = cons.SD[0];
        Stildey[index]  = cons.SD[1];
        Stildez[index]  = cons.SD[2];

        // Now we compute the difference between original & new conservatives
        // for diagnostic purposes:
        error_rho_numer += fabs(cons.rho - cons_orig.rho);
        error_tau_numer += fabs(cons.tau - cons_orig.tau);
        error_Sx_numer  += fabs(cons.SD[0] - cons_orig.SD[0]);
        error_Sy_numer  += fabs(cons.SD[1] - cons_orig.SD[1]);
        error_Sz_numer  += fabs(cons.SD[2] - cons_orig.SD[2]);
        error_rho_denom += cons_orig.tau;
        error_tau_denom += cons_orig.rho;
        error_Sx_denom  += fabs(cons_orig.SD[0]);
        error_Sy_denom  += fabs(cons_orig.SD[1]);
        error_Sz_denom  += fabs(cons_orig.SD[2]);
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
  if(CCTK_Equals(verbose, "yes")) {
    const int pointcount = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
    const CCTK_REAL rho_error = (error_rho_denom == 0) ? error_rho_numer : error_rho_numer/error_rho_denom;
    const CCTK_REAL tau_error = (error_tau_denom == 0) ? error_tau_numer : error_tau_numer/error_tau_denom;
    const CCTK_REAL Sx_error  = (error_Sx_denom == 0)  ? error_Sx_numer  : error_Sx_numer/error_Sx_denom;
    const CCTK_REAL Sy_error  = (error_Sy_denom == 0)  ? error_Sy_numer  : error_Sy_numer/error_Sy_denom;
    const CCTK_REAL Sz_error  = (error_Sz_denom == 0)  ? error_Sz_numer  : error_Sz_numer/error_Sz_denom;

    CCTK_VINFO(
        "C2P: Iter. # %d, Lev: %d NumPts= %d | Backups: %d %d %d | Fixes: VL= %d rho*= %d\n"
        "                 Averaged pts = %d Font1D %d | Failures: %d InHoriz= %d / %d | %.2f iters/gridpt\n"
        "   Error, Sum: rho %.3e, %.3e | tau %.3e, %.3e | "
        "Sx %.3e, %.3e | Sy %.3e, %.3e | Sz %.3e, %.3e\n",
        cctk_iteration, (int)GetRefinementLevel(cctkGH), pointcount, backup0,
        backup1, backup2, vel_limited_ptcount, rho_star_fix_applied,
        pointcount_avg, pointcount_Font,
        failures, failures_inhoriz, pointcount_inhoriz,
        (CCTK_REAL)n_iter / ((CCTK_REAL)(pointcount)),
        rho_error, error_rho_denom, tau_error, error_tau_denom,
        Sx_error, error_Sx_denom, Sy_error, error_Sy_denom,
        Sz_error, error_Sz_denom);
  }
}
