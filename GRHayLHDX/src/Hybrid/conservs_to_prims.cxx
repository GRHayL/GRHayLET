#include "GRHayLHDX.h"

extern "C" void GRHayLHDX_hybrid_conservs_to_prims(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_hybrid_conservs_to_prims;
  DECLARE_CCTK_PARAMETERS;

  const int imax = cctk_lsh[0];
  const int jmax = cctk_lsh[1];
  const int kmax = cctk_lsh[2];

  // Diagnostic variables.
  //int failures = 0;
  //int vel_limited_ptcount = 0;
  //int rho_star_fix_applied = 0;
  //int failures_inhoriz = 0;
  //int pointcount_inhoriz = 0;
  //int backup0 = 0;
  //int backup1 = 0;
  //int backup2 = 0;
  //int n_iter = 0;
  //int pointcount_avg = 0;
  //int pointcount_Font = 0;

  const Loop::GF3D2layout layout(cctkGH, {1, 1, 1});

//I don't think I can do this sort of diagnostic summing in carpetx
//#pragma omp parallel for reduction(+: \
//      pointcount, backup0, backup1, backup2, vel_limited_ptcount, rho_star_fix_applied, failures, failures_inhoriz, pointcount_inhoriz, n_iter, \
//      error_rho_numer, error_tau_numer, error_Sx_numer, error_Sy_numer, error_Sz_numer, error_entropy_numer, error_Ye_numer,                    \
//      error_rho_denom, error_tau_denom, error_Sx_denom, error_Sy_denom, error_Sz_denom, error_entropy_denom, error_Ye_denom) schedule(static)
  grid.loop_all<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index index(layout, p.I);

    CCTK_REAL local_failure_checker = 0;

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
    prims.BU[0] = prims.BU[1] = prims.BU[2] = 0.0;

    // Read in conservative variables from gridfunctions
    ghl_conservative_quantities cons, cons_undens;
    cons.rho   = rho_star(index); 
    cons.tau   = tau(index);
    cons.SD[0] = Stildex(index);
    cons.SD[1] = Stildey(index);
    cons.SD[2] = Stildez(index);

    int check;

    /************* Main conservative-to-primitive logic ************/
    if(cons.rho>0.0) {
      // Apply the tau floor
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
      //pointcount_avg++;
      check = 1;

      // Reload cons for consistency (hybrid applies limits to cons)
      cons.rho   = rho_star(index);
      cons.tau   = tau(index);
      cons.SD[0] = Stildex(index);
      cons.SD[1] = Stildey(index);
      cons.SD[2] = Stildez(index);

      ghl_conservative_quantities cons_neigh_avg, cons_avg;
      cons_neigh_avg.rho   = 0.0;
      cons_neigh_avg.tau   = 0.0;
      cons_neigh_avg.SD[0] = 0.0;
      cons_neigh_avg.SD[1] = 0.0;
      cons_neigh_avg.SD[2] = 0.0;

      // Set minimum index of loop to -1 unless that index is below the minimum index of this tile
      const int iavg_min = (0 > p.i-1) - 1;
      const int javg_min = (0 > p.j-1) - 1;
      const int kavg_min = (0 > p.k-1) - 1;
      // Set maximum index of loop to +1 unless that index is above the maximum index of this tile
      const int iavg_max = -(imax < p.i+2) + 1;
      const int javg_max = -(jmax < p.j+2) + 1;
      const int kavg_max = -(kmax < p.k+2) + 1;

      int n_avg = 0;
      // We compute the average of neighboring points once and
      // reuse it for the various weighted averages.
      for(int kavg=kavg_min; kavg<kavg_max; kavg++) {
        for(int javg=javg_min; javg<javg_max; javg++) {
          for(int iavg=iavg_min; iavg<iavg_max; iavg++) {
            // Skip this point
            const Loop::GF3D2index inavg(layout, p.I + iavg*p.DI[0] + javg*p.DI[1] + kavg*p.DI[2]);
            if((iavg==p.i && javg==p.j && kavg==p.k))
              continue;
            cons_neigh_avg.rho   += rho_star(inavg);
            cons_neigh_avg.tau   += tau(inavg);
            cons_neigh_avg.SD[0] += Stildex(inavg);
            cons_neigh_avg.SD[1] += Stildey(inavg);
            cons_neigh_avg.SD[2] += Stildez(inavg);
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
        //pointcount_Font++;

        ghl_apply_conservative_limits(
            ghl_params, ghl_eos, &ADM_metric,
            &prims, &cons, &diagnostics);

        ghl_undensitize_conservatives(ADM_metric.sqrt_detgamma, &cons, &cons_undens);

        check = ghl_hybrid_Font1D(
              ghl_params, ghl_eos, &ADM_metric, &metric_aux,
              &cons_undens, &prims, &diagnostics);

        if(isnan(prims.rho*prims.press*prims.eps*prims.vU[0]*prims.vU[1]*prims.vU[2]) )
          check = 1;

        if(check) {
          // We are still failing after exhausting the averaging options.
          // We'll surrender and resort to atmospheric reset...

          local_failure_checker += 100;

          ghl_set_prims_to_constant_atm(ghl_eos, &prims);

          //failures++;
          if(ADM_metric.sqrt_detgamma > ghl_params->psi6threshold) {
            //failures_inhoriz++;
            //pointcount_inhoriz++;
          }
          CCTK_VINFO("***********************************************************\n"
                     "Con2Prim and averaging backups failed! Resetting to atmosphere...\n"
                     "position = %e %e %e\n"
                     "lapse, shift = %e, %e, %e, %e\n"
                     "gij = %e, %e, %e, %e, %e, %e\n"
                     "rho_*, ~tau, ~S_{i}: %e, %e, %e, %e, %e\n"
                     "***********************************************************",
                     p.x, p.y, p.z,
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
    ghl_compute_conservs(
          &ADM_metric, &metric_aux, &prims, &cons);

    rho(index)   = prims.rho;
    press(index) = prims.press;
    eps(index)   = prims.eps;
    u0(index)    = prims.u0;
    vx(index)    = prims.vU[0];
    vy(index)    = prims.vU[1];
    vz(index)    = prims.vU[2];

    rho_star(index) = cons.rho;
    tau(index)      = cons.tau;
    Stildex(index)  = cons.SD[0];
    Stildey(index)  = cons.SD[1];
    Stildez(index)  = cons.SD[2];

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

    if(diagnostics.speed_limited) {
      local_failure_checker += 10;
      vel_limited_ptcount++;
    }
    backup0 += diagnostics.backup[0];
    backup1 += diagnostics.backup[1];
    backup2 += diagnostics.backup[2];
    n_iter += diagnostics.n_iter;
    */
    failure_checker(index) = local_failure_checker
                           + 1000*diagnostics.backup[0]
                           + 10000*diagnostics.tau_fix
                           + 100000*diagnostics.Stilde_fix;
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
  */
}
