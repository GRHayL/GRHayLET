#include "IllinoisGRMHD.h"
#include "Symmetry.h"

static inline bool check_depsdT_condition(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_conservative_quantities *restrict cons_undens,
      const ghl_primitive_quantities *restrict prims_guess,
      const CCTK_REAL x) {

  double SD[3] = {cons_undens->SD[0], cons_undens->SD[1], cons_undens->SD[2]};
  double S_squared = ghl_compute_vec2_from_vec3D(ADM_metric->gammaUU, SD);
  const double tau = MAX(cons_undens->tau, 0.99*eos->tau_atm);

  // Step 2: Enforce ceiling on S^{2} (Eq. A5 of [1])
  // Step 2.1: Compute maximum allowed value for S^{2}
  const double S_squared_max = SQR(tau + cons_undens->rho);
  if(S_squared > S_squared_max) {
    // Step 2.2: Rescale S_{i}
    const double rescale_factor = sqrt(0.9999*S_squared_max/S_squared);
    for(int i=0;i<3;i++)
      SD[i] *= rescale_factor;

    // Step 2.3: Recompute S^{2}
    S_squared = ghl_compute_vec2_from_vec3D(ADM_metric->gammaUU, SD);
  }

  // Step 3: Compute B^{2} = gamma_{ij}B^{i}B^{j}
  const double B_squared = ghl_compute_vec2_from_vec3D(ADM_metric->gammaDD, prims_guess->BU);

  // Step 4: Compute B.S = B^{i}S_{i}
  double BdotS = 0.0;
  for(int i=0;i<3;i++) BdotS += prims_guess->BU[i]*SD[i];

  const double q = cons_undens->tau/cons_undens->rho;
  const double r = S_squared/(cons_undens->rho*cons_undens->rho);
  const double s = B_squared/cons_undens->rho;
  const double t = BdotS/(pow(cons_undens->rho, 1.5));
  double Wminus2 = 1.0 - ( x*x*r + (2*x+s)*t*t ) / ( x*x*(x+s)*(x+s) );
  Wminus2        = fmin(fmax(Wminus2, params->inv_sq_max_Lorentz_factor), 1.0);
  double W       = pow(Wminus2, -0.5);

  double rho     = cons_undens->rho/W;
  double ye      = cons_undens->Y_e/cons_undens->rho;
  double eps     = W - 1.0 + (1.0-W*W)*x/W + W*(q - s + t*t/(2*x*x) + s/(2*W*W)  );
  double press, ent, depsdT, temp;

  ghl_tabulated_enforce_bounds_rho_Ye_eps(ghl_eos, &rho, &ye, &eps);
  ghl_tabulated_compute_P_S_depsdT_T_from_eps(ghl_eos, rho, ye, eps, &press, &ent, &depsdT, &temp);

  if(depsdT < 5.0000000000000002396e-05) {
    return true;
  } else {
    return false;
  }
}

void IllinoisGRMHD_tabulated_entropy_conservs_to_prims(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_IllinoisGRMHD_tabulated_entropy_conservs_to_prims;
  DECLARE_CCTK_PARAMETERS;

  const int imax = cctk_lsh[0];
  const int jmax = cctk_lsh[1];
  const int kmax = cctk_lsh[2];

  if(CCTK_EQUALS(Symmetry, "equatorial")) {
    // SET SYMMETRY GHOSTZONES ON ALL CONSERVATIVE VARIABLES!
    int ierr = 0;
    ierr += CartSymGN(cctkGH, "IllinoisGRMHD::grmhd_conservatives");
    ierr += CartSymGN(cctkGH, "IllinoisGRMHD::Ye_star");
    ierr += CartSymGN(cctkGH, "IllinoisGRMHD::ent_star");
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
  int failures_inhoriz = 0;
  int pointcount_inhoriz = 0;
  int backup0 = 0;
  int backup1 = 0;
  int backup2 = 0;
  int n_iter = 0;
  int pointcount_avg = 0;

#pragma omp parallel for reduction(+:                   \
      backup0, backup1, backup2, vel_limited_ptcount,   \
      rho_star_fix_applied, failures, failures_inhoriz, \
      pointcount_inhoriz, n_iter, pointcount_avg)       \
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
        prims.BU[0] = Bx_center[index];
        prims.BU[1] = By_center[index];
        prims.BU[2] = Bz_center[index];

        // Read in conservative variables from gridfunctions
        ghl_conservative_quantities cons, cons_undens;
        cons.rho     = rho_star[index];
        cons.tau     = tau[index];
        cons.SD[0]   = Stildex[index];
        cons.SD[1]   = Stildey[index];
        cons.SD[2]   = Stildez[index];
        cons.entropy = ent_star[index];
        cons.Y_e     = Ye_star[index];

        int check;

        /************* Main conservative-to-primitive logic ************/
        if(cons.rho>0.0) {
          ghl_undensitize_conservatives(ADM_metric.sqrt_detgamma, &cons, &cons_undens);

          ghl_guess_primitives(ghl_eos, &ADM_metric, &cons, &prims);

          const double invD = 1.0/cons_undens.rho;
          const double Bsq  = ghl_compute_vec2_from_vec3D(ADM_metric.gammaDD, prims.BU);
          const double tmp1 = cons_undens.tau + invD;
          const double tmp2 = Bsq*invD;
          const double xlow = 1.0 + tmp1 - tmp2;
          const double xup  = 2.0*(1 + tmp1) - tmp2;
          bool test1 = check_depsdT_condition(ghl_params, ghl_eos, &ADM_metric, &cons_undens, &prims, xlow);
          bool test2 = check_depsdT_condition(ghl_params, ghl_eos, &ADM_metric, &cons_undens, &prims, xup);
          const bool use_entropy = (test1 || test2);

          // Store primitive guesses (used if con2prim fails)
          ghl_primitive_quantities prims_guess = prims;
          const double Tguess[3] = {ghl_eos->T_atm, ghl_eos->T_min, ghl_eos->T_max};
          for(int iter=0; iter<3; iter++) {
            prims_guess.temperature = Tguess[iter];

            if(use_entropy) {
              prims = prims_guess;
              check = ghl_tabulated_Palenzuela1D_entropy(
                    ghl_params, ghl_eos, &ADM_metric, &metric_aux,
                    &cons_undens, &prims, &diagnostics);
              if(check) {
                prims = prims_guess;
                check = ghl_tabulated_Palenzuela1D_energy(
                      ghl_params, ghl_eos, &ADM_metric, &metric_aux,
                      &cons_undens, &prims, &diagnostics);
                diagnostics.backup[0] = true;
                if(check) {
                  diagnostics.backup[1] = true;
                  prims = prims_guess;
                  check = ghl_tabulated_Newman1D_entropy(
                        ghl_params, ghl_eos, &ADM_metric, &metric_aux,
                        &cons_undens, &prims, &diagnostics);
                  if(check) {
                    diagnostics.backup[2] = true;
                    prims = prims_guess;
                    check = ghl_tabulated_Newman1D_energy(
                          ghl_params, ghl_eos, &ADM_metric, &metric_aux,
                          &cons_undens, &prims, &diagnostics);
                  }
                }
              }
            } else {
              prims = prims_guess;
              check = ghl_tabulated_Palenzuela1D_energy(
                    ghl_params, ghl_eos, &ADM_metric, &metric_aux,
                    &cons_undens, &prims, &diagnostics);
              if(check) {
                diagnostics.backup[2] = true;
                prims = prims_guess;
                check = ghl_tabulated_Newman1D_energy(
                      ghl_params, ghl_eos, &ADM_metric, &metric_aux,
                      &cons_undens, &prims, &diagnostics);
              }
            }
          }
        } else {
          ghl_set_prims_to_constant_atm(ghl_eos, &prims);
          local_failure_checker += 1;
          rho_star_fix_applied++;
          check = 0;
        }

        if(check || isnan(prims.rho*prims.press*prims.eps*prims.vU[0]*prims.vU[1]*prims.vU[2]*
                          prims.entropy*prims.Y_e*prims.temperature)) {
          pointcount_avg++;
          check = 1;

          ghl_conservative_quantities cons_neigh_avg = { 0 };

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
                cons_neigh_avg.rho     += rho_star[inavg];
                cons_neigh_avg.tau     += tau[inavg];
                cons_neigh_avg.SD[0]   += Stildex[inavg];
                cons_neigh_avg.SD[1]   += Stildey[inavg];
                cons_neigh_avg.SD[2]   += Stildez[inavg];
                cons_neigh_avg.entropy += ent_star[inavg];
                cons_neigh_avg.Y_e     += Ye_star[inavg];
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

            ghl_conservative_quantities cons_avg;
            cons_avg.rho     = wfac*cons_neigh_avg.rho     + cfac*cons.rho;
            cons_avg.tau     = wfac*cons_neigh_avg.tau     + cfac*cons.tau;
            cons_avg.SD[0]   = wfac*cons_neigh_avg.SD[0]   + cfac*cons.SD[0];
            cons_avg.SD[1]   = wfac*cons_neigh_avg.SD[1]   + cfac*cons.SD[1];
            cons_avg.SD[2]   = wfac*cons_neigh_avg.SD[2]   + cfac*cons.SD[2];
            cons_avg.entropy = wfac*cons_neigh_avg.entropy + cfac*cons.entropy;
            cons_avg.Y_e     = wfac*cons_neigh_avg.Y_e     + cfac*cons.Y_e;

            cons_avg.rho     /= n_avg;
            cons_avg.tau     /= n_avg;
            cons_avg.SD[0]   /= n_avg;
            cons_avg.SD[1]   /= n_avg;
            cons_avg.SD[2]   /= n_avg;
            cons_avg.entropy /= n_avg;
            cons_avg.Y_e     /= n_avg;

            ghl_undensitize_conservatives(ADM_metric.sqrt_detgamma, &cons, &cons_undens);

            ghl_guess_primitives(ghl_eos, &ADM_metric, &cons, &prims);
  
            const double invD = 1.0/cons_undens.rho;
            const double Bsq  = ghl_compute_vec2_from_vec3D(ADM_metric.gammaDD, prims.BU);
            const double tmp1 = cons_undens.tau + invD;
            const double tmp2 = Bsq*invD;
            const double xlow = 1.0 + tmp1 - tmp2;
            const double xup  = 2.0*(1 + tmp1) - tmp2;
            bool test1 = check_depsdT_condition(ghl_params, ghl_eos, &ADM_metric, &cons_undens, &prims, xlow);
            bool test2 = check_depsdT_condition(ghl_params, ghl_eos, &ADM_metric, &cons_undens, &prims, xup);
            const bool use_entropy = (test1 || test2);
  
            // Store primitive guesses (used if con2prim fails)
            ghl_primitive_quantities prims_guess = prims;
            const double Tguess[3] = {ghl_eos->T_atm, ghl_eos->T_min, ghl_eos->T_max};
            for(int iter=0; iter<3; iter++) {
              prims_guess.temperature = Tguess[iter];
  
              if(use_entropy) {
                prims = prims_guess;
                check = ghl_tabulated_Palenzuela1D_entropy(
                      ghl_params, ghl_eos, &ADM_metric, &metric_aux,
                      &cons_undens, &prims, &diagnostics);
                if(check) {
                  prims = prims_guess;
                  check = ghl_tabulated_Palenzuela1D_energy(
                        ghl_params, ghl_eos, &ADM_metric, &metric_aux,
                        &cons_undens, &prims, &diagnostics);
                  diagnostics.backup[0] = true;
                  if(check) {
                    diagnostics.backup[1] = true;
                    prims = prims_guess;
                    check = ghl_tabulated_Newman1D_entropy(
                          ghl_params, ghl_eos, &ADM_metric, &metric_aux,
                          &cons_undens, &prims, &diagnostics);
                    if(check) {
                      diagnostics.backup[2] = true;
                      prims = prims_guess;
                      check = ghl_tabulated_Newman1D_energy(
                            ghl_params, ghl_eos, &ADM_metric, &metric_aux,
                            &cons_undens, &prims, &diagnostics);
                    }
                  }
                }
              } else {
                prims = prims_guess;
                check = ghl_tabulated_Palenzuela1D_energy(
                      ghl_params, ghl_eos, &ADM_metric, &metric_aux,
                      &cons_undens, &prims, &diagnostics);
                if(check) {
                  diagnostics.backup[2] = true;
                  prims = prims_guess;
                  check = ghl_tabulated_Newman1D_energy(
                        ghl_params, ghl_eos, &ADM_metric, &metric_aux,
                        &cons_undens, &prims, &diagnostics);
                }
              }
            }

            avg_weight++;
            if(isnan(prims.rho*prims.press*prims.eps*prims.vU[0]*prims.vU[1]*prims.vU[2]*
                     prims.entropy*prims.Y_e*prims.temperature) )
              check = 1;
          }
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
          } // atmospheric backup
        } // if c2p failed
        /***************************************************************/

        //--------------------------------------------------
        //---------- Primitive recovery succeeded ----------
        //--------------------------------------------------
        // Enforce limits on primitive variables and recompute conservatives.
        diagnostics.speed_limited += ghl_enforce_primitive_limits_and_compute_u0(
              ghl_params, ghl_eos, &ADM_metric, &prims);

        rho[index]         = prims.rho;
        press[index]       = prims.press;
        eps[index]         = prims.eps;
        u0[index]          = prims.u0;
        vx[index]          = prims.vU[0];
        vy[index]          = prims.vU[1];
        vz[index]          = prims.vU[2];
        entropy[index]     = prims.entropy;
        Y_e[index]         = prims.Y_e;
        temperature[index] = prims.temperature;

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
  CCTK_REAL error_ent_numer = 0;
  CCTK_REAL error_Ye_numer = 0;

  CCTK_REAL error_rho_denom = 0;
  CCTK_REAL error_tau_denom = 0;
  CCTK_REAL error_Sx_denom = 0;
  CCTK_REAL error_Sy_denom = 0;
  CCTK_REAL error_Sz_denom = 0;
  CCTK_REAL error_ent_denom = 0;
  CCTK_REAL error_Ye_denom = 0;

// We split the prim2con into a separate loop so the averaging
// loop is deterministic and doesn't have a race condition
#pragma omp parallel for reduction(+:                                   \
      error_rho_numer, error_tau_numer, error_Sx_numer, error_Sy_numer, \
      error_Sz_numer, error_rho_denom, error_tau_denom, error_Sx_denom, \
      error_Sy_denom, error_Sz_denom, error_ent_numer, error_ent_denom, \
      error_Ye_numer, error_Ye_denom) \
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
        prims.rho         = rho[index];
        prims.press       = press[index];
        prims.eps         = eps[index];
        prims.u0          = u0[index];
        prims.vU[0]       = vx[index];
        prims.vU[1]       = vy[index];
        prims.vU[2]       = vz[index];
        prims.BU[0]       = Bx_center[index];
        prims.BU[1]       = By_center[index];
        prims.BU[2]       = Bz_center[index];
        prims.entropy     = entropy[index];
        prims.Y_e         = Y_e[index];
        prims.temperature = temperature[index];

        ghl_conservative_quantities cons, cons_orig;
        cons_orig.rho     = rho_star[index];
        cons_orig.tau     = tau[index];
        cons_orig.SD[0]   = Stildex[index];
        cons_orig.SD[1]   = Stildey[index];
        cons_orig.SD[2]   = Stildez[index];
        cons_orig.entropy = ent_star[index];
        cons_orig.Y_e     = Ye_star[index];

        ghl_compute_conservs(&ADM_metric, &metric_aux, &prims, &cons);

        rho_star[index] = cons.rho;
        tau[index]      = cons.tau;
        Stildex[index]  = cons.SD[0];
        Stildey[index]  = cons.SD[1];
        Stildez[index]  = cons.SD[2];
        ent_star[index] = cons.entropy;
        Ye_star[index]  = cons.Y_e;

        // Now we compute the difference between original & new conservatives
        // for diagnostic purposes:
        error_rho_numer += fabs(cons.rho - cons_orig.rho);
        error_tau_numer += fabs(cons.tau - cons_orig.tau);
        error_Sx_numer  += fabs(cons.SD[0] - cons_orig.SD[0]);
        error_Sy_numer  += fabs(cons.SD[1] - cons_orig.SD[1]);
        error_Sz_numer  += fabs(cons.SD[2] - cons_orig.SD[2]);
        error_ent_numer += fabs(cons.entropy - cons_orig.entropy);
        error_Ye_numer  += fabs(cons.Y_e - cons_orig.Y_e);
        error_rho_denom += cons_orig.rho;
        error_tau_denom += cons_orig.tau;
        error_Sx_denom  += fabs(cons_orig.SD[0]);
        error_Sy_denom  += fabs(cons_orig.SD[1]);
        error_Sz_denom  += fabs(cons_orig.SD[2]);
        error_ent_denom += cons_orig.entropy;
        error_Ye_denom  += cons_orig.Y_e;
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
    const CCTK_REAL ent_error = (error_ent_denom == 0) ? error_ent_numer : error_ent_numer/error_ent_denom;
    const CCTK_REAL Ye_error  = (error_Ye_denom == 0)  ? error_Ye_numer  : error_Ye_numer/error_Ye_denom;

    CCTK_VINFO(
        "C2P: Iter. # %d, Lev: %d NumPts= %d | Backups: %d %d %d | Fixes: VL= %d rho*= %d\n"
        "                 Averaged pts = %d | Failures: %d InHoriz= %d / %d | %.2f iters/gridpt\n"
        "   Error, Sum: rho %.3e, %.3e | tau %.3e, %.3e | entropy %.3e, %.3e | "
        "Y_e %.3e, %.3e\n"
        "               Sx %.3e, %.3e | Sy %.3e, %.3e | Sz %.3e, %.3e\n",
        cctk_iteration, GetRefinementLevel(cctkGH), pointcount,
        backup0, backup1, backup2, vel_limited_ptcount,
        rho_star_fix_applied, pointcount_avg,
        failures, failures_inhoriz, pointcount_inhoriz,
        (CCTK_REAL)n_iter / ((CCTK_REAL)(pointcount)),
        rho_error, error_rho_denom, tau_error, error_tau_denom,
        ent_error, error_ent_denom, Ye_error, error_Ye_denom,
        Sx_error, error_Sx_denom, Sy_error, error_Sy_denom,
        Sz_error, error_Sz_denom);
  }
}
