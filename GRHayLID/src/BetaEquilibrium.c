#include "GRHayLID.h"

/* Use base chemical potentials, avoiding shifted munu in affected LS tables. */
static double GRHayLID_beta_residual(const double density, const double ye, const double temp) {
  double pressure, energy, muhat, mu_e, mu_p, mu_n;
  const ghl_error_codes_t err = ghl_tabulated_compute_P_eps_muhat_mue_mup_mun_from_T(
      ghl_eos, density, ye, temp, &pressure, &energy, &muhat, &mu_e, &mu_p, &mu_n);
  if(err != ghl_success || !isfinite(mu_e) || !isfinite(mu_p) || !isfinite(mu_n))
    CCTK_VERROR("Beta chemical-potential lookup failed at rho=%g, Ye=%g, T=%g (status %d)",
                density, ye, temp, (int)err);
  const double residual = mu_e - mu_n + mu_p;
  if(!isfinite(residual))
    CCTK_ERROR("Nonfinite beta-equilibrium residual");
  return residual;
}

/* Search each admissible Ye interval at the actual rho and T, including nodes.
 * Choose the first root in increasing Ye; a missing root is an error. */
static double GRHayLID_beta_root(const double density, const double temp, const double tolerance) {
  double lower = ghl_eos->Y_e_min;
  double f_lower = GRHayLID_beta_residual(density, lower, temp);
  if(fabs(f_lower) <= tolerance) return lower;
  for(int n=0; n<=ghl_eos->N_Ye; n++) {
    const double upper = n == ghl_eos->N_Ye ? ghl_eos->Y_e_max : ghl_eos->table_Y_e[n];
    if(upper <= lower || upper > ghl_eos->Y_e_max) continue;
    double f_upper = GRHayLID_beta_residual(density, upper, temp);
    if(fabs(f_upper) <= tolerance) return upper;
    if(signbit(f_lower) != signbit(f_upper)) {
      double a = lower, b = upper;
      for(int iteration=0; iteration<80; iteration++) {
        const double mid = a + 0.5*(b-a);
        const double f_mid = GRHayLID_beta_residual(density, mid, temp);
        if(fabs(f_mid) <= tolerance) return mid;
        if(mid == a || mid == b) break;
        if(signbit(f_mid) == signbit(f_lower)) {
          a = mid;
          f_lower = f_mid;
        } else {
          b = mid;
        }
      }
      const double f_b = GRHayLID_beta_residual(density, b, temp);
      CCTK_VERROR("Beta root did not meet residual tolerance %g MeV at rho=%g, T=%g: bisection stalled on Ye in [%.17g, %.17g] with residuals %g and %g MeV; check EOS interpolation continuity and the required residual tolerance",
                  tolerance, density, temp, a, b, f_lower, f_b);
    }
    lower = upper;
    f_lower = f_upper;
  }
  CCTK_VERROR("No beta-equilibrium root within effective Ye bounds at rho=%g, T=%g", density, temp);
  return 0.0; // CCTK_VERROR terminates.
}

void GRHayLID_BetaEquilibrium(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLID_BetaEquilibrium;
  DECLARE_CCTK_PARAMETERS;
  CCTK_INFO("Starting routine to impose neutrino free beta-equilibrium");
  GRHayLID_require_storage(cctkGH, "HydroBase::Y_e");
  GRHayLID_require_storage(cctkGH, "HydroBase::temperature");
  if(ghl_eos->eos_type != ghl_eos_tabulated)
    CCTK_ERROR("Beta equilibrium requires Tabulated EOS");
  CHECK_PARAMETER(beq_temperature);
  if(!isfinite(beq_temperature) || beq_temperature <= 0.0 ||
     beq_temperature < ghl_eos->T_min || beq_temperature > ghl_eos->T_max)
    CCTK_ERROR("beq_temperature must be positive, finite, and within effective EOS bounds");
  if(!isfinite(beq_residual_tolerance) || beq_residual_tolerance <= 0.0)
    CCTK_ERROR("beq_residual_tolerance must be positive and finite");
  GRHayLID_check_table_state(ghl_eos->rho_atm, ghl_eos->Y_e_atm, ghl_eos->T_atm, "beta atmosphere");

  for(int k=0; k<cctk_lsh[2]; k++) {
    for(int j=0; j<cctk_lsh[1]; j++) {
      for(int i=0; i<cctk_lsh[0]; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);
        if(!isfinite(rho[index]))
          CCTK_ERROR("Beta equilibrium requires finite initialized density");
        const int atmosphere = rho[index] <= 1.01*ghl_eos->rho_atm;
        const double density = atmosphere ? ghl_eos->rho_atm : rho[index];
        const double temp = atmosphere ? ghl_eos->T_atm : beq_temperature;
        if(density < ghl_eos->rho_min || density > ghl_eos->rho_max || density <= 0.0)
          CCTK_ERROR("Beta density lies outside effective EOS bounds");
        const double ye = atmosphere ? ghl_eos->Y_e_atm :
                                      GRHayLID_beta_root(density, temp, beq_residual_tolerance);
        GRHayLID_check_table_state(density, ye, temp, "beta output");
        double pressure, energy;
        const ghl_error_codes_t err = ghl_tabulated_compute_P_eps_from_T(
            ghl_eos, density, ye, temp, &pressure, &energy);
        if(err != ghl_success || !isfinite(pressure) || !isfinite(energy))
          CCTK_VERROR("Beta thermodynamic lookup failed (status %d)", (int)err);
        if(!atmosphere && fabs(GRHayLID_beta_residual(density, ye, temp)) > beq_residual_tolerance)
          CCTK_ERROR("Final beta state fails its chemical-potential residual tolerance");
        rho[index] = density;
        press[index] = pressure;
        eps[index] = energy;
        Y_e[index] = ye;
        temperature[index] = temp;
      }
    }
  }
  CCTK_INFO("Finished imposing neutrino free beta-equilibrium");
}
