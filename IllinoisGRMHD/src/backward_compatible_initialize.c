#include "IllinoisGRMHD.h"

void IllinoisGRMHD_backward_compatible_initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  CCTK_INFO(
        "Warning!! Many features in IllinoisGRMHD have been improved alongside\n"
        "the development of GRHayL. Please update your parfile to use the new features.\n"
        "Running GRHayLib with backward-compatible settings...");

  if(neos>1)
    CCTK_VERROR(
          "Error: neos %d > 1. Original IllinoisGRMHD did not properly implement piecewise polytrope,\n"
          "so this must be updated to use the modern version of the thorn for reliable behavior.", neos);

  ghl_params = (ghl_parameters *)malloc(sizeof(ghl_parameters));
  ghl_eos = (ghl_eos_parameters *)malloc(sizeof(ghl_eos_parameters));

  const int main = Noble2D;
  const int backups[3] = {Font1D, None, None};

  const bool evolve_temperature = false;
  const bool calc_primitive_guess = true;

  ghl_initialize_params(
      main, backups,
      evolve_entropy, evolve_temperature,
      calc_primitive_guess, Psi6threshold,
      GAMMA_SPEED_LIMIT, damp_lorenz,
      ghl_params);

      ghl_params->ppm_flattening_epsilon = 0.33;
      ghl_params->ppm_flattening_omega1  = 0.75;
      ghl_params->ppm_flattening_omega2  = 10.0;

      ghl_params->ppm_shock_epsilon = 0.01;
      ghl_params->ppm_shock_eta1    = 20.0;
      ghl_params->ppm_shock_eta2    = 0.05;
      ghl_params->ppm_shock_k0      = 0.1;

  const double rho_ppoly_in[1] = {0.0};
  const double Gamma_ppoly_in[1] = {gamma_th};

  ghl_con2prim_multi_method = ghl_con2prim_hybrid_multi_method;
  ghl_initialize_hybrid_eos_functions_and_params(
        rho_b_atm, rho_b_atm, rho_b_max,
        neos, rho_ppoly_in,
        Gamma_ppoly_in, K_poly,
        gamma_th, ghl_eos);
}
