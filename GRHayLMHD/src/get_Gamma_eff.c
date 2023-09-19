#include "GRHayLMHD.h"

double get_Gamma_eff_hybrid(
      const double rho_in,
      const double press_in) {
  double K, Gamma;
  ghl_hybrid_get_K_and_Gamma(ghl_eos, rho_in, &K, &Gamma);
  const double P_cold = K*pow(rho_in, Gamma);
  return ghl_eos->Gamma_th + (Gamma - ghl_eos->Gamma_th)*P_cold/press_in;
}

double get_Gamma_eff_tabulated(
      const double rho_in,
      const double press_in) {
  return 1.0;
}
