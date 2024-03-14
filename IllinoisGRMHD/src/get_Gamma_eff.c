#include "IllinoisGRMHD.h"

CCTK_REAL get_Gamma_eff_hybrid(
      const CCTK_REAL rho_in,
      const CCTK_REAL press_in) {
  CCTK_REAL K, Gamma;
  ghl_hybrid_get_K_and_Gamma(ghl_eos, rho_in, &K, &Gamma);
  const CCTK_REAL P_cold = K*pow(rho_in, Gamma);
  return ghl_eos->Gamma_th + (Gamma - ghl_eos->Gamma_th)*P_cold/press_in;
}

CCTK_REAL get_Gamma_eff_tabulated(
      const CCTK_REAL rho_in,
      const CCTK_REAL press_in) {
  return 1.0;
}
