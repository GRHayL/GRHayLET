#include "IllinoisGRMHD.h"

static ghl_error_codes_t (*original_compute_h_and_cs2)(
    const ghl_eos_parameters *restrict, ghl_primitive_quantities *restrict,
    double *restrict, double *restrict);

// The direct Gamma-law formulas avoid subtraction of the auxiliary cold curve.
ghl_error_codes_t IllinoisGRMHD_compute_h_and_cs2(
    const ghl_eos_parameters *restrict eos,
    ghl_primitive_quantities *restrict prims,
    double *restrict h, double *restrict cs2) {
  if(eos->eos_type != ghl_eos_simple)
    return original_compute_h_and_cs2(eos, prims, h, cs2);

  const double p_over_rho = prims->press/prims->rho;
  prims->eps = p_over_rho/(eos->Gamma_th - 1.0);
  *h = 1.0 + prims->eps + p_over_rho;
  *cs2 = eos->Gamma_th*p_over_rho/(*h);
  return ghl_success;
}

void IllinoisGRMHD_check_eos_support(CCTK_ARGUMENTS) {
  if(!ghl_params || !ghl_eos)
    CCTK_ERROR("GRHayL must be initialized before IllinoisGRMHD checks EOS support.");

  if(ghl_params->evolve_entropy) {
    if(ghl_eos->eos_type == ghl_eos_tabulated)
      CCTK_ERROR("Tabulated entropy evolution is unsupported until GRHayL packing, fluxes, recovery, and atmosphere use a density-weighted specific entropy current. Use evolve_entropy=no.");
    if(ghl_eos->eos_type == ghl_eos_hybrid &&
       (ghl_eos->neos != 1 || ghl_eos->Gamma_th != ghl_eos->Gamma_ppoly[0]))
      CCTK_ERROR("Hybrid entropy evolution requires neos=1 and Gamma_th=Gamma_ppoly_in[0]. Use evolve_entropy=no for a general hybrid EOS.");
  }

  // GRHayL uses this callback in source, wave-speed, and flux calculations.
  // Preserve other EOS callbacks; install once after modern or legacy startup.
  if(ghl_compute_h_and_cs2 != IllinoisGRMHD_compute_h_and_cs2) {
    original_compute_h_and_cs2 = ghl_compute_h_and_cs2;
    ghl_compute_h_and_cs2 = IllinoisGRMHD_compute_h_and_cs2;
  }
}
