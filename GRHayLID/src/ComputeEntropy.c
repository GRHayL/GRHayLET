#include "GRHayLID.h"

void GRHayLID_compute_entropy_hybrid(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLID_compute_entropy_hybrid;
  DECLARE_CCTK_PARAMETERS;

#pragma omp parallel for
  for(int k=0; k<cctk_lsh[2]; k++) {
    for(int j=0; j<cctk_lsh[1]; j++) {
      for(int i=0; i<cctk_lsh[0]; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        entropy[index] = ghl_hybrid_compute_entropy_function(ghl_eos, rho[index], press[index]);
      }
    }
  }
}

void GRHayLID_compute_entropy_tabulated(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLID_compute_entropy_tabulated;
  DECLARE_CCTK_PARAMETERS;

#pragma omp parallel for
  for(int k=0; k<cctk_lsh[2]; k++) {
    for(int j=0; j<cctk_lsh[1]; j++) {
      for(int i=0; i<cctk_lsh[0]; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        ghl_tabulated_enforce_bounds_rho_Ye_T(ghl_eos, &rho[index], &Y_e[index], &temperature[index]);
        ghl_tabulated_compute_P_eps_S_from_T(
              ghl_eos, rho[index], Y_e[index], temperature[index],
              &press[index], &eps[index], &entropy[index]);
      }
    }
  }
}
