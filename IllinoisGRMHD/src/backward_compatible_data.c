#include "IllinoisGRMHD.h"

void IllinoisGRMHD_backward_compatible_data(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_IllinoisGRMHD_backward_compatible_data;
  DECLARE_CCTK_PARAMETERS;

#pragma omp parallel for
  for(int k=0; k<cctk_lsh[2]; k++) {
    for(int j=0; j<cctk_lsh[1]; j++) {
      for(int i=0; i<cctk_lsh[0]; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);

        rho_b[index] = rho[index];
        P[index] = press[index];
        Bx[index] = Bx_center[index];
        By[index] = By_center[index];
        Bz[index] = Bz_center[index];
        psi6phi[index] = phitilde[index];
      }
    }
  }
}
