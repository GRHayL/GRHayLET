#include "GRHayLMHD.h"

void GRHayLMHD_convert_for_NRPyLeakage(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLMHD_convert_for_NRPyLeakage;
  DECLARE_CCTK_PARAMETERS;

  // Convert rho, Y_e, T, and velocities to HydroBase
  // because they are needed by NRPyLeakage
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) {
    for(int j=0;j<cctk_lsh[1];j++) {
      for(int i=0;i<cctk_lsh[0];i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);
        const int idxvx = CCTK_VECTGFINDEX3D(cctkGH, i, j, k, 0);
        const int idxvy = CCTK_VECTGFINDEX3D(cctkGH, i, j, k, 1);
        const int idxvz = CCTK_VECTGFINDEX3D(cctkGH, i, j, k, 2);

        // Read from main memory
        const CCTK_REAL invalpL      = 1.0/alp[index];
        const CCTK_REAL betaxL       = betax[index];
        const CCTK_REAL betayL       = betay[index];
        const CCTK_REAL betazL       = betaz[index];
        const CCTK_REAL rhoL         = rho_b[index];
        const CCTK_REAL vxL          = vx[index];
        const CCTK_REAL vyL          = vy[index];
        const CCTK_REAL vzL          = vz[index];

        // Write to main memory, converting to HydroBase
        rho [index] = rhoL;
        vel [idxvx] = (vxL + betaxL)*invalpL;
        vel [idxvy] = (vyL + betayL)*invalpL;
        vel [idxvz] = (vzL + betazL)*invalpL;
      }
    }
  }
}
