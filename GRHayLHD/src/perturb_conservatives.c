#include "GRHayLHD.h"

void GRHayLHD_perturb_conservatives(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLHD_perturb_conservatives;
  DECLARE_CCTK_PARAMETERS;

  const int imax = cctk_lsh[0];
  const int jmax = cctk_lsh[1];
  const int kmax = cctk_lsh[2];

  srand(random_seed); // Use srand() as rand() is thread-safe.
#pragma omp parallel for
  for(int k=0; k<kmax; k++) {
    for(int j=0; j<jmax; j++) {
      for(int i=0; i<imax; i++) {
        const int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
        rho_star[index] *= one_plus_pert(random_pert);
        tau[index]      *= one_plus_pert(random_pert);
        Stildex[index]  *= one_plus_pert(random_pert);
        Stildey[index]  *= one_plus_pert(random_pert);
        Stildez[index]  *= one_plus_pert(random_pert);
      }
    }
  }
}
