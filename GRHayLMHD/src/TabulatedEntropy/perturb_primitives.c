#include "GRHayLMHD.h"

void GRHayLMHD_tabulated_entropy_perturb_primitives(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLMHD_tabulated_entropy_perturb_primitives;
  DECLARE_CCTK_PARAMETERS;

  const int imax = cctk_lsh[0];
  const int jmax = cctk_lsh[1];
  const int kmax = cctk_lsh[2];

  srand(random_seed); // Use srand() as rand() is thread-safe.
#pragma omp parallel for
  for(int k=0; k<kmax; k++) {
    for(int j=0; j<jmax; j++) {
      for(int i=0; i<imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);
        rho[index]         *= one_plus_pert(random_pert);
        press[index]       *= one_plus_pert(random_pert);
        vx[index]          *= one_plus_pert(random_pert);
        vy[index]          *= one_plus_pert(random_pert);
        vz[index]          *= one_plus_pert(random_pert);
        entropy[index]     *= one_plus_pert(random_pert);
        Y_e[index]         *= one_plus_pert(random_pert);
        temperature[index] *= one_plus_pert(random_pert);

        phitilde[index] *= one_plus_pert(random_pert);
        Ax[index]       *= one_plus_pert(random_pert);
        Ay[index]       *= one_plus_pert(random_pert);
        Az[index]       *= one_plus_pert(random_pert);
      }
    }
  }
}
