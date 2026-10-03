#include "IllinoisGRMHD.h"

void IllinoisGRMHD_hybrid_entropy_perturb_conservatives(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_IllinoisGRMHD_hybrid_entropy_perturb_conservatives;
  DECLARE_CCTK_PARAMETERS;

  const int imax = cctk_lsh[0];
  const int jmax = cctk_lsh[1];
  const int kmax = cctk_lsh[2];

#pragma omp parallel for
  for(int k=0; k<kmax; k++) {
    for(int j=0; j<jmax; j++) {
      for(int i=0; i<imax; i++) {
        const int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
        const int gi = cctk_lbnd[0] + i;
        const int gj = cctk_lbnd[1] + j;
        const int gk = cctk_lbnd[2] + k;
        rho_star[index] *= IllinoisGRMHD_one_plus_pert(random_pert, random_seed, gi, gj, gk, 0);
        tau[index]      *= IllinoisGRMHD_one_plus_pert(random_pert, random_seed, gi, gj, gk, 1);
        Stildex[index]  *= IllinoisGRMHD_one_plus_pert(random_pert, random_seed, gi, gj, gk, 2);
        Stildey[index]  *= IllinoisGRMHD_one_plus_pert(random_pert, random_seed, gi, gj, gk, 3);
        Stildez[index]  *= IllinoisGRMHD_one_plus_pert(random_pert, random_seed, gi, gj, gk, 4);
        ent_star[index] *= IllinoisGRMHD_one_plus_pert(random_pert, random_seed, gi, gj, gk, 5);
      }
    }
  }
}
