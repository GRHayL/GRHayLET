#include "IllinoisGRMHD.h"

void IllinoisGRMHD_tabulated_perturb_primitives(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_IllinoisGRMHD_tabulated_perturb_primitives;
  DECLARE_CCTK_PARAMETERS;

  const int imax = cctk_lsh[0];
  const int jmax = cctk_lsh[1];
  const int kmax = cctk_lsh[2];

#pragma omp parallel for
  for(int k=0; k<kmax; k++) {
    for(int j=0; j<jmax; j++) {
      for(int i=0; i<imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);
        const int gi = cctk_lbnd[0] + i;
        const int gj = cctk_lbnd[1] + j;
        const int gk = cctk_lbnd[2] + k;
        rho[index]         *= IllinoisGRMHD_one_plus_pert(random_pert, random_seed, gi, gj, gk, 0);
        press[index]       *= IllinoisGRMHD_one_plus_pert(random_pert, random_seed, gi, gj, gk, 1);
        vx[index]          *= IllinoisGRMHD_one_plus_pert(random_pert, random_seed, gi, gj, gk, 2);
        vy[index]          *= IllinoisGRMHD_one_plus_pert(random_pert, random_seed, gi, gj, gk, 3);
        vz[index]          *= IllinoisGRMHD_one_plus_pert(random_pert, random_seed, gi, gj, gk, 4);
        Y_e[index]         *= IllinoisGRMHD_one_plus_pert(random_pert, random_seed, gi, gj, gk, 5);
        temperature[index] *= IllinoisGRMHD_one_plus_pert(random_pert, random_seed, gi, gj, gk, 6);

        phitilde[index] *= IllinoisGRMHD_one_plus_pert(random_pert, random_seed, gi, gj, gk, 7);
        Ax[index]       *= IllinoisGRMHD_one_plus_pert(random_pert, random_seed, gi, gj, gk, 8);
        Ay[index]       *= IllinoisGRMHD_one_plus_pert(random_pert, random_seed, gi, gj, gk, 9);
        Az[index]       *= IllinoisGRMHD_one_plus_pert(random_pert, random_seed, gi, gj, gk, 10);
      }
    }
  }
}
