/********************************
 * CONVERT ET ID TO IllinoisGRMHD
 *
 * Written in 2014 by Zachariah B. Etienne
 *
 * Sets metric & MHD variables needed
 * by IllinoisGRMHD, converting from
 * HydroBase and ADMBase.
 ********************************/

#include "GRHayLHD.h"

void convert_HydroBase_to_GRHayLHD(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_convert_HydroBase_to_GRHayLHD;
  DECLARE_CCTK_PARAMETERS;

  const int imax = cctk_lsh[0];
  const int jmax = cctk_lsh[1];
  const int kmax = cctk_lsh[2];

#pragma omp parallel for
  for (int k = 0; k < kmax; k++) {
    for (int j = 0; j < jmax; j++) {
      for (int i = 0; i < imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);
        const int ind0 = CCTK_VECTGFINDEX3D(cctkGH, i, j, k, 0);
        const int ind1 = CCTK_VECTGFINDEX3D(cctkGH, i, j, k, 1);
        const int ind2 = CCTK_VECTGFINDEX3D(cctkGH, i, j, k, 2);

        const double ETvx = vel[ind0];
        const double ETvy = vel[ind1];
        const double ETvz = vel[ind2];

        // GRHayLHD defines v^i = u^i/u^0.

        // Meanwhile, the ET/HydroBase formalism, called the Valencia
        // formalism, splits the 4 velocity into a purely spatial part
        // and a part that is normal to the spatial hypersurface:
        // u^a = G (n^a + U^a), (Eq. 14 of arXiv:1304.5544; G=W, U^a=v^a)
        // where n^a is the unit normal vector to the spatial hypersurface,
        // n_a = {-\alpha,0,0,0}, and U^a is the purely spatial part, which
        // is defined in HydroBase as the vel[] vector gridfunction.
        // Then u^a n_a = - \alpha u^0 = G n^a n_a = -G, and
        // of course \alpha u^0 = 1/sqrt(1+γ^ij u_i u_j) = \Gamma,
        // the standard Lorentz factor.

        // Note that n^i = - \beta^i / \alpha, so
        // u^a = \Gamma (n^a + U^a)
        // -> u^i = \Gamma ( U^i - \beta^i / \alpha )
        // which implies
        // v^i = u^i/u^0
        //     = \Gamma/u^0 ( U^i - \beta^i / \alpha ) <- \Gamma = \alpha u^0
        //     = \alpha ( U^i - \beta^i / \alpha )
        //     = \alpha U^i - \beta^i

        vx[index] = alp[index] * ETvx - betax[index];
        vy[index] = alp[index] * ETvy - betay[index];
        vz[index] = alp[index] * ETvz - betaz[index];
      }
    }
  }
}
