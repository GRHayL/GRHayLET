#include "GRHayLMHD.h"

extern "C" void GRHayLMHD_Convert_from_HydroBase(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    const int imin = 0, imax = cctk_lsh[0];
    const int jmin = 0, jmax = cctk_lsh[1];
    const int kmin = 0, kmax = cctk_lsh[2];

    OMPLOOP3D(imin, imax, jmin, jmax, kmin, kmax)
    {
        const CCTK_REAL alpL = alp[index];
        const CCTK_REAL gxxL = gxx[ijk];
        const CCTK_REAL gxyL = gxy[ijk];
        const CCTK_REAL gxzL = gxz[ijk];
        const CCTK_REAL gyyL = gyy[ijk];
        const CCTK_REAL gyzL = gyz[ijk];
        const CCTK_REAL gzzL = gzz[ijk];

        vx[ijk] = alpL * vel[ijkx] - betax[ijk];
        vy[ijk] = alpL * vel[ijky] - betax[ijk];
        vz[ijk] = alpL * vel[ijkz] - betax[ijk];

        A_x[ijk]  = Avec[ijkx];
        A_y[ijk]  = Avec[ijky];
        A_z[ijk]  = Avec[ijkz];
        Phi_tilde = Aphi[ijk];
    }
    ENDLOOP3D
}
