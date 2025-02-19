#include "GRHayLMHD.h"

extern "C" void GRHayLMHD_Convert_to_HydroBase(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    bool convert_enabled = Convert_to_HydroBase_every != 0;
    bool time_to_convert = (cctk_iteration % Convert_to_HydroBase_every) == 0;

    if(!convert_enabled || !time_to_convert) {
        return;
    }

    const int imin = 0, imax = cctk_lsh[0];
    const int jmin = 0, jmax = cctk_lsh[1];
    const int kmin = 0, kmax = cctk_lsh[2];

    OMPLOOP3D(imin, imax, jmin, jmax, kmin, kmax)
    {
        const CCTK_REAL inv_lapse = 1.0 / alp[index];
        const CCTK_REAL gxxL      = gxx[ijk];
        const CCTK_REAL gxyL      = gxy[ijk];
        const CCTK_REAL gxzL      = gxz[ijk];
        const CCTK_REAL gyyL      = gyy[ijk];
        const CCTK_REAL gyzL      = gyz[ijk];
        const CCTK_REAL gzzL      = gzz[ijk];
        const CCTK_REAL velx      = (vx[ijk] + betax[ijk]) * inv_lapse;
        const CCTK_REAL vely      = (vy[ijk] + betay[ijk]) * inv_lapse;
        const CCTK_REAL velz      = (vz[ijk] + betaz[ijk]) * inv_lapse;

        // clang-format off
        const CCTK_REAL vsqr = gxx * velx * velx
                             + gxy * velx * vely * 2
                             + gxz * velx * velz * 2
                             + gyy * vely * vely
                             + gyz * vely * velz * 2
                             + gzz * velz * velz;
        // clang-format on

        vel[ijkx]      = velx;
        vel[ijky]      = vely;
        vel[ijkz]      = velz;
        w_lorentz[ijk] = 1.0 / sqrt(1.0 - vsqr);
    }
    ENDLOOP3D
}
