#include "GRHayLMHD.h"

void GRHayLMHD_Convert_from_HydroBase(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    const int imin = 0, imax = cctk_lsh[0];
    const int jmin = 0, jmax = cctk_lsh[1];
    const int kmin = 0, kmax = cctk_lsh[2];

    OMPLOOP3D(imin, imax, jmin, jmax, kmin, kmax)
    {
        const CCTK_REAL alpL = alp[ijk];

        vx[ijk] = alpL * vel[ijkx] - betax[ijk];
        vy[ijk] = alpL * vel[ijky] - betay[ijk];
        vz[ijk] = alpL * vel[ijkz] - betaz[ijk];

        // A_x_tilde[ijk] = Avec[ijkx];
        // A_y_tilde[ijk] = Avec[ijky];
        // A_z_tilde[ijk] = Avec[ijkz];
        // Phi_tilde[ijk] = Aphi[ijk];
    }
    ENDLOOP3D
}
