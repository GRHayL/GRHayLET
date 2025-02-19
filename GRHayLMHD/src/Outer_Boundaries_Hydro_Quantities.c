#include "GRHayLMHD.h"

void GRHayLMHD_Outer_Boundaries_Hydro_Quantities(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // TODO: add BC options
    if(cctk_iteration == 0 || GetRefinementLevel(cctkGH)) {
        return;
    }

    // clang-format off
    const int lower[3] = {
        cctk_nghostzones[0],
        cctk_nghostzones[1],
        cctk_nghostzones[2]
    };
    const int upper[3] = {
        cctk_lsh[0] - cctk_nghostzones[0] - 1,
        cctk_lsh[1] - cctk_nghostzones[1] - 1,
        cctk_lsh[2] - cctk_nghostzones[2] - 1
    };
    // clang-format on

    for(int gzidx = 0; gzidx < cctk_nghostzones[0]; gzidx++) {
        const int imin = lower[0] - gzidx;
        const int jmin = lower[1] - gzidx;
        const int kmin = lower[2] - gzidx;
        const int imax = upper[0] + gzidx;
        const int jmax = upper[1] + gzidx;
        const int kmax = upper[2] + gzidx;

        OBLOOPX(GRHayLMHD_xmin)
        {
            rho[dst_idx]         = rho[src_idx];
            press[dst_idx]       = press[src_idx];
            eps[dst_idx]         = eps[src_idx];
            Y_e[dst_idx]         = Y_e[src_idx];
            temperature[dst_idx] = temperature[src_idx];
            entropy[dst_idx]     = entropy[src_idx];
            vx[dst_idx]          = vx[src_idx] < 0 ? vx[src_idx] : 0;
            vy[dst_idx]          = vy[src_idx];
            vz[dst_idx]          = vz[src_idx];
        }
        ENDLOOP

        OBLOOPX(GRHayLMHD_xmax)
        {
            rho[dst_idx]         = rho[src_idx];
            press[dst_idx]       = press[src_idx];
            eps[dst_idx]         = eps[src_idx];
            Y_e[dst_idx]         = Y_e[src_idx];
            temperature[dst_idx] = temperature[src_idx];
            entropy[dst_idx]     = entropy[src_idx];
            vx[dst_idx]          = vx[src_idx] > 0 ? vx[src_idx] : 0;
            vy[dst_idx]          = vy[src_idx];
            vz[dst_idx]          = vz[src_idx];
        }
        ENDLOOP

        OBLOOPY(GRHayLMHD_ymin)
        {
            rho[dst_idx]         = rho[src_idx];
            press[dst_idx]       = press[src_idx];
            eps[dst_idx]         = eps[src_idx];
            Y_e[dst_idx]         = Y_e[src_idx];
            temperature[dst_idx] = temperature[src_idx];
            entropy[dst_idx]     = entropy[src_idx];
            vx[dst_idx]          = vx[src_idx];
            vy[dst_idx]          = vy[src_idx] < 0 ? vy[src_idx] : 0;
            vz[dst_idx]          = vz[src_idx];
        }
        ENDLOOP

        OBLOOPY(GRHayLMHD_ymin)
        {
            rho[dst_idx]         = rho[src_idx];
            press[dst_idx]       = press[src_idx];
            eps[dst_idx]         = eps[src_idx];
            Y_e[dst_idx]         = Y_e[src_idx];
            temperature[dst_idx] = temperature[src_idx];
            entropy[dst_idx]     = entropy[src_idx];
            vx[dst_idx]          = vx[src_idx];
            vy[dst_idx]          = vy[src_idx] > 0 ? vy[src_idx] : 0;
            vz[dst_idx]          = vz[src_idx];
        }
        ENDLOOP

        OBLOOPZ(GRHayLMHD_zmin)
        {
            rho[dst_idx]         = rho[src_idx];
            press[dst_idx]       = press[src_idx];
            eps[dst_idx]         = eps[src_idx];
            Y_e[dst_idx]         = Y_e[src_idx];
            temperature[dst_idx] = temperature[src_idx];
            entropy[dst_idx]     = entropy[src_idx];
            vx[dst_idx]          = vx[src_idx];
            vy[dst_idx]          = vy[src_idx];
            vz[dst_idx]          = vz[src_idx] < 0 ? vz[src_idx] : 0;
        }
        ENDLOOP

        OBLOOPZ(GRHayLMHD_zmin)
        {
            rho[dst_idx]         = rho[src_idx];
            press[dst_idx]       = press[src_idx];
            eps[dst_idx]         = eps[src_idx];
            Y_e[dst_idx]         = Y_e[src_idx];
            temperature[dst_idx] = temperature[src_idx];
            entropy[dst_idx]     = entropy[src_idx];
            vx[dst_idx]          = vx[src_idx];
            vy[dst_idx]          = vy[src_idx];
            vz[dst_idx]          = vz[src_idx] > 0 ? vz[src_idx] : 0;
        }
        ENDLOOP
    }
}
