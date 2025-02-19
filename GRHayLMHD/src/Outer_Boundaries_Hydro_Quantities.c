#include "GRHayLMHD.h"

void GRHayLMHD_Outer_Boundaries_Hydro_Quantities(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // TODO: add BC options
    if(cctk_iteration == 0 || GetRefinementLevel(cctkGH)) {
        return;
    }

    const int lower[3] = {
        cctk_nghostzones[0],
        cctk_nghostzones[1],
        cctk_nghostzones[2],
    };
    const int upper[3] = {
        cctk_lsh[0] - cctk_nghostzones[0] - 1,
        cctk_lsh[1] - cctk_nghostzones[1] - 1,
        cctk_lsh[2] - cctk_nghostzones[2] - 1,
    };

    for(int gzidx = 0; gzidx < cctk_nghostzones[0]; gzidx++) {
        const int imin = lower[0] - gzidx;
        const int jmin = lower[1] - gzidx;
        const int kmin = lower[2] - gzidx;
        const int imax = upper[0] + gzidx;
        const int jmax = upper[1] + gzidx;
        const int kmax = upper[2] + gzidx;

        if(cctk_bbox[GRHayLMHD_xmin]) {
            for(int k = 0; k < cctk_lsh[2]; k++) {
                for(int j = 0; j < cctk_lsh[1]; j++) {
                    const int src_idx = CCTK_GFINDEX3D(cctkGH, imin + 1, j, k);
                    const int dst_idx = CCTK_GFINDEX3D(cctkGH, imin, j, k);

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
            }
        }

        if(cctk_bbox[GRHayLMHD_xmax]) {
            for(int k = 0; k < cctk_lsh[2]; k++) {
                for(int j = 0; j < cctk_lsh[1]; j++) {
                    const int src_idx = CCTK_GFINDEX3D(cctkGH, imax - 1, j, k);
                    const int dst_idx = CCTK_GFINDEX3D(cctkGH, imax, j, k);

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
            }
        }

        if(cctk_bbox[GRHayLMHD_ymin]) {
            for(int k = 0; k < cctk_lsh[2]; k++) {
                for(int i = 0; i < cctk_lsh[0]; i++) {
                    const int src_idx = CCTK_GFINDEX3D(cctkGH, i, jmin + 1, k);
                    const int dst_idx = CCTK_GFINDEX3D(cctkGH, i, jmin, k);

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
            }
        }

        if(cctk_bbox[GRHayLMHD_ymax]) {
            for(int k = 0; k < cctk_lsh[2]; k++) {
                for(int i = 0; i < cctk_lsh[0]; i++) {
                    const int src_idx = CCTK_GFINDEX3D(cctkGH, i, jmax - 1, k);
                    const int dst_idx = CCTK_GFINDEX3D(cctkGH, i, jmax, k);

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
            }
        }

        if(cctk_bbox[GRHayLMHD_zmin]) {
            for(int j = 0; j < cctk_lsh[1]; j++) {
                for(int i = 0; i < cctk_lsh[0]; i++) {
                    const int src_idx = CCTK_GFINDEX3D(cctkGH, i, j, kmin + 1);
                    const int dst_idx = CCTK_GFINDEX3D(cctkGH, i, j, kmin);

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
            }
        }

        if(cctk_bbox[GRHayLMHD_zmax]) {
            for(int j = 0; j < cctk_lsh[1]; j++) {
                for(int i = 0; i < cctk_lsh[0]; i++) {
                    const int src_idx = CCTK_GFINDEX3D(cctkGH, i, j, kmax - 1);
                    const int dst_idx = CCTK_GFINDEX3D(cctkGH, i, j, kmax);

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
            }
        }
    } // for(int gzidx = 0; gzidx < cctk_nghostzones[0]; gzidx++)
}
