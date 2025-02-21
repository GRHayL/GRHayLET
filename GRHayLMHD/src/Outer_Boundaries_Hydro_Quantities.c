#include "GRHayLMHD.h"

#define GRHAYLMHD_APPLY_BCS(_src, _dst)    \
    rho[_dst]         = rho[_src];         \
    press[_dst]       = press[_src];       \
    eps[_dst]         = eps[_src];         \
    Y_e[_dst]         = Y_e[_src];         \
    temperature[_dst] = temperature[_src]; \
    entropy[_dst]     = entropy[_src];     \
    vx[_dst]          = vx[_src];          \
    vy[_dst]          = vy[_src];          \
    vz[_dst]          = vz[_src];

static void Prim2Con(CCTK_ARGUMENTS, const int i, const int j, const int k)
{
    DECLARE_CCTK_PARAMETERS;
    const int ijk  = CCTK_GFINDEX3D(cctkGH, i, j, k);
    const int ijkx = CCTK_GFINDEX4D(cctkGH, i, j, k, 0);
    const int ijky = CCTK_GFINDEX4D(cctkGH, i, j, k, 1);
    const int ijkz = CCTK_GFINDEX4D(cctkGH, i, j, k, 2);
    GRHayLMHD_Prim2Con_SinglePoint(CCTK_PASS_CTOC, ijk, ijkx, ijky, ijkz);
}

void GRHayLMHD_Outer_Boundaries_Hydro_Quantities(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // TODO: add BC options
    if(cctk_iteration == 0 || GetRefinementLevel(cctkGH)) {
        return;
    }

    const int lower[3] = {
        cctk_nghostzones[0] - 1,
        cctk_nghostzones[1] - 1,
        cctk_nghostzones[2] - 1,
    };
    const int upper[3] = {
        cctk_lsh[0] - cctk_nghostzones[0],
        cctk_lsh[1] - cctk_nghostzones[1],
        cctk_lsh[2] - cctk_nghostzones[2],
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

                    GRHAYLMHD_APPLY_BCS(src_idx, dst_idx);
                    if(vx[dst_idx] > 0) {
                        vx[dst_idx] = 0;
                    }

                    Prim2Con(CCTK_PASS_CTOC, imin, j, k);
                }
            }
        }

        if(cctk_bbox[GRHayLMHD_xmax]) {
            for(int k = 0; k < cctk_lsh[2]; k++) {
                for(int j = 0; j < cctk_lsh[1]; j++) {
                    const int src_idx = CCTK_GFINDEX3D(cctkGH, imax - 1, j, k);
                    const int dst_idx = CCTK_GFINDEX3D(cctkGH, imax, j, k);

                    GRHAYLMHD_APPLY_BCS(src_idx, dst_idx);
                    if(vx[dst_idx] < 0) {
                        vx[dst_idx] = 0;
                    }

                    Prim2Con(CCTK_PASS_CTOC, imax, j, k);
                }
            }
        }

        if(cctk_bbox[GRHayLMHD_ymin]) {
            for(int k = 0; k < cctk_lsh[2]; k++) {
                for(int i = 0; i < cctk_lsh[0]; i++) {
                    const int src_idx = CCTK_GFINDEX3D(cctkGH, i, jmin + 1, k);
                    const int dst_idx = CCTK_GFINDEX3D(cctkGH, i, jmin, k);

                    GRHAYLMHD_APPLY_BCS(src_idx, dst_idx);
                    if(vy[dst_idx] > 0) {
                        vy[dst_idx] = 0;
                    }

                    Prim2Con(CCTK_PASS_CTOC, i, jmin, k);
                }
            }
        }

        if(cctk_bbox[GRHayLMHD_ymax]) {
            for(int k = 0; k < cctk_lsh[2]; k++) {
                for(int i = 0; i < cctk_lsh[0]; i++) {
                    const int src_idx = CCTK_GFINDEX3D(cctkGH, i, jmax - 1, k);
                    const int dst_idx = CCTK_GFINDEX3D(cctkGH, i, jmax, k);

                    GRHAYLMHD_APPLY_BCS(src_idx, dst_idx);
                    if(vy[dst_idx] < 0) {
                        vy[dst_idx] = 0;
                    }

                    Prim2Con(CCTK_PASS_CTOC, i, jmax, k);
                }
            }
        }

        if(cctk_bbox[GRHayLMHD_zmin]) {
            for(int j = 0; j < cctk_lsh[1]; j++) {
                for(int i = 0; i < cctk_lsh[0]; i++) {
                    const int src_idx = CCTK_GFINDEX3D(cctkGH, i, j, kmin + 1);
                    const int dst_idx = CCTK_GFINDEX3D(cctkGH, i, j, kmin);

                    GRHAYLMHD_APPLY_BCS(src_idx, dst_idx);
                    if(vz[dst_idx] > 0) {
                        vz[dst_idx] = 0;
                    }

                    Prim2Con(CCTK_PASS_CTOC, i, j, kmin);
                }
            }
        }

        if(cctk_bbox[GRHayLMHD_zmax]) {
            for(int j = 0; j < cctk_lsh[1]; j++) {
                for(int i = 0; i < cctk_lsh[0]; i++) {
                    const int src_idx = CCTK_GFINDEX3D(cctkGH, i, j, kmax - 1);
                    const int dst_idx = CCTK_GFINDEX3D(cctkGH, i, j, kmax);

                    GRHAYLMHD_APPLY_BCS(src_idx, dst_idx);
                    if(vz[dst_idx] < 0) {
                        vz[dst_idx] = 0;
                    }

                    Prim2Con(CCTK_PASS_CTOC, i, j, kmax);
                }
            }
        }
    } // for(int gzidx = 0; gzidx < cctk_nghostzones[0]; gzidx++)
}
