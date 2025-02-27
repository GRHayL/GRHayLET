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

#define BCLOOP(_bc_dir, _a, _amin, _amax, _b, _bmin, _bmax) \
    if(cctk_bbox[_bc_dir]) {                                \
        for(int _b = _bmin; _b < _bmax; _b++) {             \
            for(int _a = _amin; _a < _amax; _a++) {

#define ENDBCLOOP \
    }             \
    }             \
    }

#define BCLOOP_XY(bc_dir) BCLOOP(bc_dir, i, 0, cctk_lsh[0], j, 0, cctk_lsh[1])
#define BCLOOP_XZ(bc_dir) BCLOOP(bc_dir, i, 0, cctk_lsh[0], k, 0, cctk_lsh[2])
#define BCLOOP_YZ(bc_dir) BCLOOP(bc_dir, j, 0, cctk_lsh[1], k, 0, cctk_lsh[2])

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

        BCLOOP_YZ(GRHayLMHD_xmin)
        {
            const int src_idx = CCTK_GFINDEX3D(cctkGH, imin + 1, j, k);
            const int dst_idx = CCTK_GFINDEX3D(cctkGH, imin, j, k);

            GRHAYLMHD_APPLY_BCS(src_idx, dst_idx);
            if(vx[dst_idx] > 0) {
                vx[dst_idx] = 0;
            }

            GRHayLMHD_Prim2Con_SinglePoint(CCTK_PASS_CTOC, dst_idx);
        }
        ENDBCLOOP

        BCLOOP_YZ(GRHayLMHD_xmax)
        {
            const int src_idx = CCTK_GFINDEX3D(cctkGH, imax - 1, j, k);
            const int dst_idx = CCTK_GFINDEX3D(cctkGH, imax, j, k);

            GRHAYLMHD_APPLY_BCS(src_idx, dst_idx);
            if(vx[dst_idx] < 0) {
                vx[dst_idx] = 0;
            }

            GRHayLMHD_Prim2Con_SinglePoint(CCTK_PASS_CTOC, dst_idx);
        }
        ENDBCLOOP

        BCLOOP_XZ(GRHayLMHD_ymin)
        {
            const int src_idx = CCTK_GFINDEX3D(cctkGH, i, jmin + 1, k);
            const int dst_idx = CCTK_GFINDEX3D(cctkGH, i, jmin, k);

            GRHAYLMHD_APPLY_BCS(src_idx, dst_idx);
            if(vy[dst_idx] > 0) {
                vy[dst_idx] = 0;
            }

            GRHayLMHD_Prim2Con_SinglePoint(CCTK_PASS_CTOC, dst_idx);
        }
        ENDBCLOOP

        BCLOOP_XZ(GRHayLMHD_ymax)
        {
            const int src_idx = CCTK_GFINDEX3D(cctkGH, i, jmax - 1, k);
            const int dst_idx = CCTK_GFINDEX3D(cctkGH, i, jmax, k);

            GRHAYLMHD_APPLY_BCS(src_idx, dst_idx);
            if(vy[dst_idx] < 0) {
                vy[dst_idx] = 0;
            }

            GRHayLMHD_Prim2Con_SinglePoint(CCTK_PASS_CTOC, dst_idx);
        }
        ENDBCLOOP

        BCLOOP_XY(GRHayLMHD_zmin)
        {
            const int src_idx = CCTK_GFINDEX3D(cctkGH, i, j, kmin + 1);
            const int dst_idx = CCTK_GFINDEX3D(cctkGH, i, j, kmin);

            GRHAYLMHD_APPLY_BCS(src_idx, dst_idx);
            if(vz[dst_idx] > 0) {
                vz[dst_idx] = 0;
            }

            GRHayLMHD_Prim2Con_SinglePoint(CCTK_PASS_CTOC, dst_idx);
        }
        ENDBCLOOP

        BCLOOP_XY(GRHayLMHD_zmax)
        {
            const int src_idx = CCTK_GFINDEX3D(cctkGH, i, j, kmax - 1);
            const int dst_idx = CCTK_GFINDEX3D(cctkGH, i, j, kmax);

            GRHAYLMHD_APPLY_BCS(src_idx, dst_idx);
            if(vz[dst_idx] < 0) {
                vz[dst_idx] = 0;
            }

            GRHayLMHD_Prim2Con_SinglePoint(CCTK_PASS_CTOC, dst_idx);
        }
        ENDBCLOOP
    } // for(int gzidx = 0; gzidx < cctk_nghostzones[0]; gzidx++)
}
