#include "GRHayLMHD.h"

#define GRHAYLMHD_APPLY_EM_BCS(_ijk, _ijk_pm1, _ijk_pm2)               \
    A_x_tilde[_ijk] = 2.0 * A_x_tilde[_ijk_pm1] - A_x_tilde[_ijk_pm2]; \
    A_y_tilde[_ijk] = 2.0 * A_y_tilde[_ijk_pm1] - A_y_tilde[_ijk_pm2]; \
    A_z_tilde[_ijk] = 2.0 * A_z_tilde[_ijk_pm1] - A_z_tilde[_ijk_pm2]; \
    Phi_tilde[_ijk] = 2.0 * Phi_tilde[_ijk_pm1] - Phi_tilde[_ijk_pm2];

void GRHayLMHD_Outer_Boundaries_EM_FourPotential(CCTK_ARGUMENTS)
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
                    const int ijk   = CCTK_GFINDEX3D(cctkGH, imin, j, k);
                    const int ip1jk = CCTK_GFINDEX3D(cctkGH, imin + 1, j, k);
                    const int ip2jk = CCTK_GFINDEX3D(cctkGH, imin + 2, j, k);
                    GRHAYLMHD_APPLY_EM_BCS(ijk, ip1jk, ip2jk);
                }
            }
        }

        if(cctk_bbox[GRHayLMHD_xmax]) {
            for(int k = 0; k < cctk_lsh[2]; k++) {
                for(int j = 0; j < cctk_lsh[1]; j++) {
                    const int ijk   = CCTK_GFINDEX3D(cctkGH, imax, j, k);
                    const int im1jk = CCTK_GFINDEX3D(cctkGH, imax - 1, j, k);
                    const int im2jk = CCTK_GFINDEX3D(cctkGH, imax - 2, j, k);
                    GRHAYLMHD_APPLY_EM_BCS(ijk, im1jk, im2jk);
                }
            }
        }

        if(cctk_bbox[GRHayLMHD_ymin]) {
            for(int k = 0; k < cctk_lsh[2]; k++) {
                for(int i = 0; i < cctk_lsh[0]; i++) {
                    const int ijk   = CCTK_GFINDEX3D(cctkGH, i, jmin, k);
                    const int ijp1k = CCTK_GFINDEX3D(cctkGH, i, jmin + 1, k);
                    const int ijp2k = CCTK_GFINDEX3D(cctkGH, i, jmin + 2, k);
                    GRHAYLMHD_APPLY_EM_BCS(ijk, ijp1k, ijp2k);
                }
            }
        }

        if(cctk_bbox[GRHayLMHD_ymax]) {
            for(int k = 0; k < cctk_lsh[2]; k++) {
                for(int i = 0; i < cctk_lsh[0]; i++) {
                    const int ijk   = CCTK_GFINDEX3D(cctkGH, i, jmax, k);
                    const int ijm1k = CCTK_GFINDEX3D(cctkGH, i, jmax - 1, k);
                    const int ijm2k = CCTK_GFINDEX3D(cctkGH, i, jmax - 2, k);
                    GRHAYLMHD_APPLY_EM_BCS(ijk, ijm1k, ijm2k);
                }
            }
        }

        if(cctk_bbox[GRHayLMHD_zmin]) {
            for(int j = 0; j < cctk_lsh[1]; j++) {
                for(int i = 0; i < cctk_lsh[0]; i++) {
                    const int ijk   = CCTK_GFINDEX3D(cctkGH, i, j, kmin);
                    const int ijkp1 = CCTK_GFINDEX3D(cctkGH, i, j, kmin + 1);
                    const int ijkp2 = CCTK_GFINDEX3D(cctkGH, i, j, kmin + 2);
                    GRHAYLMHD_APPLY_EM_BCS(ijk, ijkp1, ijkp2);
                }
            }
        }

        if(cctk_bbox[GRHayLMHD_zmax]) {
            for(int j = 0; j < cctk_lsh[1]; j++) {
                for(int i = 0; i < cctk_lsh[0]; i++) {
                    const int ijk   = CCTK_GFINDEX3D(cctkGH, i, j, kmax);
                    const int ijkm1 = CCTK_GFINDEX3D(cctkGH, i, j, kmax - 1);
                    const int ijkm2 = CCTK_GFINDEX3D(cctkGH, i, j, kmax - 2);
                    GRHAYLMHD_APPLY_EM_BCS(ijk, ijkm1, ijkm2);
                }
            }
        }
    }
}
