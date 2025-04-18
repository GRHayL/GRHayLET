#include "IllinoisGRMHD.h"

void IllinoisGRMHD_A_flux_rhs(
      const cGH *restrict cctkGH,
      const int flux_dir,
      /*const*/ CCTK_REAL **in_prims_r,
      /*const*/ CCTK_REAL **in_prims_l,
      /*const*/ CCTK_REAL **cmin,
      /*const*/ CCTK_REAL **cmax,
      CCTK_REAL *restrict A_rhs) {

  const int xdir = (flux_dir==0);
  const int ydir = (flux_dir==1);
  const int zdir = (flux_dir==2);

  // These are used to determine which components of v,
  // B, cmin, and cmax are used in the computation.
  const int dir1_offset = (flux_dir+1)%3, dir2_offset = (flux_dir+2)%3;

  // This offsets the index by +1 in the perpendicular directions
  const int v_offset[3] = { !xdir, !ydir, !zdir };

  // This offsets the index by +1 in the permuted direction (x<-y<-z)
  const int B1_offset[3] = { ydir, zdir, xdir };

  // This offsets the index by +1 in the permuted direction (x->y->z)
  const int B2_offset[3] = { zdir, xdir, ydir };

  const int imin = cctkGH->cctk_nghostzones[0];
  const int imax = cctkGH->cctk_lsh[0]-cctkGH->cctk_nghostzones[0];
  const int jmin = cctkGH->cctk_nghostzones[1];
  const int jmax = cctkGH->cctk_lsh[1]-cctkGH->cctk_nghostzones[1];
  const int kmin = cctkGH->cctk_nghostzones[2];
  const int kmax = cctkGH->cctk_lsh[2]-cctkGH->cctk_nghostzones[2];

#pragma omp parallel for
  for(int k=kmin; k<kmax; k++) {
    for(int j=jmin; j<jmax; j++) {
      for(int i=imin; i<imax; i++) {
        const int index    = CCTK_GFINDEX3D(cctkGH,i,j,k);
        const int index_v  = CCTK_GFINDEX3D(cctkGH,i+v_offset[0], j+v_offset[1], k+v_offset[2]);
        const int index_B1 = CCTK_GFINDEX3D(cctkGH,i+B1_offset[0],j+B1_offset[1],k+B1_offset[2]);
        const int index_B2 = CCTK_GFINDEX3D(cctkGH,i+B2_offset[0],j+B2_offset[1],k+B2_offset[2]);

        ghl_HLL_vars vars;

        vars.v1rr = in_prims_r[VXR+dir1_offset][index_v];
        vars.v1rl = in_prims_l[VXR+dir1_offset][index_v];
        vars.v1lr = in_prims_r[VXL+dir1_offset][index_v];
        vars.v1ll = in_prims_l[VXL+dir1_offset][index_v];

        vars.v2rr = in_prims_r[VXR+dir2_offset][index_v];
        vars.v2rl = in_prims_l[VXR+dir2_offset][index_v];
        vars.v2lr = in_prims_r[VXL+dir2_offset][index_v];
        vars.v2ll = in_prims_l[VXL+dir2_offset][index_v];

        vars.B1r = in_prims_r[BX_STAGGER+dir1_offset][index_B1];
        vars.B1l = in_prims_l[BX_STAGGER+dir1_offset][index_B1];

        vars.B2r = in_prims_r[BX_STAGGER+dir2_offset][index_B2];
        vars.B2l = in_prims_l[BX_STAGGER+dir2_offset][index_B2];

        /*
          Note that cmax/cmin (\alpha^{\pm}  as defined in Del Zanna et al) is at a slightly DIFFERENT
          point (e.g., (i+1/2,j,k) instead of (i+1/2,j+1/2,k) for F3). Yuk Tung Liu discussed this point
          with M. Shibata, who found that the effect is negligible.
        */
        vars.c1_min = cmin[dir1_offset][index_B2];
        vars.c1_max = cmax[dir1_offset][index_B2];
        vars.c2_min = cmin[dir2_offset][index_B1];
        vars.c2_max = cmax[dir2_offset][index_B1];

        A_rhs[index] = ghl_HLL_flux_with_Btilde(&vars);
      }
    }
  }
}
