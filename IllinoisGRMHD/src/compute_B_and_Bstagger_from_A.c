#include "IllinoisGRMHD.h"

void IllinoisGRMHD_compute_B_and_Bstagger_from_A(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_IllinoisGRMHD_compute_B_and_Bstagger_from_A;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dxi = 1.0/CCTK_DELTA_SPACE(0);
  const CCTK_REAL dyi = 1.0/CCTK_DELTA_SPACE(1);
  const CCTK_REAL dzi = 1.0/CCTK_DELTA_SPACE(2);

//TODO: Do we have to do this every time? IDK how symmetries work in ET
  const CCTK_REAL gridfunc_syms_Ax[3]       = {-1,  1,  Sym_Bz};
  const CCTK_REAL gridfunc_syms_Ay[3]       = { 1, -1,  Sym_Bz};
  const CCTK_REAL gridfunc_syms_Az[3]       = { 1,  1, -Sym_Bz};
  const CCTK_REAL gridfunc_syms_phitilde[3] = { 1,  1,  1};
  IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, x, y, z, Ax,       gridfunc_syms_Ax,       0, 1, 1);
  IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, x, y, z, Ay,       gridfunc_syms_Ay,       1, 0, 1);
  IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, x, y, z, Az,       gridfunc_syms_Az,       1, 1, 0);
  IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, x, y, z, phitilde, gridfunc_syms_phitilde, 1, 1, 1);

#pragma omp parallel for
  for(int k=0; k<cctk_lsh[2]; k++) {
    for(int j=0; j<cctk_lsh[1]; j++) {
      for(int i=0; i<cctk_lsh[0]; i++) {
        // Look Mom, no if() statements!
        const int shiftedim1 = (i-1)*(i!=0); // This way, i=0 yields shiftedim1=0 and shiftedi=1, used below for our COPY boundary condition.
        const int shiftedi   = shiftedim1+1;

        const int shiftedjm1 = (j-1)*(j!=0);
        const int shiftedj   = shiftedjm1+1;

        const int shiftedkm1 = (k-1)*(k!=0);
        const int shiftedk   = shiftedkm1+1;

        int index,indexim1,indexjm1,indexkm1;

        const int actual_index = CCTK_GFINDEX3D(cctkGH, i, j, k);

        // For the lower boundaries, the following applies a "copy"
        //    boundary condition on Bi_stagger where needed.
        //    E.g., Bx_stagger(i,jmin,k) = Bx_stagger(i,jmin+1,k)
        //    We find the copy BC works better than extrapolation.
        // For the upper boundaries, we do the following copy:
        //    E.g., Bx_stagger(imax+1,j,k) = Bx_stagger(imax,j,k)
        // Note that B_stagger is densitized, so we do not divide by psi^6.
        /**************/
        /* Bx_stagger */
        /**************/

        index    = CCTK_GFINDEX3D(cctkGH, i, shiftedj,   shiftedk);
        indexjm1 = CCTK_GFINDEX3D(cctkGH, i, shiftedjm1, shiftedk);
        indexkm1 = CCTK_GFINDEX3D(cctkGH, i, shiftedj,   shiftedkm1);
        // Set Bx_stagger = \partial_y A_z - partial_z A_y
        // "Grid" Ax(i,j,k) is actually Ax(i,j+1/2,k+1/2)
        // "Grid" Ay(i,j,k) is actually Ay(i+1/2,j,k+1/2)
        // "Grid" Az(i,j,k) is actually Az(i+1/2,j+1/2,k)
        // Therefore, the 2nd order derivative \partial_z A_y at (i+1/2,j,k) is:
        //          ["Grid" Ay(i,j,k) - "Grid" Ay(i,j,k-1)]/dZ
        Bx_stagger[actual_index] = (Az[index]-Az[indexjm1])*dyi - (Ay[index]-Ay[indexkm1])*dzi;

        /**************/
        /* By_stagger */
        /**************/

        index    = CCTK_GFINDEX3D(cctkGH, shiftedi,   j, shiftedk);
        indexim1 = CCTK_GFINDEX3D(cctkGH, shiftedim1, j, shiftedk);
        indexkm1 = CCTK_GFINDEX3D(cctkGH, shiftedi,   j, shiftedkm1);
        // Set By_stagger = \partial_z A_x - \partial_x A_z
        By_stagger[actual_index] = (Ax[index]-Ax[indexkm1])*dzi - (Az[index]-Az[indexim1])*dxi;

        /**************/
        /* Bz_stagger */
        /**************/

        index    = CCTK_GFINDEX3D(cctkGH,shiftedi,   shiftedj,   k);
        indexim1 = CCTK_GFINDEX3D(cctkGH,shiftedim1, shiftedj,   k);
        indexjm1 = CCTK_GFINDEX3D(cctkGH,shiftedi,   shiftedjm1, k);
        // Set Bz_stagger = \partial_x A_y - \partial_y A_x
        Bz_stagger[actual_index] = (Ay[index]-Ay[indexim1])*dxi - (Ax[index]-Ax[indexjm1])*dyi;
      }
    }
  }

#pragma omp parallel for
  for(int k=0; k<cctk_lsh[2]; k++) {
    for(int j=0; j<cctk_lsh[1]; j++) {
      for(int i=0; i<cctk_lsh[0]; i++) {
        // Look Mom, no if() statements!
        const int shiftedim1 = (i-1)*(i!=0); // This way, i=0 yields shiftedim1=0 and shiftedi=1, used below for our COPY boundary condition.
        const int shiftedi   = shiftedim1+1;

        const int shiftedjm1 = (j-1)*(j!=0);
        const int shiftedj   = shiftedjm1+1;

        const int shiftedkm1 = (k-1)*(k!=0);
        const int shiftedk   = shiftedkm1+1;

        int index,indexim1,indexjm1,indexkm1;

        const int actual_index = CCTK_GFINDEX3D(cctkGH, i, j, k);

        // Unlike B_stagger, B_center is undensitized. Therefore, we
        // need to divide by psi^6 = \sqrt{gamma}
        const CCTK_REAL sqrt_detgamma = sqrt(fabs(
                gxx[actual_index] * (gyy[actual_index]*gzz[actual_index] - gyz[actual_index]*gyz[actual_index])
              + gxy[actual_index] * (gyz[actual_index]*gxz[actual_index] - gxy[actual_index]*gzz[actual_index])
              + gxz[actual_index] * (gxy[actual_index]*gyz[actual_index] - gyy[actual_index]*gxz[actual_index])));

        // For the lower boundaries, the following applies a "copy"
        //    boundary condition on Bi and Bi_stagger where needed.
        //    E.g., Bx_center(imin,j,k) = Bx_center(imin+1,j,k)
        //    We find the copy BC works better than extrapolation.
        /*************/
        /* Bx_center */
        /*************/
        index = CCTK_GFINDEX3D(cctkGH, shiftedi, j, k);
        indexim1 = CCTK_GFINDEX3D(cctkGH, shiftedim1, j, k);
        // Set Bx_center = 0.5 ( Bx_stagger + Bx_stagger_im1 )
        // "Grid" Bx_stagger(i,j,k) is actually Bx_stagger(i+1/2,j,k)
        Bx_center[actual_index] = 0.5 * ( Bx_stagger[index] + Bx_stagger[indexim1] )/sqrt_detgamma;

        /*************/
        /* By_center */
        /*************/
        index = CCTK_GFINDEX3D(cctkGH, i, shiftedj, k);
        indexjm1 = CCTK_GFINDEX3D(cctkGH, i, shiftedjm1, k);
        // Set By_center = 0.5 ( By_stagger + By_stagger_im1 )
        // "Grid" By_stagger(i,j,k) is actually By_stagger(i,j+1/2,k)
        By_center[actual_index] = 0.5 * ( By_stagger[index] + By_stagger[indexjm1] )/sqrt_detgamma;

        /*************/
        /* Bz_center */
        /*************/
        index = CCTK_GFINDEX3D(cctkGH, i, j, shiftedk);
        indexkm1 = CCTK_GFINDEX3D(cctkGH, i, j, shiftedkm1);
        // Set Bz_center = 0.5 ( Bz_stagger + Bz_stagger_im1 )
        // "Grid" Bz_stagger(i,j,k) is actually Bz_stagger(i,j+1/2,k)
        Bz_center[actual_index] = 0.5 * ( Bz_stagger[index] + Bz_stagger[indexkm1] )/sqrt_detgamma;
      }
    }
  }

  // Finish up by setting symmetry ghostzones on Bx, By, Bz, and their staggered variants.
  const CCTK_REAL gridfunc_syms_Bx[3] = {-1,  1, -Sym_Bz};
  const CCTK_REAL gridfunc_syms_By[3] = { 1, -1, -Sym_Bz};
  const CCTK_REAL gridfunc_syms_Bz[3] = { 1,  1,  Sym_Bz};
  IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, x, y, z,  Bx_center,  gridfunc_syms_Bx, 0, 0, 0);
  IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, x, y, z,  Bx_stagger, gridfunc_syms_Bx, 1, 0, 0);
  IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, x, y, z,  By_center,  gridfunc_syms_By, 0, 0, 0);
  IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, x, y, z,  By_stagger, gridfunc_syms_By, 0, 1, 0);
  IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, x, y, z,  Bz_center,  gridfunc_syms_Bz, 0, 0, 0);
  IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, x, y, z,  Bz_stagger, gridfunc_syms_Bz, 0, 0, 1);
}
