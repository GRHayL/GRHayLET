#include "GRHayLMHD.h"

// This function is simple, but requires bookkeeping. When we access a point
// (i, j, k), in reality we are accessing staggered fields at the following
// locations:
//
// Staggered magnetic (vector) potential:
// A_x: (i      , j + 1/2, k + 1/2)
// A_y: (i + 1/2, j      , k + 1/2)
// A_z: (i + 1/2, j + 1/2, k      )
//
// Staggered magnetic fields:
// Bx : (i + 1/2, j      , k      )
// By : (i      , j + 1/2, k      )
// Bz : (i      , j      , k + 1/2)
//
// These are chosen to make our lives as easy as possible. For example, we
// compute the derivative of A_x along the y-direction using:
//
// partial_y A_x = ( A_x(i, j + 1/2, k + 1/2) - A_x(i, j - 1/2, k + 1/2) ) / dy ,
//
// Similarly, we compute the derivative of A_y along the x-direction using:
//
// partial_x A_y = ( A_y(i + 1/2, j, k + 1/2) - A_y(i - 1/2, j, k + 1/2) ) / dx ,
//
// Note that both derivatives are evaluated at the point (i, j, k + 1/2),
// exactly what we need to compute Bz_stagger:
//
// Bz_stagger = partial_y A_x - partial_x A_y .
//
// Once we have the staggered magnetic fields, then we just average them
// to obtain the magnetic fields at (i, j, k), e.g.,
//
// Bx(i, j, k) = 0.5 * ( Bx_stagger(i + 1/2, j, k) + Bx_stagger(i - 1/2, j, k) ) .
void GRHayLMHD_Compute_Magnetic_Fields(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    const CCTK_REAL inv_dx = 1.0 / CCTK_DELTA_SPACE(0);
    const CCTK_REAL inv_dy = 1.0 / CCTK_DELTA_SPACE(1);
    const CCTK_REAL inv_dz = 1.0 / CCTK_DELTA_SPACE(2);

    const int imin = 0, imax = cctk_lsh[0];
    const int jmin = 0, jmax = cctk_lsh[1];
    const int kmin = 0, kmax = cctk_lsh[2];

    OMPLOOP3D(imin, imax, jmin, jmax, kmin, kmax)
    {
        const int im1 = (i != 0) * (i - 1);
        const int jm1 = (j != 0) * (j - 1);
        const int km1 = (k != 0) * (k - 1);

        const int im1jk = CCTK_GFINDEX3D(cctkGH, im1, j, k);
        const int ijm1k = CCTK_GFINDEX3D(cctkGH, i, jm1, k);
        const int ijkm1 = CCTK_GFINDEX3D(cctkGH, i, j, km1);

        const int ip1jk = CCTK_GFINDEX3D(cctkGH, im1 + 1, j, k);
        const int ijp1k = CCTK_GFINDEX3D(cctkGH, i, jm1 + 1, k);
        const int ijkp1 = CCTK_GFINDEX3D(cctkGH, i, j, km1 + 1);

        const CCTK_REAL partial_y_A_x_tilde = (A_x_tilde[ijp1k] - A_x_tilde[ijm1k]) * inv_dy;
        const CCTK_REAL partial_z_A_x_tilde = (A_x_tilde[ijkp1] - A_x_tilde[ijkm1]) * inv_dz;
        const CCTK_REAL partial_x_A_y_tilde = (A_y_tilde[ip1jk] - A_y_tilde[im1jk]) * inv_dx;
        const CCTK_REAL partial_z_A_y_tilde = (A_y_tilde[ijkp1] - A_y_tilde[ijkm1]) * inv_dz;
        const CCTK_REAL partial_x_A_z_tilde = (A_z_tilde[ip1jk] - A_z_tilde[im1jk]) * inv_dx;
        const CCTK_REAL partial_y_A_z_tilde = (A_z_tilde[ijp1k] - A_z_tilde[ijm1k]) * inv_dy;

        Bx_tilde_stagger[ijk] = partial_y_A_z_tilde - partial_z_A_y_tilde;
        By_tilde_stagger[ijk] = partial_z_A_x_tilde - partial_x_A_z_tilde;
        Bz_tilde_stagger[ijk] = partial_x_A_y_tilde - partial_y_A_x_tilde;
    }
    ENDLOOP3D

    OMPLOOP3D(imin, imax, jmin, jmax, kmin, kmax)
    {
        const int im1 = (i != 0) * (i - 1);
        const int jm1 = (j != 0) * (j - 1);
        const int km1 = (k != 0) * (k - 1);

        const int im1jk = CCTK_GFINDEX3D(cctkGH, im1, j, k);
        const int ijm1k = CCTK_GFINDEX3D(cctkGH, i, jm1, k);
        const int ijkm1 = CCTK_GFINDEX3D(cctkGH, i, j, km1);

        const int ip1jk = CCTK_GFINDEX3D(cctkGH, im1 + 1, j, k);
        const int ijp1k = CCTK_GFINDEX3D(cctkGH, i, jm1 + 1, k);
        const int ijkp1 = CCTK_GFINDEX3D(cctkGH, i, j, km1 + 1);

        Bvec[ijkx] = 0.5 * ( Bx_tilde_stagger[ip1jk] + Bx_tilde_stagger[im1jk] );
        Bvec[ijky] = 0.5 * ( By_tilde_stagger[ijp1k] + By_tilde_stagger[ijm1k] );
        Bvec[ijkz] = 0.5 * ( Bz_tilde_stagger[ijkp1] + Bz_tilde_stagger[ijkm1] );
    }
    ENDLOOP3D
}
