#ifndef GRHAYLHD_H_
#define GRHAYLHD_H_

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "GRHayLib.h"

// This is used to perturb data for testing
#define one_plus_pert(perturb) (1 + (perturb*(double)rand() / RAND_MAX))

// The inner two points of the interpolation function use
// the value of A_in, and the outer two points use A_out.
#define A_out -0.0625
#define A_in  0.5625
// Interpolates to the +1/2 face of point Var
#define COMPUTE_FCVAL(Varm1,Var,Varp1,Varp2) (A_out*(Varm1) + A_in*(Var) + A_in*(Varp1) + A_out*(Varp2))

/*
   Computes derivative factor face(+1/2) - face(-1/2):
   Let A = A_out, B = A_in. Let Var at index i be f[i]. Then,
   dx*deriv = Af[-1] + Bf[0] + Bf[1] + Af[2] - (Af[-2] + Bf[-1] + Bf[0] + Af[1])
            = Af[-1] - Bf[-1] + Bf[1] - Af[1] + Af[2] - Af[-2]
            = (B-A)(f[1] - f[-1]) + A(f[2] - f[-2])
*/ 
#define COMPUTE_DERIV(Varm2,Varm1,Varp1,Varp2) ((A_in - A_out)*(Varp1 - Varm1) + A_out*(Varp2 - Varm2))

void GRHayLHD_interpolate_metric_to_face(
      const cGH *cctkGH,
      const int i, const int j, const int k,
      const int flux_dirn,
      const CCTK_REAL *restrict lapse,
      const CCTK_REAL *restrict betax,
      const CCTK_REAL *restrict betay,
      const CCTK_REAL *restrict betaz,
      const CCTK_REAL *restrict gxx,
      const CCTK_REAL *restrict gxy,
      const CCTK_REAL *restrict gxz,
      const CCTK_REAL *restrict gyy,
      const CCTK_REAL *restrict gyz,
      const CCTK_REAL *restrict gzz,
      metric_quantities *restrict metric);

void GRHayLHD_compute_metric_derivs(
      const cGH *cctkGH,
      const int i, const int j, const int k,
      const int flux_dirn,
      const CCTK_REAL dxi,
      const CCTK_REAL *restrict lapse,
      const CCTK_REAL *restrict betax,
      const CCTK_REAL *restrict betay,
      const CCTK_REAL *restrict betaz,
      const CCTK_REAL *restrict gxx,
      const CCTK_REAL *restrict gxy,
      const CCTK_REAL *restrict gxz,
      const CCTK_REAL *restrict gyy,
      const CCTK_REAL *restrict gyz,
      const CCTK_REAL *restrict gzz,
      metric_quantities *restrict metric_derivs);

#endif // GRHAYLHD_H_
