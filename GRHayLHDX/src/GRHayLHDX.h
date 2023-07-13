#ifndef GRHAYLHDX_H_
#define GRHAYLHDX_H_

#include "loop_device.hxx"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "GRHayLib.h"

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

extern "C" void GRHayLHDX_interpolate_metric_to_face(
      const Loop::GF3D2index indm1,
      const Loop::GF3D2index index,
      const Loop::GF3D2index indp1,
      const Loop::GF3D2index indp2,
      const Loop::GF3D2<const CCTK_REAL> lapse,
      const Loop::GF3D2<const CCTK_REAL> betax,
      const Loop::GF3D2<const CCTK_REAL> betay,
      const Loop::GF3D2<const CCTK_REAL> betaz,
      const Loop::GF3D2<const CCTK_REAL> gxx,
      const Loop::GF3D2<const CCTK_REAL> gxy,
      const Loop::GF3D2<const CCTK_REAL> gxz,
      const Loop::GF3D2<const CCTK_REAL> gyy,
      const Loop::GF3D2<const CCTK_REAL> gyz,
      const Loop::GF3D2<const CCTK_REAL> gzz,
      ghl_metric_quantities *restrict metric);

extern "C" void GRHayLHDX_compute_metric_derivs(
      const CCTK_REAL dxi,
      const Loop::GF3D2index indm2,
      const Loop::GF3D2index indm1,
      const Loop::GF3D2index indp1,
      const Loop::GF3D2index indp2,
      const Loop::GF3D2<const CCTK_REAL> lapse,
      const Loop::GF3D2<const CCTK_REAL> betax,
      const Loop::GF3D2<const CCTK_REAL> betay,
      const Loop::GF3D2<const CCTK_REAL> betaz,
      const Loop::GF3D2<const CCTK_REAL> gxx,
      const Loop::GF3D2<const CCTK_REAL> gxy,
      const Loop::GF3D2<const CCTK_REAL> gxz,
      const Loop::GF3D2<const CCTK_REAL> gyy,
      const Loop::GF3D2<const CCTK_REAL> gyz,
      const Loop::GF3D2<const CCTK_REAL> gzz,
      ghl_metric_quantities *restrict metric_derivs);

#endif // GRHAYLHDX_H_
