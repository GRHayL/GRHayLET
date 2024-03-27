#ifndef ILLINOISGRMHDX_H_
#define ILLINOISGRMHDX_H_

#include "loop_device.hxx"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "GRHayLib.h"

enum recon_indices{
      BX_STAGGER, BY_STAGGER, BZ_STAGGER,
      VXR, VYR, VZR, VXL,VYL, VZL, MAXNUMVARS};

//template <int flux_dir>
//void IllinoisGRMHDX_reconstruction_loop(
//      const cGH *restrict cctkGH,
//      const int num_B,
//      const int num_vel,
//      const int B_indices[3],
//      const int vel_indices[3],
//      Loop::GF3D2<const CCTK_REAL> v_flux_dir,
//      const std::array<Loop::GF3D2<const CCTK_REAL>, 3> Bstag,
//      const std::array<Loop::GF3D2<CCTK_REAL>, 3> Bstag_r,
//      const std::array<Loop::GF3D2<CCTK_REAL>, 3> Bstag_l,
//      const std::array<Loop::GF3D2<CCTK_REAL>, 6> vrl,
//      const std::array<Loop::GF3D2<CCTK_REAL>, 6> vrl_r,
//      const std::array<Loop::GF3D2<CCTK_REAL>, 6> vrl_l);

template <int flux_dir>
void IllinoisGRMHDX_A_flux_rhs(
      const cGH *restrict cctkGH,
      Loop::GF3D2<CCTK_REAL> A_rhs);

// The inner two points of the interpolation function use
// the value of A_in, and the outer two points use A_out.
#define A_out -0.0625
#define A_in  0.5625
//Interpolates to the -1/2 face of point Var
#define COMPUTE_FCVAL(Varm2,Varm1,Var,Varp1) (A_out*(Varm2) + A_in*(Varm1) + A_in*(Var) + A_out*(Varp1))

// Computes 4th-order derivative
#define B_out -1.0/12.0
#define B_in  2.0/3.0
#define COMPUTE_DERIV(Varm2,Varm1,Varp1,Varp2) (B_in*(Varp1 - Varm1) + B_out*(Varp2 - Varm2))

extern "C" void IllinoisGRMHDX_interpolate_metric_to_face(
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

extern "C" void IllinoisGRMHDX_compute_metric_derivs(
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

#endif // ILLINOISGRMHDX_H_
