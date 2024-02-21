#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "GRHayLib.h"


// Computes 4th-order derivative
#define B_out -1.0/12.0
#define B_in  2.0/3.0
#define COMPUTE_DERIV(Varm2,Varm1,Varp1,Varp2) (B_in*(Varp1 - Varm1) + B_out*(Varp2 - Varm2))


static void GRHayLHD_compute_metric_derivs(
      const cGH *cctkGH,
      const int i, const int j, const int k,
      const int flux_dir,
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
      ghl_metric_quantities *restrict metric_derivs) {

  const int xdir = (flux_dir == 0);
  const int ydir = (flux_dir == 1);
  const int zdir = (flux_dir == 2);

  const int indm2  = CCTK_GFINDEX3D(cctkGH, i-2*xdir, j-2*ydir, k-2*zdir);
  const int indm1  = CCTK_GFINDEX3D(cctkGH, i-xdir,   j-ydir,   k-zdir);
  const int indp1  = CCTK_GFINDEX3D(cctkGH, i+xdir,   j+ydir,   k+zdir);
  const int indp2  = CCTK_GFINDEX3D(cctkGH, i+2*xdir, j+2*ydir, k+2*zdir);

  const CCTK_REAL d_lapse = dxi*COMPUTE_DERIV(lapse[indm2], lapse[indm1], lapse[indp1], lapse[indp2]);
  const CCTK_REAL d_betax = dxi*COMPUTE_DERIV(betax[indm2], betax[indm1], betax[indp1], betax[indp2]);
  const CCTK_REAL d_betay = dxi*COMPUTE_DERIV(betay[indm2], betay[indm1], betay[indp1], betay[indp2]);
  const CCTK_REAL d_betaz = dxi*COMPUTE_DERIV(betaz[indm2], betaz[indm1], betaz[indp1], betaz[indp2]);

  const CCTK_REAL d_gxx = dxi*COMPUTE_DERIV(gxx[indm2], gxx[indm1], gxx[indp1], gxx[indp2]);
  const CCTK_REAL d_gxy = dxi*COMPUTE_DERIV(gxy[indm2], gxy[indm1], gxy[indp1], gxy[indp2]);
  const CCTK_REAL d_gxz = dxi*COMPUTE_DERIV(gxz[indm2], gxz[indm1], gxz[indp1], gxz[indp2]);
  const CCTK_REAL d_gyy = dxi*COMPUTE_DERIV(gyy[indm2], gyy[indm1], gyy[indp1], gyy[indp2]);
  const CCTK_REAL d_gyz = dxi*COMPUTE_DERIV(gyz[indm2], gyz[indm1], gyz[indp1], gyz[indp2]);
  const CCTK_REAL d_gzz = dxi*COMPUTE_DERIV(gzz[indm2], gzz[indm1], gzz[indp1], gzz[indp2]);

  ghl_initialize_metric(
        d_lapse,
        d_betax, d_betay, d_betaz,
        d_gxx, d_gxy, d_gxz,
        d_gyy, d_gyz, d_gzz,
        metric_derivs);
}


void GRHayLHD_tabulated_evaluate_sources_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const int imin = cctk_nghostzones[0];
  const int jmin = cctk_nghostzones[1];
  const int kmin = cctk_nghostzones[2];
  const int imax = cctk_lsh[0] - cctk_nghostzones[0];
  const int jmax = cctk_lsh[1] - cctk_nghostzones[1];
  const int kmax = cctk_lsh[2] - cctk_nghostzones[2];
  const CCTK_REAL dxi = 1.0/CCTK_DELTA_SPACE(0);
  const CCTK_REAL dyi = 1.0/CCTK_DELTA_SPACE(1);
  const CCTK_REAL dzi = 1.0/CCTK_DELTA_SPACE(2);

#pragma omp parallel for
  for(int k=kmin; k<kmax; k++) {
    for(int j=jmin; j<jmax; j++) {
      for(int i=imin; i<imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j ,k);

        // These variables have no source terms
        rho_star_rhs[index] = 0.0;
        Ye_star_rhs[index]  = 0.0;

        ghl_metric_quantities ADM_metric;
        ghl_initialize_metric(
              alp[index],
              betax[index], betay[index], betaz[index],
              gxx[index], gxy[index], gxz[index],
              gyy[index], gyz[index], gzz[index],
              &ADM_metric);

        ghl_extrinsic_curvature curv;
        ghl_initialize_extrinsic_curvature(
              kxx[index], kxy[index], kxz[index],
              kyy[index], kyz[index], kzz[index],
              &curv);

        ghl_primitive_quantities prims;
        prims.BU[0] = prims.BU[1] = prims.BU[2] = 0.0;
        prims.rho   = rho[index];
        prims.press = press[index];
        prims.vU[0] = vx[index];
        prims.vU[1] = vy[index];
        prims.vU[2] = vz[index];
        prims.Y_e   = Y_e[index];
        prims.temperature = temperature[index];

        const int speed_limited CCTK_ATTRIBUTE_UNUSED = ghl_limit_v_and_compute_u0(ghl_params, &ADM_metric, &prims);

        ghl_metric_quantities ADM_metric_derivs_x;
        GRHayLHD_compute_metric_derivs(
              cctkGH, i, j, k,
              0, dxi, alp,
              betax, betay, betaz,
              gxx, gxy, gxz,
              gyy, gyz, gzz,
              &ADM_metric_derivs_x);

        ghl_metric_quantities ADM_metric_derivs_y;
        GRHayLHD_compute_metric_derivs(
              cctkGH, i, j, k,
              1, dyi, alp,
              betax, betay, betaz,
              gxx, gxy, gxz,
              gyy, gyz, gzz,
              &ADM_metric_derivs_y);

        ghl_metric_quantities ADM_metric_derivs_z;
        GRHayLHD_compute_metric_derivs(
              cctkGH, i, j, k,
              2, dzi, alp,
              betax, betay, betaz,
              gxx, gxy, gxz,
              gyy, gyz, gzz,
              &ADM_metric_derivs_z);

        ghl_conservative_quantities cons_source;
        ghl_calculate_source_terms(
              ghl_eos, &prims, &ADM_metric,
              &ADM_metric_derivs_x,
              &ADM_metric_derivs_y,
              &ADM_metric_derivs_z,
              &curv, &cons_source);

        tau_rhs    [index] = cons_source.tau;
        Stildex_rhs[index] = cons_source.SD[0];
        Stildey_rhs[index] = cons_source.SD[1];
        Stildez_rhs[index] = cons_source.SD[2];
      }
    }
  }
}

