#include "GRHayLMHD.h"

void GRHayLMHD_tabulated_evaluate_sources_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLMHD_evaluate_sources_rhs;
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
        prims.rho   = rho_b[index];
        prims.press = pressure[index];
        prims.vU[0] = vx[index];
        prims.vU[1] = vy[index];
        prims.vU[2] = vz[index];
        prims.BU[0] = Bx_center[index];
        prims.BU[1] = By_center[index];
        prims.BU[2] = Bz_center[index];
        prims.Y_e = Y_e[index];
        prims.temperature = temperature[index];
        const int speed_limited CCTK_ATTRIBUTE_UNUSED = ghl_limit_v_and_compute_u0(ghl_params, &ADM_metric, &prims);

        ghl_metric_quantities ADM_metric_derivs_x;
        GRHayLMHD_compute_metric_derivs(
              cctkGH, i, j, k,
              0, dxi, alp,
              betax, betay, betaz,
              gxx, gxy, gxz,
              gyy, gyz, gzz,
              &ADM_metric_derivs_x);

        ghl_metric_quantities ADM_metric_derivs_y;
        GRHayLMHD_compute_metric_derivs(
              cctkGH, i, j, k,
              1, dyi, alp,
              betax, betay, betaz,
              gxx, gxy, gxz,
              gyy, gyz, gzz,
              &ADM_metric_derivs_y);

        ghl_metric_quantities ADM_metric_derivs_z;
        GRHayLMHD_compute_metric_derivs(
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
