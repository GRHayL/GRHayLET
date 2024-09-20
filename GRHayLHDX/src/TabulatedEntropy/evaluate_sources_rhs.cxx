#include "GRHayLHDX.h"

void GRHayLHDX_tabulated_entropy_evaluate_sources_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_tabulated_entropy_evaluate_sources_rhs;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dxi = 1.0/CCTK_DELTA_SPACE(0);
  const CCTK_REAL dyi = 1.0/CCTK_DELTA_SPACE(1);
  const CCTK_REAL dzi = 1.0/CCTK_DELTA_SPACE(2);

  const Loop::GF3D2layout layout(cctkGH, {1, 1, 1});

  grid.loop_int<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index index(layout, p.I);

    // These variables have no source terms
    rho_star_rhs(index) = 0.0;
    ent_star_rhs(index) = 0.0;
    Ye_star_rhs(index)  = 0.0;

    ghl_metric_quantities ADM_metric;
    ghl_initialize_metric(
          ccc_lapse(index),
          ccc_betax(index), ccc_betay(index), ccc_betaz(index),
          ccc_gxx(index), ccc_gxy(index), ccc_gxz(index),
          ccc_gyy(index), ccc_gyz(index), ccc_gzz(index),
          &ADM_metric);

    ghl_extrinsic_curvature curv;
    ghl_initialize_extrinsic_curvature(
          ccc_kxx(index), ccc_kxy(index), ccc_kxz(index),
          ccc_kyy(index), ccc_kyz(index), ccc_kzz(index),
          &curv);

    ghl_primitive_quantities prims;
    prims.BU[0] = prims.BU[1] = prims.BU[2] = 0.0;
    prims.rho         = rho(index);
    prims.press       = press(index);
    prims.vU[0]       = vx(index);
    prims.vU[1]       = vy(index);
    prims.vU[2]       = vz(index);
    prims.entropy     = entropy(index);
    prims.Y_e         = Ye(index);
    prims.temperature = temperature(index);

    bool speed_limited;
    ghl_error_codes_t error = ghl_limit_v_and_compute_u0(ghl_params, &ADM_metric, &prims, &speed_limited);
    if(error)
      ghl_read_error_codes(error);

    const Loop::GF3D2index indm2x(layout, p.I - 2*p.DI[0]);
    const Loop::GF3D2index indm1x(layout, p.I -   p.DI[0]);
    const Loop::GF3D2index indp1x(layout, p.I +   p.DI[0]);
    const Loop::GF3D2index indp2x(layout, p.I + 2*p.DI[0]);

    ghl_metric_quantities ADM_metric_derivs_x;
    GRHayLHDX_compute_metric_derivs(
          dxi, indm2x, indm1x, indp1x, indp2x,
          ccc_lapse, ccc_betax, ccc_betay, ccc_betaz,
          ccc_gxx, ccc_gxy, ccc_gxz,
          ccc_gyy, ccc_gyz, ccc_gzz,
          &ADM_metric_derivs_x);

    const Loop::GF3D2index indm2y(layout, p.I - 2*p.DI[1]);
    const Loop::GF3D2index indm1y(layout, p.I -   p.DI[1]);
    const Loop::GF3D2index indp1y(layout, p.I +   p.DI[1]);
    const Loop::GF3D2index indp2y(layout, p.I + 2*p.DI[1]);

    ghl_metric_quantities ADM_metric_derivs_y;
    GRHayLHDX_compute_metric_derivs(
          dyi, indm2y, indm1y, indp1y, indp2y,
          ccc_lapse, ccc_betax, ccc_betay, ccc_betaz,
          ccc_gxx, ccc_gxy, ccc_gxz,
          ccc_gyy, ccc_gyz, ccc_gzz,
          &ADM_metric_derivs_y);

    const Loop::GF3D2index indm2z(layout, p.I - 2*p.DI[2]);
    const Loop::GF3D2index indm1z(layout, p.I -   p.DI[2]);
    const Loop::GF3D2index indp1z(layout, p.I +   p.DI[2]);
    const Loop::GF3D2index indp2z(layout, p.I + 2*p.DI[2]);

    ghl_metric_quantities ADM_metric_derivs_z;
    GRHayLHDX_compute_metric_derivs(
          dzi, indm2z, indm1z, indp1z, indp2z,
          ccc_lapse, ccc_betax, ccc_betay, ccc_betaz,
          ccc_gxx, ccc_gxy, ccc_gxz,
          ccc_gyy, ccc_gyz, ccc_gzz,
          &ADM_metric_derivs_z);

    ghl_conservative_quantities cons_source;
    ghl_calculate_source_terms(
          ghl_eos, &prims, &ADM_metric,
          &ADM_metric_derivs_x,
          &ADM_metric_derivs_y,
          &ADM_metric_derivs_z,
          &curv, &cons_source);

    tau_rhs(index) = cons_source.tau;
    Stildex_rhs(index) = cons_source.SD[0];
    Stildey_rhs(index) = cons_source.SD[1];
    Stildez_rhs(index) = cons_source.SD[2];
  }); // ccc loop interior
}
