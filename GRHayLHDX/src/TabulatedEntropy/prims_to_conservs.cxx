#include "GRHayLHDX.h"

extern "C" void GRHayLHDX_tabulated_entropy_prims_to_conservs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_tabulated_entropy_prims_to_conservs;
  DECLARE_CCTK_PARAMETERS;

  constexpr std::array<int, Loop::dim> indextype = {1, 1, 1};
  const Loop::GF3D2layout layout(cctkGH, indextype);

  grid.loop_all<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index index(layout, p.I);

    ghl_metric_quantities ADM_metric;
    ghl_initialize_metric(
          ccc_lapse(index),
          ccc_betax(index), ccc_betay(index), ccc_betaz(index),
          ccc_gxx(index), ccc_gxy(index), ccc_gxz(index),
          ccc_gyy(index), ccc_gyz(index), ccc_gzz(index),
          &ADM_metric);

    ghl_ADM_aux_quantities metric_aux;
    ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

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

    ghl_conservative_quantities cons;
    const int speed_limited CCTK_ATTRIBUTE_UNUSED = ghl_enforce_primitive_limits_and_compute_u0(
          ghl_params, ghl_eos, &ADM_metric, &prims);
    ghl_compute_conservs(
          &ADM_metric, &metric_aux, &prims, &cons);

    rho(index)         = prims.rho;
    press(index)       = prims.press;
    eps(index)         = prims.eps;
    u0(index)          = prims.u0;
    vx(index)          = prims.vU[0];
    vy(index)          = prims.vU[1];
    vz(index)          = prims.vU[2];
    entropy(index)     = prims.entropy;
    Ye(index)          = prims.Y_e;
    temperature(index) = prims.temperature;

    rho_star(index) = cons.rho;
    tau(index)      = cons.tau;
    Stildex(index)  = cons.SD[0];
    Stildey(index)  = cons.SD[1];
    Stildez(index)  = cons.SD[2];
    ent_star(index) = cons.entropy;
    Ye_star(index)  = cons.Y_e;
  }); // ccc loop everywhere
}
