/********************************
 * CONVERT ET ID TO IllinoisGRMHD
 *
 * Written in 2014 by Zachariah B. Etienne
 *
 * Sets metric & MHD variables needed
 * by IllinoisGRMHD, converting from
 * HydroBase and ADMBase.
 ********************************/

#include "GRHayLHDX.h"

extern "C" void GRHayLHDX_prims_to_conservs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_prims_to_conservs;
  DECLARE_CCTK_PARAMETERS;

  constexpr std::array<int, Loop::dim> indextype = {1, 1, 1};
  const Loop::GF3D2layout layout(cctkGH, indextype);

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
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
    ghl_initialize_primitives(
          rho(index), press(index), eps(index),
          vx(index), vy(index), vz(index),
          0.0, 0.0, 0.0,
          entropy(index), Ye(index), temperature(index),
          &prims);

    ghl_conservative_quantities cons;
    //This applies inequality fixes on the conservatives
    const int speed_limited CCTK_ATTRIBUTE_UNUSED = ghl_enforce_primitive_limits_and_compute_u0(
          ghl_params, ghl_eos, &ADM_metric, &prims);
    //This computes the conservatives from the new primitives
    ghl_compute_conservs(
          &ADM_metric, &metric_aux, &prims, &cons);

    CCTK_REAL dummy1, dummy2, dummy3;
    CCTK_REAL dummy4, dummy5, dummy6;
    ghl_return_primitives(
          &prims,
          &rho(index), &press(index), &eps(index),
          &vx(index), &vy(index), &vz(index),
          &dummy1, &dummy2, &dummy3,
          &entropy(index), &Ye(index), &temperature(index));

    ghl_return_conservatives(
          &cons,
          &rho_star(index), &tau(index),
          &Stildex(index), &Stildey(index), &Stildez(index),
          &ent_star(index), &Ye_star(index));
  }); // ccc loop everywhere
}
