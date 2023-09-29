#include "GRHayLIDX.h"

void GRHayLIDX_compute_entropy_hybrid(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLIDX_compute_entropy_hybrid;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D2layout layout(cctkGH, {1, 1, 1});

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index index(layout, p.I);
    entropy(index) = ghl_hybrid_compute_entropy_function(ghl_eos, rho(index), press(index));
  });
}

void GRHayLIDX_compute_entropy_tabulated(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLIDX_compute_entropy_tabulated;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D2layout layout(cctkGH, {1, 1, 1});

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index index(layout, p.I);

    double P_local, eps_local, S_local;
    ghl_tabulated_enforce_bounds_rho_Ye_T(ghl_eos, &rho(index), &Y_e(index), &temperature(index));
    ghl_tabulated_compute_P_eps_S_from_T(
          ghl_eos, rho(index), Ye(index), temperature(index),
          &P_local, &eps_local, &S_local);

    entropy(index) = S_local;
  });
}
