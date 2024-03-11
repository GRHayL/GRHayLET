#include "GRHayLHDX.h"

extern "C" void GRHayLHDX_tabulated_entropy_perturb_primitives(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_tabulated_entropy_perturb_primitives;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D2layout layout(cctkGH, {1, 1, 1});

  srand(random_seed); // Use srand() as rand() is thread-safe.
  grid.loop_all<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index index(layout, p.I);
    rho(index)         *= one_plus_pert(random_pert);
    press(index)       *= one_plus_pert(random_pert);
    vx(index)          *= one_plus_pert(random_pert);
    vy(index)          *= one_plus_pert(random_pert);
    vz(index)          *= one_plus_pert(random_pert);
    entropy(index)     *= one_plus_pert(random_pert);
    Ye(index)          *= one_plus_pert(random_pert);
    temperature(index) *= one_plus_pert(random_pert);
  }); // ccc loop everywhere
}
