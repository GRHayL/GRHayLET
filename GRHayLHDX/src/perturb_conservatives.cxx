#include "GRHayLHDX.h"

extern "C" void GRHayLHDX_perturb_conservatives(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_perturb_conservatives;
  DECLARE_CCTK_PARAMETERS;

  srand(random_seed); // Use srand() as rand() is thread-safe.
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index index(layout, p.I);
    rho_star(index) *= one_plus_pert(random_pert);
    tau(index)      *= one_plus_pert(random_pert);
    Stildex(index)  *= one_plus_pert(random_pert);
    Stildey(index)  *= one_plus_pert(random_pert);
    Stildez(index)  *= one_plus_pert(random_pert);
  });
}
