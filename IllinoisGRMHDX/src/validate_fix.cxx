#include "IllinoisGRMHDX.hxx"

extern "C" void GRHayLHDX_validate_entstar(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_validate_entstar;
  DECLARE_CCTK_PARAMETERS;

  constexpr std::array<int, Loop::dim> ccctype = {1, 1, 1};
  const Loop::GF3D2layout ccc_layout(cctkGH, ccctype);

  grid.loop_all<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index index(ccc_layout, p.I);

    ent_star(index) = 0;
    ent_star_rhs(index) = 0;
  }); // ccc loop everywhere
}

extern "C" void GRHayLHDX_validate_yestar(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_validate_yestar;
  DECLARE_CCTK_PARAMETERS;

  constexpr std::array<int, Loop::dim> ccctype = {1, 1, 1};
  const Loop::GF3D2layout ccc_layout(cctkGH, ccctype);

  grid.loop_all<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index index(ccc_layout, p.I);

    Ye_star(index) = 0;
    Ye_star_rhs(index) = 0;
  }); // ccc loop everywhere
}
