#include "GRHayLIDX.h"

extern "C" void GRHayLIDX_ConstantDensitySphere(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLIDX_ConstantDensitySphere;
  DECLARE_CCTK_PARAMETERS;

  if(!CCTK_EQUALS(EOS_type, "tabulated"))
    CCTK_VERROR("ConstantDensitySphere initial data is only defined for tabulated EOS. Please change GRHayLib::EOS_type to \"tabulated\" in the parfile.");

  CHECK_PARAMETER(ConstantDensitySphere_sphere_radius);
  CHECK_PARAMETER(ConstantDensitySphere_rho_interior);
  CHECK_PARAMETER(ConstantDensitySphere_Y_e_interior);
  CHECK_PARAMETER(ConstantDensitySphere_T_interior);
  CHECK_PARAMETER(ConstantDensitySphere_rho_exterior);
  CHECK_PARAMETER(ConstantDensitySphere_Y_e_exterior);
  CHECK_PARAMETER(ConstantDensitySphere_T_exterior);

  CCTK_INFO("Beginning ConstantDensitySphere initial data");

  // Compute hydro quantities inside and outside the sphere
  CCTK_REAL P_interior, eps_interior, S_interior;
  ghl_tabulated_compute_P_eps_S_from_T(
        ghl_eos,
        ConstantDensitySphere_rho_interior,
        ConstantDensitySphere_Y_e_interior,
        ConstantDensitySphere_T_interior,
        &P_interior,
        &eps_interior,
        &S_interior );

  CCTK_REAL P_exterior, eps_exterior, S_exterior;
  ghl_tabulated_compute_P_eps_S_from_T(
        ghl_eos,
        ConstantDensitySphere_rho_exterior,
        ConstantDensitySphere_Y_e_exterior,
        ConstantDensitySphere_T_exterior,
        &P_exterior,
        &eps_exterior,
        &S_exterior );

  const Loop::GF3D2layout layout(cctkGH, {1, 1, 1});

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index index(layout, p.I);

    const double r = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
    velx(index) = 0;
    vely(index) = 0;
    velz(index) = 0;
    if( r > ConstantDensitySphere_sphere_radius ) {
      // Outside the sphere
      rho(index)         = ConstantDensitySphere_rho_exterior;
      Ye(index)          = ConstantDensitySphere_Y_e_exterior;
      temperature(index) = ConstantDensitySphere_T_exterior;
      press(index)       = P_exterior;
      eps(index)         = eps_exterior;
      entropy(index)     = S_exterior;
    }
    else {
      // Inside the sphere
      rho(index)         = ConstantDensitySphere_rho_interior;
      Ye(index)          = ConstantDensitySphere_Y_e_interior;
      temperature(index) = ConstantDensitySphere_T_interior;
      press(index)       = P_interior;
      eps(index)         = eps_interior;
      entropy(index)     = S_interior;
    }
  });

  CCTK_INFO("All done!");
}
