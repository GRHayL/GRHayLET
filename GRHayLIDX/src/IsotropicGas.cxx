#include "GRHayLIDX.h"

void GRHayLIDX_IsotropicGas(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLIDX_IsotropicGas;
  DECLARE_CCTK_PARAMETERS;

  if(!CCTK_EQUALS(EOS_type, "Tabulated"))
    CCTK_ERROR("IsotropicGas initial data is only defined for tabulated EOS. Please change GRHayLib::EOS_type to \"Tabulated\" in the parfile.");

  CHECK_PARAMETER(IsotropicGas_rho);
  CHECK_PARAMETER(IsotropicGas_Y_e);
  CHECK_PARAMETER(IsotropicGas_temperature);

  CCTK_INFO("Beginning IsotropicGas initial data");

  CCTK_REAL IsotropicGas_press, IsotropicGas_eps;
  ghl_tabulated_compute_P_eps_from_T(
        ghl_eos,
        IsotropicGas_rho,
        IsotropicGas_Y_e,
        IsotropicGas_temperature,
        &IsotropicGas_press,
        &IsotropicGas_eps);

  const Loop::GF3D2layout layout(cctkGH, {1, 1, 1});

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index index(layout, p.I);

    velx(index) = 0;
    vely(index) = 0;
    velz(index) = 0;

    rho(index)         = IsotropicGas_rho;
    Ye(index)          = IsotropicGas_Y_e;
    temperature(index) = IsotropicGas_temperature;
    press(index)       = IsotropicGas_press;
    eps(index)         = IsotropicGas_eps;
  });

  CCTK_INFO("All done!");
}
