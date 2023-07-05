#include "GRHayLID.h"

/*
 *
 * (c) 2021 Leo Werneck
 *
 * This is the thorn's driver function, responsible
 * for setting the initial data to that of an isotropic
 * gas of constant density, temperature, and electron
 * fraction in Minkowski space.
 */
void GRHayLID_IsotropicGas(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLID_IsotropicGas;
  DECLARE_CCTK_PARAMETERS;

  if(!CCTK_EQUALS(EOS_type, "tabulated"))
    CCTK_VERROR("IsotropicGas initial data is only defined for tabulated EOS. Please change GRHayLib::EOS_type to \"tabulated\" in the parfile.");
  if(!CCTK_EQUALS(initial_Y_e, "GRHayLID"))
    CCTK_VERROR("To use IsotropicGas initial data, please add HydroBase::initial_Y_e=\"GRHayLID\" to the parfile.");
  if(!CCTK_EQUALS(initial_temperature, "GRHayLID"))
    CCTK_VERROR("To use IsotropicGas initial data, please add HydroBase::initial_temperature=\"GRHayLID\" to the parfile.");

  CHECK_PARAMETER(IsotropicGas_rho);
  CHECK_PARAMETER(IsotropicGas_Y_e);
  CHECK_PARAMETER(IsotropicGas_temperature);

  CCTK_INFO("Beginning IsotropicGas initial data");

  CCTK_REAL IsotropicGas_press, IsotropicGas_eps, IsotropicGas_entropy;
  ghl_tabulated_compute_P_eps_S_from_T(
        ghl_eos,
        IsotropicGas_rho,
        IsotropicGas_Y_e,
        IsotropicGas_temperature,
        &IsotropicGas_press,
        &IsotropicGas_eps,
        &IsotropicGas_entropy);

#pragma omp parallel for
  for(int k=0; k<cctk_lsh[2]; k++) {
    for(int j=0; j<cctk_lsh[1]; j++) {
      for(int i=0; i<cctk_lsh[0]; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        const int ind4x = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0);
        const int ind4y = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1);
        const int ind4z = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2);

        Y_e[index] = IsotropicGas_Y_e;
        eps[index] = IsotropicGas_eps;
        press[index] = IsotropicGas_press;
        rho[index] = IsotropicGas_rho;
        temperature[index] = IsotropicGas_temperature;

        vel[ind4x] = 0;
        vel[ind4y] = 0;
        vel[ind4z] = 0;

        if(CCTK_EQUALS(initial_entropy, "GRHayLID")) entropy[index] = IsotropicGas_entropy;
      }
    }
  }

  CCTK_INFO("All done!");
}
