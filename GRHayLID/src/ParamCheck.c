#include "GRHayLID.h"

void GRHayLID_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLID_ParamCheck;
  DECLARE_CCTK_PARAMETERS;

  const int one_d = CCTK_EQUALS(initial_hydro, "HydroTest1D");
  const int tabulated_id = CCTK_EQUALS(initial_hydro, "IsotropicGas") ||
                           CCTK_EQUALS(initial_hydro, "ConstantDensitySphere");
  const int avec = CCTK_EQUALS(initial_Avec, "GRHayLID");
  const int bvec = CCTK_EQUALS(initial_Bvec, "GRHayLID");
  if(avec != bvec || ((avec || bvec) && (!one_d || !initialize_magnetic_quantities)))
    CCTK_ERROR("GRHayLID magnetics require HydroTest1D, initialize_magnetic_quantities=yes, and both initial_Avec and initial_Bvec=GRHayLID");
  if(one_d && initialize_magnetic_quantities && !(avec && bvec))
    CCTK_ERROR("Disable initialize_magnetic_quantities or select both GRHayLID magnetic destinations");
  if(one_d && !(CCTK_EQUALS(EOS_type, "Simple") || CCTK_EQUALS(EOS_type, "Hybrid")))
    CCTK_ERROR("HydroTest1D requires EOS_type=Simple or Hybrid");
  if(tabulated_id && (!CCTK_EQUALS(EOS_type, "Tabulated") ||
     !CCTK_EQUALS(initial_Y_e, "GRHayLID") || !CCTK_EQUALS(initial_temperature, "GRHayLID")))
    CCTK_ERROR("Gas/sphere require Tabulated EOS and initial_Y_e=initial_temperature=GRHayLID");
  if(!tabulated_id && !impose_beta_equilibrium &&
     (CCTK_EQUALS(initial_Y_e, "GRHayLID") || CCTK_EQUALS(initial_temperature, "GRHayLID")))
    CCTK_ERROR("GRHayLID Ye/temperature selectors require gas, sphere, or impose_beta_equilibrium");
  if(impose_beta_equilibrium && !CCTK_EQUALS(EOS_type, "Tabulated"))
    CCTK_ERROR("Beta equilibrium requires EOS_type=Tabulated");
  if((impose_beta_equilibrium || CCTK_EQUALS(initial_entropy, "GRHayLID")) &&
     CCTK_EQUALS(initial_hydro, "none"))
    CCTK_ERROR("Standalone beta/entropy requires an initialized hydro producer in HydroBase_Initial");
  if(CCTK_EQUALS(initial_entropy, "GRHayLID")) {
    if(CCTK_EQUALS(EOS_type, "Tabulated")) {
      if(!impose_beta_equilibrium && (CCTK_EQUALS(initial_Y_e, "none") ||
                                    CCTK_EQUALS(initial_temperature, "none")))
        CCTK_ERROR("Tabulated entropy requires initialized Ye and temperature producers, or beta equilibrium");
    } else if(CCTK_EQUALS(EOS_type, "Hybrid") || CCTK_EQUALS(EOS_type, "Simple")) {
      if(!allow_native_entropy_proxy)
        CCTK_ERROR("Simple/Hybrid entropy is a GRHayL native proxy, not physical specific entropy; set allow_native_entropy_proxy=yes only for compatible consumers");
    } else {
      CCTK_ERROR("GRHayLID entropy requires EOS_type=Simple, Hybrid, or Tabulated");
    }
  }
  if(one_d && CCTK_EQUALS(initial_data_1D, "sound wave") &&
     (!isfinite(wave_amplitude) || wave_amplitude < 0.0 || wave_amplitude >= 1.0))
    CCTK_ERROR("Sound-wave velocity amplitude must be finite and in [0,1)");
}
