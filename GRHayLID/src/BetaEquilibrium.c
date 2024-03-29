#include "GRHayLID.h"

void GRHayLID_BetaEquilibrium( CCTK_ARGUMENTS ) {

  CCTK_INFO("Starting routine to impose neutrino free beta-equilibrium");

  DECLARE_CCTK_ARGUMENTS_GRHayLID_BetaEquilibrium;
  DECLARE_CCTK_PARAMETERS;

  if(ghl_eos->eos_type != ghl_eos_tabulated)
    CCTK_ERROR("GRHayL can only impose beta-equilibrium if tabulated EOS is used");

  if(beq_temperature < ghl_eos->table_T_min || beq_temperature > ghl_eos->table_T_max)
    CCTK_VERROR("Parameter beq_temperature (%g) exceeds table bounds [%g, %g]",
                beq_temperature, ghl_eos->table_T_min, ghl_eos->table_T_max);

  CHECK_PARAMETER(beq_temperature);

  ghl_tabulated_compute_Ye_of_rho_beq_constant_T(beq_temperature, ghl_eos);

  for(int k=0; k<cctk_lsh[2]; k++) {
    for(int j=0; j<cctk_lsh[1]; j++) {
      for(int i=0; i<cctk_lsh[0]; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);

        const double rhoL  = rho[index];

        if( rhoL <= 1.01*ghl_eos->rho_atm ) {
          press      [index] = ghl_eos->press_atm;
          eps        [index] = ghl_eos->eps_atm;
          Y_e        [index] = ghl_eos->Y_e_atm;
          temperature[index] = ghl_eos->T_atm;
        } else {
          const double tempL = beq_temperature;
          const double YeL   = ghl_tabulated_compute_Ye_from_rho(ghl_eos, rhoL);

          double pressL, epsL;
          ghl_tabulated_compute_P_eps_from_T(ghl_eos, rhoL, YeL, tempL, &pressL, &epsL);

          press      [index] = pressL;
          eps        [index] = epsL;
          Y_e        [index] = YeL;
          temperature[index] = tempL;
        }
      }
    }
  }
  CCTK_INFO("Finished imposing neutrino free beta-equilibrium");
}
