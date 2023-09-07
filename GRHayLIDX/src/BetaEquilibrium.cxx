#include "GRHayLIDX.h"

void GRHayLIDX_BetaEquilibrium( CCTK_ARGUMENTS ) {

  CCTK_INFO("Starting routine to impose neutrino free beta-equilibrium");

  DECLARE_CCTK_ARGUMENTSX_GRHayLIDX_BetaEquilibrium;
  DECLARE_CCTK_PARAMETERS;

  if( ghl_eos->eos_type != ghl_eos_tabulated )
    CCTK_ERROR("GRHayL can only impose beta-equilibrium if tabulated EOS is used");

  if( beq_temperature < ghl_eos->table_T_min || beq_temperature > ghl_eos->table_T_max )
    CCTK_VERROR("Parameter beq_temperature (%g) exceeds table bounds [%g, %g]",
                beq_temperature, ghl_eos->table_T_min, ghl_eos->table_T_max);

  CHECK_PARAMETER(beq_temperature);

  double *Ye_of_lr;
  ghl_tabulated_compute_Ye_of_rho_beq_constant_T(ghl_eos, beq_temperature, &Ye_of_lr);

  const Loop::GF3D2layout layout(cctkGH, {1, 1, 1});

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index index(layout, p.I);

    const double rhoL  = rho(index);

    if( rhoL <= 1.01*ghl_eos->rho_atm ) {
      press      (index) = ghl_eos->press_atm;
      eps        (index) = ghl_eos->eps_atm;
      Ye         (index) = ghl_eos->Y_e_atm;
      temperature(index) = ghl_eos->T_atm;
    }
    else {
      const double tempL = beq_temperature;
      const double YeL   = ghl_tabulated_get_Ye_from_rho(
                              ghl_eos->N_rho, ghl_eos->table_logrho, Ye_of_lr, rhoL);

      double pressL, epsL;
      ghl_tabulated_compute_P_eps_from_T(ghl_eos, rhoL, YeL, tempL, &pressL, &epsL);

      press      (index) = pressL;
      eps        (index) = epsL;
      Ye         (index) = YeL;
      temperature(index) = tempL;
    }
  });

  free(Ye_of_lr);
  CCTK_INFO("Finished imposing neutrino free beta-equilibrium");
}
