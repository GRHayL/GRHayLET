#include "IGM.h"

void IllinoisGRMHD_support_deprecated(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  char valbuf[100];
  bool uses_deprecated = false;
  int ierr;

  if( fabs(Psi6threshold - 1e100) > 1e-10 ) {
    uses_deprecated = true;
    CCTK_WARN(CCTK_WARN_COMPLAIN, "WARNING: parameter IllinoisGRMHD::Psi6threshold is deprecated.\n"
              "Please switch to using GRHayLib::Psi6threshold.\n"
              "Setting GRHayLib::Psi6threshold to match IllinoisGRMHD::Psi6threshold.");
    sprintf (valbuf,"%.17f", Psi6threshold);
    ierr = CCTK_ParameterSet("Psi6threshold", "GRHayLib", valbuf);
  }

  if( fabs(rho_b_atm - 1e100) > 1e-10 ) {
    uses_deprecated = true;
    CCTK_WARN(CCTK_WARN_COMPLAIN, "WARNING: parameter IllinoisGRMHD::rho_b_atm is deprecated.\n"
              "Please switch to using GRHayLib::rho_b_atm.\n"
              "and GRHayLib::rho_b_min.\n"
              "Setting GRHayLib::rho_b_atm and GRHayLib::rho_b_min to match IllinoisGRMHD::rho_b_atm.");
    sprintf (valbuf,"%.17f", rho_b_atm);
    ierr = CCTK_ParameterSet("rho_b_atm", "GRHayLib", valbuf);
    ierr = CCTK_ParameterSet("rho_b_min", "GRHayLib", valbuf);
  }

  if( rho_b_max > 0 ) {
    uses_deprecated = true;
    CCTK_WARN(CCTK_WARN_COMPLAIN, "WARNING: parameter IllinoisGRMHD::rho_b_max is deprecated.\n"
              "Please switch to using GRHayLib::rho_b_max.\n"
              "Setting GRHayLib::rho_b_max to match IllinoisGRMHD::rho_b_max.");
    sprintf (valbuf,"%.17f", rho_b_max);
    ierr = CCTK_ParameterSet("rho_b_max", "GRHayLib", valbuf);
  }

  if( neos != 1 ) {
    uses_deprecated = true;
    CCTK_WARN(CCTK_WARN_COMPLAIN, "WARNING: parameter IllinoisGRMHD::neos is deprecated.\n"
              "Please switch to using GRHayLib::neos.\n"
              "Setting GRHayLib::neos to match IllinoisGRMHD::neos.");
    sprintf (valbuf,"%d", neos);
    ierr = CCTK_ParameterSet("neos", "GRHayLib", valbuf);
  }

  if( gamma_th > 0 ) {
    uses_deprecated = true;
    CCTK_WARN(CCTK_WARN_COMPLAIN, "WARNING: parameter IllinoisGRMHD::gamma_th is deprecated.\n"
              "Please switch to using GRHayLib::Gamma_th.\n"
              "Setting GRHayLib::Gamma_th to match IllinoisGRMHD::gamma_th.");
    sprintf (valbuf,"%.17f", gamma_th);
    ierr = CCTK_ParameterSet("Gamma_th", "GRHayLib", valbuf);
  }
    
  if( fabs(GAMMA_SPEED_LIMIT - 10.0) < 1e-10 ) {
    uses_deprecated = true;
    CCTK_WARN(CCTK_WARN_COMPLAIN, "WARNING: parameter IllinoisGRMHD::GAMMA_SPEED_LIMIT is deprecated.\n"
              "Please switch to using GRHayLib::max_lorenz_factor.\n"
              "Setting GRHayLib::max_lorenz_factor to match IllinoisGRMHD::GAMMA_SPEED_LIMIT.");
    sprintf (valbuf,"%.17f", GAMMA_SPEED_LIMIT);
    ierr = CCTK_ParameterSet("max_lorenz_factor", "GRHayLib", valbuf);
  }

//Still needs to check Gamma, but IDK how that's set in IGM originally

  if( K_poly != 1 ) {
    uses_deprecated = true;
    CCTK_WARN(CCTK_WARN_COMPLAIN, "WARNING: parameter IllinoisGRMHD::K_poly is deprecated.\n"
              "Please switch to using GRHayLib::k_ppoly0.\n"
              "Setting GRHayLib::k_ppoly0 to match IllinoisGRMHD::K_poly.");
    sprintf (valbuf,"%.17f", K_poly);
    ierr = CCTK_ParameterSet("k_ppoly0", "GRHayLib", valbuf);
  }
    
  if( damp_lorenz != 0.0 ) {
    uses_deprecated = true;
    CCTK_WARN(CCTK_WARN_COMPLAIN, "WARNING: parameter IllinoisGRMHD::damp_lorenz is deprecated.\n"
              "Please switch to using GRHayLib::Lorenz_damping_factor.\n"
              "Setting GRHayLib::Lorenz_damping_factor to match IllinoisGRMHD::damp_lorenz.");
    sprintf (valbuf,"%.17f", damp_lorenz);
    ierr = CCTK_ParameterSet("Lorenz_damping_factor", "GRHayLib", valbuf);
  }
    
  if(conserv_to_prims_debug == 1 )
    CCTK_WARN(CCTK_WARN_COMPLAIN, "WARNING: output provided by IllinoisGRMHD::conserv_to_prims_debug no longer supported.");

  if(CCTK_IsThornActive("ID_converter_ILGRMHD")) {
    uses_deprecated = true;
    CCTK_VWARN(CCTK_WARN_COMPLAIN, "###############################################\n"
               "#                   WARNING!                  #\n"
               "###############################################\n"
               "#                   WARNING!                  #\n"
               "###############################################\n"
               "#                   WARNING!                  #\n"
               "###############################################\n"
               "#    The thorn ID_converter_ILGRMHD is now    #\n"
               "# deprecated, and it will be removed entirely #\n"
               "#          in the ET_2024_05 release.         #\n"
               "###############################################\n");

    const CCTK_REAL *rand_ptr = CCTK_ParameterGet("random_pert", "ID_converter_ILGRMHD", NULL);
    if ( (*rand_ptr) != 0 ) {
      CCTK_VWARN(CCTK_WARN_COMPLAIN, "WARNING: parameter ID_converter_ILGRMHD::random_pert is deprecated.\n"
                 "Please switch to using IllinoisGRMHD::perturb_initial_data.\n"
                 "Setting IllinoisGRMHD::perturb_initial_data to match ID_converter_ILGRMHD::random_pert.");
      sprintf (valbuf,"%.17f", (*rand_ptr));
      ierr = CCTK_ParameterSet("perturb_initial_data", "IllinoisGRMHD", valbuf);
    }

    const CCTK_INT *seed_ptr = CCTK_ParameterGet("random_seed", "ID_converter_ILGRMHD", NULL);
    if ( (*seed_ptr) != 0 ) {
      CCTK_VWARN(CCTK_WARN_COMPLAIN, "WARNING: parameter ID_converter_ILGRMHD::random_seed is deprecated.\n"
                 "Please switch to using IllinoisGRMHD::random_seed.\n"
                 "Setting IllinoisGRMHD::random_seed to match ID_converter_ILGRMHD::random_seed.");
      sprintf (valbuf,"%d", (*seed_ptr));
      ierr = CCTK_ParameterSet("random_seed", "IllinoisGRMHD", valbuf);
    }

    const bool *phr_ptr = CCTK_ParameterGet("pure_hydro_run", "ID_converter_ILGRMHD", NULL);
    if (*phr_ptr == true)
      CCTK_VWARN(CCTK_WARN_COMPLAIN, "WARNING: parameter ID_converter_ILGRMHD::pure_hydro_run is deprecated.\n"
                 "Note that previous code would not have set A=0 properly, so the behavior\n"
                 "is identical without this parameter. For true hydro-only code,\n"
                 "switch to using the GRHayLHD thorn.");

    const CCTK_REAL *Gamma_old_ptr = CCTK_ParameterGet("Gamma_Initial", "ID_converter_ILGRMHD", NULL);
    const CCTK_REAL *Gamma_new_ptr = CCTK_ParameterGet("Gamma_ppoly_in[0]", "GRHayLib", NULL);
    if( fabs((*Gamma_old_ptr) - (*Gamma_new_ptr)) > 1e-10 )
      CCTK_VERROR("IllinoisGRMHD does not support different Gamma values! This parameter file has\n"
                  "    ID_converter_ILGRMHD::Gamma_Initial != GRHayLib::Gamma_ppoly_in[0].\n"
                  "Remove Gamma_Initial from the parfile and just use GRHayLib::Gamma_ppoly_in[0].");

    ierr = CCTK_ParameterSet("k_ppoly0", "GRHayLib", valbuf);
    const CCTK_REAL *K_old_ptr = CCTK_ParameterGet("K_Initial", "ID_converter_ILGRMHD", NULL);
    const CCTK_REAL *K_new_ptr = CCTK_ParameterGet("k_ppoly0", "GRHayLib", NULL);
    if( fabs((*K_old_ptr) - (*K_new_ptr)) > 1e-10 )
      CCTK_VERROR("IllinoisGRMHD does not support different K values! This parameter file has\n"
                  "    ID_converter_ILGRMHD::K_Initial != GRHayLib::k_ppoly0.\n"
                  "Remove K_Initial from the parfile and just use GRHayLib::k_ppoly0.");
  }

  if(CCTK_IsThornActive("Convert_to_HydroBase")) {
    uses_deprecated = true;
    CCTK_VWARN(CCTK_WARN_COMPLAIN, "###############################################\n"
               "#                   WARNING!                  #\n"
               "###############################################\n"
               "#                   WARNING!                  #\n"
               "###############################################\n"
               "#                   WARNING!                  #\n"
               "###############################################\n"
               "#    The thorn Convert_to_HydroBase is now    #\n"
               "# deprecated, and it will be removed entirely #\n"
               "#          in the ET_2024_05 release.         #\n"
               "###############################################\n");

    const CCTK_INT *cth_ptr = CCTK_ParameterGet("Convert_to_HydroBase_every", "Convert_to_HydroBase", NULL);
    if (*cth_ptr != 0) {
      CCTK_VWARN(CCTK_WARN_COMPLAIN, "WARNING: parameter Convert_to_HydroBase::Convert_to_HydroBase_every is deprecated.\n"
                 "Please switch to using IllinoisGRMHD::Convert_to_HydroBase_every.\n"
                 "Setting IllinoisGRMHD::Convert_to_HydroBase_every to match Convert_to_HydroBase::Convert_to_HydroBase_every.");
      sprintf (valbuf,"%d", (*cth_ptr));
      ierr = CCTK_ParameterSet("Convert_to_HydroBase_every", "IllinoisGRMHD", valbuf);
    }
  }

  if(uses_deprecated)
    CCTK_INFO(
"Notice: You are currently using a parfile with deprecated features.\n"
"A large update to IllinoisGRMHD refactored the core code of IllinoisGRMHD\n"
"into a modular library (GRHayL), which is provided by the GRHayLib thorn.\n"
"It also simplified/condensed many repeated pieces of code and parameters,\n"
"eliminating the ID_converter_ILGRMHD and Convert_to_HydroBase thorns.\n"
"Deprecated parameters have warnings for each parameter which explain\n"
"what the equivalent parameter is in the new IllinoisGRMHD. Deprecated\n"
"features will be **removed** in the next release, so please begin\n"
"updating your parfiles in preparation for this change. The full\n"
"parfile changelog is provided in repos/GRHayLET/IllinoisGRMHD/README.md.\n");
}

void IllinoisGRMHD_support_deprecated_tau_atm(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  if( fabs(tau_atm - 1e100) > 1e-10 ) {
    CCTK_WARN(CCTK_WARN_COMPLAIN, "WARNING: parameter IllinoisGRMHD::tau_atm is deprecated.\n"
              "The value of tau_atm is automatically computed by GRHayLib.\n"
              "Setting ghl_eos->tau_atm to match IllinoisGRMHD::tau_atm.");
    ghl_eos->tau_atm = tau_atm;
  }
}
