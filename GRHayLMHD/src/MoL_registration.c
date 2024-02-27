//--------------------------------------------------------------------------
// Register with the time stepper
// (MoL thorn, found in arrangements/CactusBase/MoL)
// To understand this, read documentation in arrangements/CactusBase/MoL/doc
//--------------------------------------------------------------------------

#include "GRHayLMHD.h"
#include "Symmetry.h"

void GRHayLMHD_RegisterVars(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0, group, rhs;

  //***********************************************
  // Register evolution & RHS gridfunction variables

  /* Ax and Ax_rhs */
  group = CCTK_VarIndex("GRHayLMHD::Ax");
  rhs = CCTK_VarIndex("GRHayLMHD::Ax_rhs");
  ierr += MoLRegisterEvolved(group, rhs);

  /* Ay and Ay_rhs */
  group = CCTK_VarIndex("GRHayLMHD::Ay");
  rhs = CCTK_VarIndex("GRHayLMHD::Ay_rhs");
  ierr += MoLRegisterEvolved(group, rhs);

  /* Az and Az_rhs */
  group = CCTK_VarIndex("GRHayLMHD::Az");
  rhs = CCTK_VarIndex("GRHayLMHD::Az_rhs");
  ierr += MoLRegisterEvolved(group, rhs);

  /* phitilde and phitilde_rhs */
  group = CCTK_VarIndex("GRHayLMHD::phitilde");
  rhs = CCTK_VarIndex("GRHayLMHD::phitilde_rhs");
  ierr += MoLRegisterEvolved(group, rhs);

  /* ALL OTHER EVOLVED VARIABLES (rho_star,tau,Stilde_x,Stilde_y,Stilde_z) */
  group = CCTK_GroupIndex("GRHayLMHD::grmhd_conservatives");
  rhs = CCTK_GroupIndex("GRHayLMHD::grmhd_conservatives_rhs");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  if(ghl_params->evolve_entropy) {
    group = CCTK_GroupIndex("GRHayLMHD::ent_star");
    rhs = CCTK_GroupIndex("GRHayLMHD::ent_star_rhs");
    ierr += MoLRegisterEvolvedGroup(group, rhs);
  }

  if(ghl_eos->eos_type == ghl_eos_tabulated) {
    group = CCTK_GroupIndex("GRHayLMHD::Ye_star");
    rhs = CCTK_GroupIndex("GRHayLMHD::Ye_star_rhs");
    ierr += MoLRegisterEvolvedGroup(group, rhs);
  }

  if (ierr) CCTK_ERROR("Problems registering with MoL");
  //***********************************************

  //***********************************************
  // Next register ADMBase variables needed by
  //    GRHayLMHD as SaveAndRestore, so that
  //    they are not set to NaN at the start of
  //    each timestep (requiring that they be
  //    e.g., recomputed from BSSN variables
  //    in the BSSN solver, like Baikal or
  //    ML_BSSN)
  ierr += MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex("admbase::lapse"));
  ierr += MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex("admbase::shift"));
  ierr += MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex("admbase::metric"));
  ierr += MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex("admbase::curv"));
  if (ierr) CCTK_ERROR("Problems registering with MoLRegisterSaveAndRestoreGroup");
}
