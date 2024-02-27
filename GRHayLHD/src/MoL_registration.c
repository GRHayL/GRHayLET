//--------------------------------------------------------------------------
// Register with the time stepper
// (MoL thorn, found in arrangements/CactusBase/MoL)
// To understand this, read documentation in arrangements/CactusBase/MoL/doc
//--------------------------------------------------------------------------

#include "GRHayLHD.h"
#include "Symmetry.h"

void GRHayLHD_RegisterVars(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0, group, rhs;

  //***********************************************
  // Register evolution & RHS gridfunction variables
  group = CCTK_GroupIndex("GRHayLHD::grhd_conservatives");
  rhs = CCTK_GroupIndex("GRHayLHD::grhd_conservatives_rhs");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  if(ghl_params->evolve_entropy) {
    group = CCTK_GroupIndex("GRHayLHD::ent_star");
    rhs = CCTK_GroupIndex("GRHayLHD::ent_star_rhs");
    ierr += MoLRegisterEvolvedGroup(group, rhs);
  }

  if(ghl_eos->eos_type == ghl_eos_tabulated) {
    group = CCTK_GroupIndex("GRHayLHD::Ye_star");
    rhs = CCTK_GroupIndex("GRHayLHD::Ye_star_rhs");
    ierr += MoLRegisterEvolvedGroup(group, rhs);
  }

  if (ierr) CCTK_ERROR("Problems registering with MoL");
  //***********************************************

  //***********************************************
  // Next register ADMBase variables needed by
  //    GRHayLHD as SaveAndRestore, so that
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
