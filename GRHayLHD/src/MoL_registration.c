//--------------------------------------------------------------------------
// Register with the time stepper
// (MoL thorn, found in arrangements/CactusBase/MoL)
// To understand this, read documentation in arrangements/CactusBase/MoL/doc
//--------------------------------------------------------------------------

#include "GRHayLHD.h"
#include "Symmetry.h"

void GRHayLHD_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_GRHayLHD_RegisterVars;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0, var, rhs;

  //***********************************************
  // Register evolution & RHS gridfunction variables
  var = CCTK_GroupIndex("GRHayLHD::grmhd_conservatives");
  rhs = CCTK_GroupIndex("GRHayLHD::grmhd_conservatives_rhs");
  ierr += MoLRegisterEvolvedGroup(var, rhs);

  ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex("HydroBase::rho"));
  ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex("HydroBase::press"));
  ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex("HydroBase::eps"));

  if(ghl_params->evolve_entropy) {
    var = CCTK_GroupIndex("GRHayLHD::ent_star");
    rhs = CCTK_GroupIndex("GRHayLHD::ent_star_rhs");
    ierr += MoLRegisterEvolvedGroup(var, rhs);
    ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex("HydroBase::entropy"));
  }

  if(ghl_eos->eos_type == ghl_eos_tabulated) {
    var = CCTK_GroupIndex("GRHayLHD::Ye_star");
    rhs = CCTK_GroupIndex("GRHayLHD::Ye_star_rhs");
    ierr += MoLRegisterEvolvedGroup(var, rhs);
    ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex("HydroBase::Y_e"));
    ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex("HydroBase::temperature"));
  }

  if(stress_energy_at_RHS) {
    MoLRegisterConstrainedGroup(CCTK_GroupIndex("TmunuBase::stress_energy_scalar"));
    MoLRegisterConstrainedGroup(CCTK_GroupIndex("TmunuBase::stress_energy_vector"));
    MoLRegisterConstrainedGroup(CCTK_GroupIndex("TmunuBase::stress_energy_tensor"));
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
  ierr += MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex("ADMBase::lapse"));
  ierr += MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex("ADMBase::shift"));
  ierr += MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex("ADMBase::metric"));
  ierr += MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex("ADMBase::curv"));
  if (ierr) CCTK_ERROR("Problems registering with MoLRegisterSaveAndRestoreGroup");
}
