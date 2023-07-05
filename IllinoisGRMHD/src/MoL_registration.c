//--------------------------------------------------------------------------
// Register with the time stepper
// (MoL thorn, found in arrangements/CactusBase/MoL)
// To understand this, read documentation in arrangements/CactusBase/MoL/doc
//--------------------------------------------------------------------------

#include "IGM.h"
#include "Symmetry.h"

void IllinoisGRMHD_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_IllinoisGRMHD_RegisterVars;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0, var, rhs;

  //***********************************************
  // Register evolution & RHS gridfunction variables

  /* Ax and Ax_rhs */
  var = CCTK_VarIndex("IllinoisGRMHD::Ax");
  rhs = CCTK_VarIndex("IllinoisGRMHD::Ax_rhs");
  ierr += MoLRegisterEvolved(var, rhs);

  /* Ay and Ay_rhs */
  var = CCTK_VarIndex("IllinoisGRMHD::Ay");
  rhs = CCTK_VarIndex("IllinoisGRMHD::Ay_rhs");
  ierr += MoLRegisterEvolved(var, rhs);

  /* Az and Az_rhs */
  var = CCTK_VarIndex("IllinoisGRMHD::Az");
  rhs = CCTK_VarIndex("IllinoisGRMHD::Az_rhs");
  ierr += MoLRegisterEvolved(var, rhs);

  /* phitilde and phitilde_rhs */
  var = CCTK_VarIndex("IllinoisGRMHD::phitilde");
  rhs = CCTK_VarIndex("IllinoisGRMHD::phitilde_rhs");
  ierr += MoLRegisterEvolved(var, rhs);

  /* ALL OTHER EVOLVED VARIABLES (rho_star,tau,Stilde_x,Stilde_y,Stilde_z) */
  var = CCTK_GroupIndex("IllinoisGRMHD::grmhd_conservatives");
  rhs = CCTK_GroupIndex("IllinoisGRMHD::grmhd_conservatives_rhs");
  ierr += MoLRegisterEvolvedGroup(var, rhs);

  if (ierr) CCTK_ERROR("Problems registering with MoL");
  //***********************************************

  //***********************************************
  // Next register ADMBase variables needed by
  //    IllinoisGRMHD as SaveAndRestore, so that
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
