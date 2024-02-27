/*
  Set the symmetries for the GRHayLHD variables
*/

#include "GRHayLHD.h"
#include "Symmetry.h"

void GRHayLHD_InitSymBound(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLHD_InitSymBound;
  DECLARE_CCTK_PARAMETERS;

  if ((cctk_nghostzones[0]<3 || cctk_nghostzones[1]<3 || cctk_nghostzones[2]<3))
    CCTK_VERROR("ERROR: The version of PPM in this thorn requires 3 ghostzones. You only have (%d,%d,%d) ghostzones!",cctk_nghostzones[0],cctk_nghostzones[1],cctk_nghostzones[2]);

  if(cctk_iteration==0) {
    CCTK_VINFO("Setting Symmetry = %s... at iteration = %d",Symmetry,cctk_iteration);

    int sym[3];

    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymGN(cctkGH, sym, "GRHayLHD::grhd_conservatives");
    SetCartSymGN(cctkGH, sym, "GRHayLHD::grhd_velocities");
    if(ghl_params->evolve_entropy)
      SetCartSymGN(cctkGH, sym, "GRHayLHD::ent_star");
    if(ghl_eos->eos_type == ghl_eos_tabulated)
      SetCartSymGN(cctkGH, sym, "GRHayLHD::Ye_star");

    if(CCTK_EQUALS(Symmetry, "equatorial")) {
      sym[2] = -1;
      SetCartSymVN(cctkGH, sym, "GRHayLHD::Stilde_z");
      SetCartSymVN(cctkGH, sym, "GRHayLHD::vz");
    } else if (!CCTK_EQUALS(Symmetry, "none")) {
      CCTK_ERROR("GRHayLHD_initsymbound: Should not be here; picked an impossible symmetry.");
    }
  }
}

