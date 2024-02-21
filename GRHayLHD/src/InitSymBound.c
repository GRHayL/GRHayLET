/*
  Set the symmetries for the GRHayLHD variables
*/

#include "GRHayLHD.h"
#include "Symmetry.h"

void GRHayLHD_InitSymBound(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_GRHayLHD_InitSymBound;
  DECLARE_CCTK_PARAMETERS;

  if ((cctk_nghostzones[0]<3 || cctk_nghostzones[1]<3 || cctk_nghostzones[2]<3))
    CCTK_VERROR("ERROR: The version of PPM in this thorn requires 3 ghostzones. You only have (%d,%d,%d) ghostzones!",cctk_nghostzones[0],cctk_nghostzones[1],cctk_nghostzones[2]);

  if(cctk_iteration==0) {
    CCTK_VINFO("Setting Symmetry = %s... at iteration = %d",Symmetry,cctk_iteration);

    int sym[3];

    if(CCTK_EQUALS(Symmetry,"none")) {
      /* FIRST SET NO SYMMETRY OPTION */
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      SetCartSymGN(cctkGH,sym,"GRHayLHD::grmhd_conservatives");
      SetCartSymGN(cctkGH,sym,"GRHayLHD::grmhd_velocities");
    } else if(CCTK_EQUALS(Symmetry,"equatorial")) {
      /* THEN SET EQUATORIAL SYMMETRY OPTION */
      // Set default to no symmetry, which is correct for scalars and most vectors:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      SetCartSymGN(cctkGH,sym,"GRHayLHD::grmhd_conservatives");
      SetCartSymGN(cctkGH,sym,"GRHayLHD::grmhd_velocities");
    } else {
      CCTK_ERROR("GRHayLHD_initsymbound: Should not be here; picked an impossible symmetry.");
    }
  }
}


