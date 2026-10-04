#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "Symmetry.h"

void NRPyLeakageET_InitSym(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int sym[3] = {1,1,1};
  if(SetCartSymGN(cctkGH,sym,"NRPyLeakageET::NRPyLeakageET_opacities") < 0 ||
     SetCartSymGN(cctkGH,sym,"NRPyLeakageET::NRPyLeakageET_optical_depths") < 0)
    CCTK_ERROR("Could not register leakage scalar parity");
}

void NRPyLeakageET_SelectDriverBoundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  const char *groups[2] = {"NRPyLeakageET::NRPyLeakageET_optical_depths",
                          "NRPyLeakageET::NRPyLeakageET_opacities"};
  for(int g=0;g<2;g++)
    if(Driver_SelectGroupForBC(cctkGH,CCTK_ALL_FACES,1,-1,groups[g],"none") < 0)
      CCTK_VERROR("Could not select driver boundaries for %s",groups[g]);
}
