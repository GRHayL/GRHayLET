#include "GRHayLHD.h"

void GRHayLHD_no_tabulated(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;

  if(CCTK_Equals(EOS_type, "Tabulated"))
    CCTK_ERROR("Tabulated support is not currently available in the release version of GRHayLHD.");
}
