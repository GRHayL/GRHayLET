#include <stdio.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "Symmetry.h"

void GRHayLM1_InitSym(CCTK_ARGUMENTS){
  DECLARE_CCTK_ARGUMENTS_GRHayLM1_set_symmetry
  DECLARE_CCTK_PARAMETERS

  int size = 64;
  char vname[size]; 

  int const sym_scal[3] = {1,1,1};
  int const sym_vecx[3] = {-1,1,1};
  int const sym_vecy[3] = {1,-1,1};
  int const sym_vecz[3] = {1,1,-1};

  for (int i = 0; i < nspecies; ++i){
    snprintf(vname, size, "GRHayLM1::rN[%d]", i);
    SetCartSymVN(cctkGH, sym_scal, vname);

    snprintf(vname, size, "GRHayLM1::rE[%d]", i);
    SetCartSymVN(cctkGH, sym_scal, vname);

    snprintf(vname, size, "GRHayLM1::rFx[%d]", i);
    SetCartSymVN(cctkGH, sym_vecx, vname);

    snprintf(vname, size, "GRHayLM1::rFy[%d]", i);
    SetCartSymVN(cctkGH, sym_vecy, vname);

    snprintf(vname, size, "GRHayLM1::rFz[%d]", i);
    SetCartSymVN(cctkGH, sym_vecz, vname);


  }

}
