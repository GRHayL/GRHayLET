#include <string.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

void GRHayLM1_InitializeToZero(CCTK_ARGUMENTS){
	DECLARE_CCTK_PARAMETERS;
	DECLARE_CCTK_ARGUMENTS;

	size_t GFsize = cctkGH->cctk_ash[0]*cctkGH->cctk_ash[1]*cctkGH->cctk_ash[2];
	size_t size = GFsize*GRHayLM1_nspecies*sizeof(CCTK_REAL);
	//size_t isize = GFsize*sizeof(CCTK_INT);

	//Set M1 variables to zero.
	memset(GRHayLM1_rN, 0, size);
	memset(GRHayLM1_rE, 0, size);
	memset(GRHayLM1_rFx, 0, size);
	memset(GRHayLM1_rFy, 0, size);
	memset(GRHayLM1_rFz, 0, size);
	//memset(GRHayLM1_mask, 0, isize);

}
