#ifndef NRPYLEAKAGEET_H_
#define NRPYLEAKAGEET_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "GRHayLib.h"

#ifdef __cplusplus
extern "C" {
#endif

int  NRPyLeakageET_ProcessOwnsData();
void NRPyLeakageET_optical_depths_initialize_to_zero(CCTK_ARGUMENTS);
void NRPyLeakageET_copy_opacities_and_optical_depths_to_previous_time_levels(CCTK_ARGUMENTS);
void NRPyLeakageET_compute_optical_depth_change(CCTK_ARGUMENTS, const int it);
void NRPyLeakageET_CopyOpticalDepthsToAux(CCTK_ARGUMENTS);
void NRPyLeakageET_copy_optical_depths_from_previous_time_level(CCTK_ARGUMENTS);

void NRPyLeakageET_compute_neutrino_opacities(CCTK_ARGUMENTS);
void NRPyLeakageET_compute_neutrino_luminosities(CCTK_ARGUMENTS);
void NRPyLeakageET_compute_neutrino_opacities_and_add_source_terms_to_MHD_rhss(CCTK_ARGUMENTS);
void NRPyLeakageET_optical_depths_PathOfLeastResistance(CCTK_ARGUMENTS);

static inline int NRPyLeakageET_opacities_finite(const ghl_neutrino_opacities *kappa) {
  for(int m=0;m<2;m++)
    if(!robust_isfinite(kappa->nue[m]) || !robust_isfinite(kappa->anue[m]) || !robust_isfinite(kappa->nux[m]))
      return 0;
  return 1;
}

#ifdef __cplusplus
} // extern "C"
#endif

#endif // NRPYLEAKAGEET_H_
