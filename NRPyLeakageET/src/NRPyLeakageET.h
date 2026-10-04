#ifndef NRPYLEAKAGEET_H_
#define NRPYLEAKAGEET_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "GRHayLib.h"

/* GRHayL's API is double precision. */
#ifndef CCTK_REAL_PRECISION_8
#error "NRPyLeakageET requires Cactus REAL_PRECISION=8 (double)"
#endif

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

// Returns 1 and stores det(gamma_ij) in *gdet_out if gamma_ij is finite and
// positive definite (checked on the computed leading principal minors);
// otherwise returns 0 and leaves *gdet_out untouched. Callers own their
// diagnostic and control flow.
static inline int NRPyLeakageET_spatial_metric_valid(
      const CCTK_REAL gxx, const CCTK_REAL gxy, const CCTK_REAL gxz,
      const CCTK_REAL gyy, const CCTK_REAL gyz, const CCTK_REAL gzz,
      CCTK_REAL *gdet_out) {
  const CCTK_REAL gdet = gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz
                       - gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz;
  if(!robust_isfinite(gxx) || !robust_isfinite(gxy) || !robust_isfinite(gxz) ||
     !robust_isfinite(gyy) || !robust_isfinite(gyz) || !robust_isfinite(gzz) ||
     !robust_isfinite(gdet) || gxx <= 0 || gxx*gyy - gxy*gxy <= 0 || gdet <= 0)
    return 0;
  *gdet_out = gdet;
  return 1;
}

#ifdef __cplusplus
} // extern "C"
#endif

#endif // NRPYLEAKAGEET_H_
