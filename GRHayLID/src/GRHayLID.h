#ifndef GRHAYLID_H_
#define GRHAYLID_H_

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "GRHayLib.h"
#include <math.h>
#include <float.h>

/* GRHayLib also passes parameter arrays to double-pointer APIs. */
#ifndef CCTK_REAL_PRECISION_8
#error "GRHayLID/GRHayLib require Cactus REAL_PRECISION=8 (double)"
#endif

static inline void GRHayLID_require_storage(const cGH *gh, const char *group) {
  if(CCTK_QueryGroupStorage(gh, group) <= 0)
    CCTK_VERROR("GRHayLID requires active storage for %s", group);
}

static inline void GRHayLID_check_flat_metric(
    const double xx, const double xy, const double xz,
    const double yy, const double yz, const double zz) {
  const double tol = 64.0*DBL_EPSILON;
  if(!isfinite(xx) || !isfinite(xy) || !isfinite(xz) ||
     !isfinite(yy) || !isfinite(yz) || !isfinite(zz) ||
     fabs(xx-1.0) > tol || fabs(yy-1.0) > tol || fabs(zz-1.0) > tol ||
     fabs(xy) > tol || fabs(xz) > tol || fabs(yz) > tol)
    CCTK_ERROR("GRHayLID test data requires a Cartesian identity spatial metric supplied by ADMBase");
}

static inline void GRHayLID_check_table_state(
    const double density, const double ye, const double temp, const char *region) {
  if(!isfinite(density) || !isfinite(ye) || !isfinite(temp) ||
     density <= 0.0 || temp <= 0.0 ||
     density < ghl_eos->rho_min || density > ghl_eos->rho_max ||
     ye < ghl_eos->Y_e_min || ye > ghl_eos->Y_e_max ||
     temp < ghl_eos->T_min || temp > ghl_eos->T_max)
    CCTK_VERROR("Invalid %s EOS state: rho=%g, Ye=%g, T=%g (outside effective bounds)",
                region, density, ye, temp);
}

#define CHECK_PARAMETER(par) if(par==-1) CCTK_VERROR("Please set %s::%s in your parfile",CCTK_THORNSTRING,#par);

#endif // GRHAYLID_H_
