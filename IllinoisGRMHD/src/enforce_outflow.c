#include "IllinoisGRMHD.h"
#include <float.h>

static double norm2(const ghl_metric_quantities *metric, const double q[3]) {
  double result = 0.0;
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++) result += metric->gammaDD[i][j]*q[i]*q[j];
  return result;
}

// Enforce all physical faces at once, including edges and corners. sign[d] is
// +1 on an upper face, -1 on a lower face, and zero on an unconstrained axis.
// Find the minimum metric-norm q=v+beta in the coordinate outflow halfspaces,
// then shorten the segment from that feasible point to the requested velocity.
// Both endpoints satisfy the signs, so this preserves them while limiting W.
bool IllinoisGRMHD_enforce_outflow(
    const ghl_parameters *params, const ghl_metric_quantities *metric,
    const int sign[3], ghl_primitive_quantities *prims) {
  double target[3], closest[3] = {0};
  for(int d=0; d<3; d++) {
    if(sign[d]*prims->vU[d] < 0.0) prims->vU[d] = 0.0;
    target[d] = prims->vU[d] + metric->betaU[d];
  }
  const double bound = metric->lapse*metric->lapse*(1.0 - params->inv_sq_max_Lorentz_factor);
  // Leave room for the norm and lapse division to round at the limiting surface.
  const double safe_bound = bound*(1.0 - 16.0*DBL_EPSILON);
  double best = INFINITY;
  // In three dimensions there are at most eight possible active face sets.
  for(int mask=0; mask<8; mask++) {
    int active[3], n=0;
    bool valid = true;
    for(int d=0; d<3; d++) if(mask & (1<<d)) {
      if(sign[d] == 0) valid = false;
      active[n++] = d;
    }
    if(!valid) continue;
    double a[3][4] = {{0}}, q[3] = {0};
    for(int i=0; i<n; i++) {
      for(int j=0; j<n; j++) a[i][j] = metric->gammaUU[active[i]][active[j]];
      a[i][n] = metric->betaU[active[i]];
    }
    // Positive-definite principal submatrices permit elimination without pivoting.
    for(int i=0; i<n; i++) {
      const double pivot = a[i][i];
      if(!(pivot > 0.0) || !isfinite(pivot)) return false;
      for(int j=i; j<=n; j++) a[i][j] /= pivot;
      for(int k=0; k<n; k++) if(k != i) {
        const double factor = a[k][i];
        for(int j=i; j<=n; j++) a[k][j] -= factor*a[i][j];
      }
    }
    for(int i=0; i<3; i++)
      for(int j=0; j<n; j++) q[i] += metric->gammaUU[i][active[j]]*a[j][n];
    // Set the constrained components exactly to avoid roundoff sign violations.
    for(int j=0; j<n; j++) q[active[j]] = metric->betaU[active[j]];
    for(int d=0; d<3; d++)
      if(sign[d]*(q[d]-metric->betaU[d]) < 0.0) valid = false;
    const double value = norm2(metric, q);
    if(valid && value < best) {
      best = value;
      for(int d=0; d<3; d++) closest[d] = q[d];
    }
  }
  if(!isfinite(best) || best > bound) return false; // No allowed outflow state.
  double q[3];
  for(int d=0; d<3; d++) q[d] = target[d];
  if(norm2(metric, q) > safe_bound) {
    // Monotone intersection along the segment from the minimum-norm point.
    // Bisection avoids cancellation when target and closest nearly coincide.
    double low=0.0, high=1.0;
    const double search_bound = fmax(best, safe_bound);
    for(int iter=0; iter<64; iter++) {
      const double mid = 0.5*(low+high);
      for(int d=0; d<3; d++) q[d] = closest[d] + mid*(target[d]-closest[d]);
      if(norm2(metric, q) <= search_bound) low=mid; else high=mid;
    }
    for(int d=0; d<3; d++) q[d] = closest[d] + low*(target[d]-closest[d]);
  }
  for(int d=0; d<3; d++) prims->vU[d] = q[d]-metric->betaU[d];
  const double v2 = norm2(metric, q)/(metric->lapse*metric->lapse);
  if(!isfinite(v2) || v2 > 1.0-params->inv_sq_max_Lorentz_factor) return false;
  prims->u0 = 1.0/(metric->lapse*sqrt(1.0-v2));
  return isfinite(prims->u0);
}
