#ifndef GRHAYLMHD_H_
#define GRHAYLMHD_H_

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "GRHayLib.h"

bool IllinoisGRMHD_enforce_outflow(
      const ghl_parameters *params, const ghl_metric_quantities *metric,
      const int sign[3], ghl_primitive_quantities *prims);

enum recon_indices{
      BX_STAGGER, BY_STAGGER, BZ_STAGGER,
      VXR, VYR, VZR, VXL,VYL, VZL, MAXNUMVARS};

// This is used to perturb data for testing. It is a counter-based generator: the value
// depends only on the seed, the global grid index (gi,gj,gk), and the per-variable slot,
// so it does not depend on the number of OpenMP threads, the loop order, or the domain
// decomposition. Returns 1 + perturb*u with u uniform in [0,1).
static inline unsigned long long IllinoisGRMHD_mix64(unsigned long long x) {
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}

static inline CCTK_REAL IllinoisGRMHD_one_plus_pert(
      const CCTK_REAL perturb, const int seed,
      const int gi, const int gj, const int gk, const int slot) {
  unsigned long long h = IllinoisGRMHD_mix64((unsigned long long)(unsigned int)seed);
  h = IllinoisGRMHD_mix64(h ^ (unsigned long long)(unsigned int)gi);
  h = IllinoisGRMHD_mix64(h ^ (unsigned long long)(unsigned int)gj);
  h = IllinoisGRMHD_mix64(h ^ (unsigned long long)(unsigned int)gk);
  h = IllinoisGRMHD_mix64(h ^ (unsigned long long)(unsigned int)slot);
  const CCTK_REAL u = (CCTK_REAL)(h >> 11) * (1.0/9007199254740992.0);
  return 1.0 + perturb*u;
}

// The inner two points of the interpolation function use
// the value of A_in, and the outer two points use A_out.
#define A_out -0.0625
#define A_in  0.5625
//Interpolates to the -1/2 face of point Var
#define COMPUTE_FCVAL(Varm2,Varm1,Var,Varp1) (A_out*(Varm2) + A_in*(Varm1) + A_in*(Var) + A_out*(Varp1))

// Computes 4th-order derivative
#define B_out -1.0/12.0
#define B_in  2.0/3.0
#define COMPUTE_DERIV(Varm2,Varm1,Varp1,Varp2) (B_in*(Varp1 - Varm1) + B_out*(Varp2 - Varm2))

void IllinoisGRMHD_interpolate_metric_to_face(
      const cGH *cctkGH,
      const int i, const int j, const int k,
      const int flux_dirn,
      const CCTK_REAL *restrict lapse,
      const CCTK_REAL *restrict betax,
      const CCTK_REAL *restrict betay,
      const CCTK_REAL *restrict betaz,
      const CCTK_REAL *restrict gxx,
      const CCTK_REAL *restrict gxy,
      const CCTK_REAL *restrict gxz,
      const CCTK_REAL *restrict gyy,
      const CCTK_REAL *restrict gyz,
      const CCTK_REAL *restrict gzz,
      ghl_metric_quantities *restrict metric);

void IllinoisGRMHD_compute_metric_derivs(
      const cGH *cctkGH,
      const int i, const int j, const int k,
      const int flux_dirn,
      const CCTK_REAL dxi,
      const CCTK_REAL *restrict lapse,
      const CCTK_REAL *restrict betax,
      const CCTK_REAL *restrict betay,
      const CCTK_REAL *restrict betaz,
      const CCTK_REAL *restrict gxx,
      const CCTK_REAL *restrict gxy,
      const CCTK_REAL *restrict gxz,
      const CCTK_REAL *restrict gyy,
      const CCTK_REAL *restrict gyz,
      const CCTK_REAL *restrict gzz,
      ghl_metric_quantities *restrict metric_derivs);

void IllinoisGRMHD_set_symmetry_gzs_staggered(
      const cGH *cctkGH,
      const CCTK_REAL *X,
      const CCTK_REAL *Y,
      const CCTK_REAL *Z,
      CCTK_REAL *gridfunc,
      const CCTK_REAL *gridfunc_syms,
      const int stagger_x,  //TODO: unused
      const int stagger_y,  //TODO: unused
      const int stagger_z);

/******** Helper functions for the RHS calculations *************/

void IllinoisGRMHD_reconstruction_loop(
      const cGH *restrict cctkGH,
      const int flux_dir,
      const int num_vars,
      const int *restrict var_indices,
      const CCTK_REAL *pressure,
      const CCTK_REAL *v_flux,
      const CCTK_REAL **in_prims,
      CCTK_REAL **out_prims_r,
      CCTK_REAL **out_prims_l);

// The const are commented out because C does not support implicit typecasting of types when
// they are more than 1 level removed from the top pointer. i.e. I can pass the argument with
// type "CCTK_REAL *" for an argument expecting "const CCTK_REAL *" because this is only 1 level
// down (pointer to CCTK_REAL -> pointer to const CCTK_REAL). It will not do
// pointer to pointer to CCTK_REAL -> pointer to pointer to const CCTK_REAL. I saw comments
// suggesting this may become part of the C23 standard, so I guess you can uncomment this
// in like 10 years.
void IllinoisGRMHD_A_flux_rhs(
      const cGH *restrict cctkGH,
      const int A_dir,
      /*const*/ CCTK_REAL **out_prims_r,
      /*const*/ CCTK_REAL **out_prims_l,
      /*const*/ CCTK_REAL **cmin,
      /*const*/ CCTK_REAL **cmax,
      CCTK_REAL *restrict A_rhs);

/****************************************************************/
#endif // GRHAYLMHD_H_
