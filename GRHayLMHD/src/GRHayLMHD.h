#ifndef GRHAYLMHD_H_
#define GRHAYLMHD_H_

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "GRHayLib.h"

enum recon_indices{
      BX_STAGGER, BY_STAGGER, BZ_STAGGER,
      VXR, VYR, VZR, VXL,VYL, VZL, MAXNUMVARS};

// This is used to perturb data for testing
#define one_plus_pert(perturb) (1 + (perturb*(double)rand() / RAND_MAX))

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

void GRHayLMHD_interpolate_metric_to_face(
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

void GRHayLMHD_compute_metric_derivs(
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

void GRHayLMHD_set_symmetry_gzs_staggered(
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

void GRHayLMHD_reconstruction_loop(const cGH *restrict cctkGH, const int flux_dir, const int num_vars,
                         const int *restrict var_indices,
                         const double *rho_b,
                         const double *pressure,
                         const double *v_flux,
                         const CCTK_REAL **in_prims,
                         CCTK_REAL **out_prims_r,
                         CCTK_REAL **out_prims_l);

void GRHayLMHD_calculate_MHD_dirn_rhs(
      const cGH *restrict cctkGH,
      const int flux_dir,
      const CCTK_REAL *restrict dX,
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
      const CCTK_REAL *restrict rho_b,
      const CCTK_REAL *restrict pressure,
      const CCTK_REAL *restrict eps,
      const CCTK_REAL *restrict ent,
      const CCTK_REAL *restrict Ye,
      const CCTK_REAL *restrict temp,
      const CCTK_REAL *vx,
      const CCTK_REAL *vy,
      const CCTK_REAL *vz,
      const CCTK_REAL **B_center,
      const CCTK_REAL *restrict B_stagger,
      CCTK_REAL **vel_r,
      CCTK_REAL **vel_l,
      CCTK_REAL *restrict cmin,
      CCTK_REAL *restrict cmax,
      CCTK_REAL *restrict rho_star_flux,
      CCTK_REAL *restrict tau_flux,
      CCTK_REAL *restrict Stildex_flux,
      CCTK_REAL *restrict Stildey_flux,
      CCTK_REAL *restrict Stildez_flux,
      CCTK_REAL *restrict ent_star_flux,
      CCTK_REAL *restrict Ye_star_flux,
      CCTK_REAL *restrict rho_star_rhs,
      CCTK_REAL *restrict tau_rhs,
      CCTK_REAL *restrict Stildex_rhs,
      CCTK_REAL *restrict Stildey_rhs,
      CCTK_REAL *restrict Stildez_rhs,
      CCTK_REAL *restrict ent_star_rhs,
      CCTK_REAL *restrict Ye_star_rhs);

// The const are commented out because C does not support implicit typecasting of types when
// they are more than 1 level removed from the top pointer. i.e. I can pass the argument with
// type "CCTK_REAL *" for an argument expecting "const CCTK_REAL *" because this is only 1 level
// down (pointer to CCTK_REAL -> pointer to const CCTK_REAL). It will not do
// pointer to pointer to CCTK_REAL -> pointer to pointer to const CCTK_REAL. I saw comments
// suggesting this may become part of the C23 standard, so I guess you can uncomment this
// in like 10 years.
void GRHayLMHD_A_flux_rhs(
      const cGH *restrict cctkGH,
      const int A_dir,
      /*const*/ CCTK_REAL **out_prims_r,
      /*const*/ CCTK_REAL **out_prims_l,
      /*const*/ CCTK_REAL **cmin,
      /*const*/ CCTK_REAL **cmax,
      CCTK_REAL *restrict A_rhs);

/****************************************************************/
#endif // GRHAYLMHD_H_
