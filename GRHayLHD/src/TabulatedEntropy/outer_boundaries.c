/*******************************************************
 * Outer boundaries are handled as follows:
 * (-1) Update RHS quantities, leave RHS quantities zero on all outer ghostzones
 *      (including outer AMR refinement, processor, and outer boundaries)
 * ( 0) Let MoL update all evolution variables
 * ( 1) Apply outer boundary conditions (BCs) on A_{\mu}
 * ( 2) Compute B^i from A_i everywhere, synchronize B^i
 * ( 3) Call con2prim to get primitives on interior pts
 * ( 4) Apply outer BCs on {P,rho,vx,vy,vz}.
 * ( 5) (optional) set conservatives on outer boundary.
 *******************************************************/

#include "GRHayLHD.h"

void GRHayLHD_tabulated_entropy_enforce_primitive_limits_and_compute_conservs(const cGH* cctkGH, const int index, ghl_primitive_quantities *restrict prims);

/*******************************************************
 * Apply outer boundary conditions on {P,rho,vx,vy,vz}
 * It is better to apply BCs on primitives than conservs,
 * because small errors in conservs can be greatly
 * amplified in con2prim, sometimes leading to unphysical
 * primitives & unnecessary fixes.
 *******************************************************/
void GRHayLHD_tabulated_entropy_outer_boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLHD_tabulated_entropy_outer_boundaries;
  DECLARE_CCTK_PARAMETERS;

  if(CCTK_EQUALS(Matter_BC, "frozen")) return;

  const bool do_outflow = CCTK_EQUALS(Matter_BC, "outflow");

  const bool Symmetry_none = CCTK_EQUALS(Symmetry, "none") ? true : false;

  // Don't apply approximate outer boundary conditions on initial data, which
  // should be defined everywhere, or on levels != [coarsest level].
  if(cctk_iteration==0 || GetRefinementLevel(cctkGH)!=0) return;

  if(cctk_nghostzones[0]!=cctk_nghostzones[1] || cctk_nghostzones[0]!=cctk_nghostzones[2])
    CCTK_ERROR("ERROR: GRHayLHD outer BC driver does not support unequal number of ghostzones in different directions!");
  for(int which_bdry_pt=0;which_bdry_pt<cctk_nghostzones[0];which_bdry_pt++) {

    /* XMIN & XMAX */
    // i=imax=outer boundary
    if(cctk_bbox[1]) {
      const int imax = cctk_lsh[0] - cctk_nghostzones[0] + which_bdry_pt;
#pragma omp parallel for
      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int j=0; j<cctk_lsh[1]; j++) {
          const int index = CCTK_GFINDEX3D(cctkGH,imax, j, k);
          const int indm1 = CCTK_GFINDEX3D(cctkGH,imax-1, j, k);

          ghl_primitive_quantities prims;
          prims.BU[0] = prims.BU[1] = prims.BU[2] = 0.0;
          prims.rho         = rho[indm1];
          prims.press       = press[indm1];
          prims.vU[0]       = (do_outflow && vx[indm1] < 0.0) ? 0 : vx[indm1];
          prims.vU[1]       = vy[indm1];
          prims.vU[2]       = vz[indm1];
          prims.entropy     = entropy[indm1];
          prims.Y_e         = Y_e[indm1];
          prims.temperature = temperature[indm1];

          GRHayLHD_tabulated_entropy_enforce_primitive_limits_and_compute_conservs(cctkGH, index, &prims);
        }
      }
    }
    // i=imin=outer boundary
    if(cctk_bbox[0]) {
      const int imin = cctk_nghostzones[0] - which_bdry_pt - 1;
#pragma omp parallel for
      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int j=0; j<cctk_lsh[1]; j++) {
          const int index = CCTK_GFINDEX3D(cctkGH, imin, j, k);
          const int indp1 = CCTK_GFINDEX3D(cctkGH, imin+1, j, k);

          ghl_primitive_quantities prims;
          prims.BU[0] = prims.BU[1] = prims.BU[2] = 0.0;
          prims.rho         = rho[indp1];
          prims.press       = press[indp1];
          prims.vU[0]       = (do_outflow && vx[indp1] > 0.0) ? 0 : vx[indp1];
          prims.vU[1]       = vy[indp1];
          prims.vU[2]       = vz[indp1];
          prims.entropy     = entropy[indp1];
          prims.Y_e         = Y_e[indp1];
          prims.temperature = temperature[indp1];

          GRHayLHD_tabulated_entropy_enforce_primitive_limits_and_compute_conservs(cctkGH, index, &prims);
        }
      }
    }

    /* YMIN & YMAX */
    // j=jmax=outer boundary
    if(cctk_bbox[3]) {
      const int jmax = cctk_lsh[1] - cctk_nghostzones[1] + which_bdry_pt;
#pragma omp parallel for
      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int i=0; i<cctk_lsh[0]; i++) {
          const int index = CCTK_GFINDEX3D(cctkGH, i, jmax, k);
          const int indm1 = CCTK_GFINDEX3D(cctkGH, i, jmax-1, k);

          ghl_primitive_quantities prims;
          prims.BU[0] = prims.BU[1] = prims.BU[2] = 0.0;
          prims.rho         = rho[indm1];
          prims.press       = press[indm1];
          prims.vU[0]       = vx[indm1];
          prims.vU[1]       = (do_outflow && vy[indm1] < 0.0) ? 0 : vy[indm1];
          prims.vU[2]       = vz[indm1];
          prims.entropy     = entropy[indm1];
          prims.Y_e         = Y_e[indm1];
          prims.temperature = temperature[indm1];

          GRHayLHD_tabulated_entropy_enforce_primitive_limits_and_compute_conservs(cctkGH, index, &prims);
        }
      }
    }
    // j=jmin=outer boundary
    if(cctk_bbox[2]) {
      const int jmin = cctk_nghostzones[1] - which_bdry_pt - 1;
#pragma omp parallel for
      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int i=0; i<cctk_lsh[0]; i++) {
          const int index = CCTK_GFINDEX3D(cctkGH, i, jmin, k);
          const int indp1 = CCTK_GFINDEX3D(cctkGH, i, jmin+1, k);

          ghl_primitive_quantities prims;
          prims.BU[0] = prims.BU[1] = prims.BU[2] = 0.0;
          prims.rho         = rho[indp1];
          prims.press       = press[indp1];
          prims.vU[0]       = vx[indp1];
          prims.vU[1]       = (do_outflow && vy[indp1] > 0.0) ? 0 : vy[indp1];
          prims.vU[2]       = vz[indp1];
          prims.entropy     = entropy[indp1];
          prims.Y_e         = Y_e[indp1];
          prims.temperature = temperature[indp1];

          GRHayLHD_tabulated_entropy_enforce_primitive_limits_and_compute_conservs(cctkGH, index, &prims);
        }
      }
    }

    /* ZMIN & ZMAX */
    // k=kmax=outer boundary
    if(cctk_bbox[5]) {
      const int kmax = cctk_lsh[2] - cctk_nghostzones[2] + which_bdry_pt;
#pragma omp parallel for
      for(int j=0; j<cctk_lsh[1]; j++) {
        for(int i=0; i<cctk_lsh[0]; i++) {
          const int index = CCTK_GFINDEX3D(cctkGH, i, j, kmax);
          const int indm1 = CCTK_GFINDEX3D(cctkGH, i, j, kmax-1);

          ghl_primitive_quantities prims;
          prims.BU[0] = prims.BU[1] = prims.BU[2] = 0.0;
          prims.rho         = rho[indm1];
          prims.press       = press[indm1];
          prims.vU[0]       = vx[indm1];
          prims.vU[1]       = vy[indm1];
          prims.vU[2]       = (do_outflow && vz[indm1] < 0.0) ? 0 : vz[indm1];
          prims.entropy     = entropy[indm1];
          prims.Y_e         = Y_e[indm1];
          prims.temperature = temperature[indm1];

          GRHayLHD_tabulated_entropy_enforce_primitive_limits_and_compute_conservs(cctkGH, index, &prims);
        }
      }
    }
    // k=kmin=outer boundary
    if((cctk_bbox[4]) && Symmetry_none) {
      const int kmin = cctk_nghostzones[2] - which_bdry_pt - 1;
#pragma omp parallel for
      for(int j=0; j<cctk_lsh[1]; j++) {
        for(int i=0; i<cctk_lsh[0]; i++) {
          const int index = CCTK_GFINDEX3D(cctkGH, i, j, kmin);
          const int indp1 = CCTK_GFINDEX3D(cctkGH, i, j, kmin+1);

          ghl_primitive_quantities prims;
          prims.BU[0] = prims.BU[1] = prims.BU[2] = 0.0;
          prims.rho         = rho[indp1];
          prims.press       = press[indp1];
          prims.vU[0]       = vx[indp1];
          prims.vU[1]       = vy[indp1];
          prims.vU[2]       = (do_outflow && vz[indp1] > 0.0) ? 0 : vz[indp1];
          prims.entropy     = entropy[indp1];
          prims.Y_e         = Y_e[indp1];
          prims.temperature = temperature[indp1];

          GRHayLHD_tabulated_entropy_enforce_primitive_limits_and_compute_conservs(cctkGH, index, &prims);
        }
      }
    }
  }
}

void GRHayLHD_tabulated_entropy_enforce_primitive_limits_and_compute_conservs(const cGH* cctkGH, const int index, ghl_primitive_quantities *restrict prims) {
  // We cheat here by using the argument list of the scheduled function
  // instead of explicitly passing all these grid functions.
  DECLARE_CCTK_ARGUMENTS_GRHayLHD_tabulated_entropy_outer_boundaries;

  ghl_metric_quantities ADM_metric;
  ghl_enforce_detgtij_and_initialize_ADM_metric(
        alp[index],
        betax[index], betay[index], betaz[index],
        gxx[index], gxy[index], gxz[index],
        gyy[index], gyz[index], gzz[index],
        &ADM_metric);

  ghl_ADM_aux_quantities metric_aux;
  ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

  bool speed_limited;
  ghl_error_codes_t error = ghl_enforce_primitive_limits_and_compute_u0(
        ghl_params, ghl_eos, &ADM_metric, prims, &speed_limited);
  if(error)
    ghl_read_error_codes(error);

  ghl_conservative_quantities cons;
  ghl_compute_conservs(
        &ADM_metric, &metric_aux, prims, &cons);

  rho[index]         = prims->rho;
  press[index]       = prims->press;
  eps[index]         = prims->eps;
  u0[index]          = prims->u0;
  vx[index]          = prims->vU[0];
  vy[index]          = prims->vU[1];
  vz[index]          = prims->vU[2];
  entropy[index]     = prims->entropy;
  Y_e[index]         = prims->Y_e;
  temperature[index] = prims->temperature;

  rho_star[index] = cons.rho;
  tau[index]      = cons.tau;
  Stildex[index]  = cons.SD[0];
  Stildey[index]  = cons.SD[1];
  Stildez[index]  = cons.SD[2];
  ent_star[index] = cons.entropy;
  Ye_star[index]  = cons.Y_e;
}
