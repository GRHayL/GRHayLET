/*******************************************************
 * Outer boundaries are handled as follows:
 * (-1) Update RHS quantities, leave RHS quantities zero on all outer ghostzones (including outer AMR refinement, processor, and outer boundaries)
 * ( 0) Let MoL update all evolution variables
 * ( 1) Apply outer boundary conditions (BCs) on A_{\mu}
 * ( 2) Compute B^i from A_i everywhere, synchronize B^i
 * ( 3) Call con2prim to get primitives on interior pts
 * ( 4) Apply outer BCs on {P,rho_b,vx,vy,vz}.
 * ( 5) (optional) set conservatives on outer boundary.
 *******************************************************/

#include "GRHayLMHD.h"

void GRHayLMHD_enforce_primitive_limits_and_compute_conservs(const cGH* cctkGH, const int index, ghl_primitive_quantities *restrict prims);

/*********************************************
 * Apply outer boundary conditions on A_{\mu}
 ********************************************/
void GRHayLMHD_outer_boundaries_on_A_mu(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLMHD_outer_boundaries_on_A_mu;
  DECLARE_CCTK_PARAMETERS;

  if(CCTK_EQUALS(EM_BC,"frozen")) return;

  bool Symmetry_none=false; if(CCTK_EQUALS(Symmetry,"none")) Symmetry_none=true;

  const int levelnumber = GetRefinementLevel(cctkGH);

  // Don't apply approximate outer boundary conditions on initial data, which should be defined everywhere, or on levels != [coarsest level].
  if(cctk_iteration==0 || levelnumber!=0) return;

  if(cctk_nghostzones[0]!=cctk_nghostzones[1] || cctk_nghostzones[0]!=cctk_nghostzones[2])
    CCTK_VERROR("ERROR: GRHayLMHD outer BC driver does not support unequal number of ghostzones in different directions!");
  for(int which_bdry_pt=0; which_bdry_pt<cctk_nghostzones[0]; which_bdry_pt++) {

//TODO: These are in blocks of 4 already; this is screaming for some SIMD
    // for cctk_nghostzones==3, max indices go {cctk_lsh-3,cctk_lsh-2,cctk_lsh-1}; outer bdry pt is at cctk_lsh-1
    // for cctk_nghostzones==3, minimum indices go {2,1,0}
    if(cctk_bbox[1]) {
      const int imax=cctk_lsh[0]-cctk_nghostzones[0]+which_bdry_pt;
#pragma omp parallel for
      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int j=0; j<cctk_lsh[1]; j++) {
          const int indm2 = CCTK_GFINDEX3D(cctkGH,imax-2, j, k);
          const int indm1 = CCTK_GFINDEX3D(cctkGH,imax-1, j, k);
          const int index = CCTK_GFINDEX3D(cctkGH,imax, j, k);

          phitilde[index] = 2.0 * phitilde[indm1] - phitilde[indm2];
          Ax[index] = 2.0 * Ax[indm1] - Ax[indm2];
          Ay[index] = 2.0 * Ay[indm1] - Ay[indm2];
          Az[index] = 2.0 * Az[indm1] - Az[indm2];
        }
      }
    }
    if(cctk_bbox[3]) {
      const int jmax=cctk_lsh[1]-cctk_nghostzones[1]+which_bdry_pt;
#pragma omp parallel for
      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int i=0; i<cctk_lsh[0]; i++) {
          const int indm2 = CCTK_GFINDEX3D(cctkGH, i, jmax-2, k);
          const int indm1 = CCTK_GFINDEX3D(cctkGH, i, jmax-1, k);
          const int index = CCTK_GFINDEX3D(cctkGH, i, jmax, k);

          phitilde[index] = 2.0 * phitilde[indm1] - phitilde[indm2];
          Ax[index] = 2.0 * Ax[indm1] - Ax[indm2];
          Ay[index] = 2.0 * Ay[indm1] - Ay[indm2];
          Az[index] = 2.0 * Az[indm1] - Az[indm2];
        }
      }
    }
    if(cctk_bbox[5]) {
      const int kmax=cctk_lsh[2]-cctk_nghostzones[2]+which_bdry_pt;
#pragma omp parallel for
      for(int j=0; j<cctk_lsh[1]; j++) {
        for(int i=0; i<cctk_lsh[0]; i++) {
          const int indm2 = CCTK_GFINDEX3D(cctkGH, i, j, kmax-2);
          const int indm1 = CCTK_GFINDEX3D(cctkGH, i, j, kmax-1);
          const int index = CCTK_GFINDEX3D(cctkGH, i, j, kmax);

          phitilde[index] = 2.0 * phitilde[indm1] - phitilde[indm2];
          Ax[index] = 2.0 * Ax[indm1] - Ax[indm2];
          Ay[index] = 2.0 * Ay[indm1] - Ay[indm2];
          Az[index] = 2.0 * Az[indm1] - Az[indm2];
        }
      }
    }

    if(cctk_bbox[0]) {
      const int imin=cctk_nghostzones[0]-which_bdry_pt-1;
#pragma omp parallel for
      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int j=0; j<cctk_lsh[1]; j++) {
          const int index = CCTK_GFINDEX3D(cctkGH, imin, j, k);
          const int indp1 = CCTK_GFINDEX3D(cctkGH, imin+1, j, k);
          const int indp2 = CCTK_GFINDEX3D(cctkGH, imin+2, j, k);

          phitilde[index] = 2.0 * phitilde[indp1] - phitilde[indp2];
          Ax[index] = 2.0 * Ax[indp1] - Ax[indp2];
          Ay[index] = 2.0 * Ay[indp1] - Ay[indp2];
          Az[index] = 2.0 * Az[indp1] - Az[indp2];
        }
      }
    }
    if(cctk_bbox[2]) {
      const int jmin=cctk_nghostzones[1]-which_bdry_pt-1;
#pragma omp parallel for
      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int i=0; i<cctk_lsh[0]; i++) {
          const int index = CCTK_GFINDEX3D(cctkGH, i, jmin, k);
          const int indp1 = CCTK_GFINDEX3D(cctkGH, i, jmin+1, k);
          const int indp2 = CCTK_GFINDEX3D(cctkGH, i, jmin+2, k);

          phitilde[index] = 2.0 * phitilde[indp1] - phitilde[indp2];
          Ax[index] = 2.0 * Ax[indp1] - Ax[indp2];
          Ay[index] = 2.0 * Ay[indp1] - Ay[indp2];
          Az[index] = 2.0 * Az[indp1] - Az[indp2];
        }
      }
    }
    if((cctk_bbox[4]) && Symmetry_none) {
      const int kmin=cctk_nghostzones[2]-which_bdry_pt-1;
#pragma omp parallel for
      for(int j=0; j<cctk_lsh[1]; j++) {
        for(int i=0; i<cctk_lsh[0]; i++) {
          const int index = CCTK_GFINDEX3D(cctkGH, i, j, kmin);
          const int indp1 = CCTK_GFINDEX3D(cctkGH, i, j, kmin+1);
          const int indp2 = CCTK_GFINDEX3D(cctkGH, i, j, kmin+2);

          phitilde[index] = 2.0 * phitilde[indp1] - phitilde[indp2];
          Ax[index] = 2.0 * Ax[indp1] - Ax[indp2];
          Ay[index] = 2.0 * Ay[indp1] - Ay[indp2];
          Az[index] = 2.0 * Az[indp1] - Az[indp2];
        }
      }
    }
  }
}

/*******************************************************
 * Apply outer boundary conditions on {P,rho_b,vx,vy,vz}
 * It is better to apply BCs on primitives than conservs,
 * because small errors in conservs can be greatly
 * amplified in con2prim, sometimes leading to unphysical
 * primitives & unnecessary fixes.
 *******************************************************/
void GRHayLMHD_outer_boundaries_on_P_rho_b_vx_vy_vz(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLMHD_outer_boundaries_on_P_rho_b_vx_vy_vz;
  DECLARE_CCTK_PARAMETERS;

  if(CCTK_EQUALS(Matter_BC,"frozen")) return;

  const bool Symmetry_none = CCTK_EQUALS(Symmetry,"none") ? true : false;

  // Don't apply approximate outer boundary conditions on initial data, which should be defined everywhere, or on levels != [coarsest level].
  if(cctk_iteration==0 || GetRefinementLevel(cctkGH)!=0) return;

  const double poison = 0.0/0.0;
  double dummy1, dummy2, dummy3;

  if(cctk_nghostzones[0]!=cctk_nghostzones[1] || cctk_nghostzones[0]!=cctk_nghostzones[2])
    CCTK_VERROR("ERROR: GRHayLMHD outer BC driver does not support unequal number of ghostzones in different directions!");
  for(int which_bdry_pt=0;which_bdry_pt<cctk_nghostzones[0];which_bdry_pt++) {

    /* XMIN & XMAX */
    // i=imax=outer boundary
    if(cctk_bbox[1]) {
      const int imax=cctk_lsh[0]-cctk_nghostzones[0]+which_bdry_pt;
#pragma omp parallel for
      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int j=0; j<cctk_lsh[1]; j++) {
          const int index = CCTK_GFINDEX3D(cctkGH,imax, j, k);
          const int indm1 = CCTK_GFINDEX3D(cctkGH,imax-1, j, k);

          const double vtmp = vx[indm1] < 0.0 ? 0 : vx[indm1];
          ghl_primitive_quantities prims;
          ghl_initialize_primitives(
                rho_b[indm1], pressure[indm1], eps[indm1],
                vtmp, vy[indm1], vz[indm1],
                Bx_center[index], By_center[index], Bz_center[index],
                poison, poison, poison, &prims);

          GRHayLMHD_enforce_primitive_limits_and_compute_conservs(cctkGH, index, &prims);

          ghl_return_primitives(
                &prims,
                &rho_b[index], &pressure[index], &eps[index],
                &vx[index], &vy[index], &vz[index],
                &Bx_center[index], &By_center[index], &Bz_center[index],
                &dummy1, &dummy2, &dummy3);
        }
      }
    }
    // i=imin=outer boundary
    if(cctk_bbox[0]) {
      const int imin=cctk_nghostzones[0]-which_bdry_pt-1;
#pragma omp parallel for
      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int j=0; j<cctk_lsh[1]; j++) {
          const int index = CCTK_GFINDEX3D(cctkGH, imin, j, k);
          const int indp1 = CCTK_GFINDEX3D(cctkGH, imin+1, j, k);

          const double vtmp = vx[indp1] > 0.0 ? 0 : vx[indp1];
          ghl_primitive_quantities prims;
          ghl_initialize_primitives(
                rho_b[indp1], pressure[indp1], eps[indp1],
                vtmp, vy[indp1], vz[indp1],
                Bx_center[index], By_center[index], Bz_center[index],
                poison, poison, poison, &prims);

          GRHayLMHD_enforce_primitive_limits_and_compute_conservs(cctkGH, index, &prims);

          ghl_return_primitives(
                &prims,
                &rho_b[index], &pressure[index], &eps[index],
                &vx[index], &vy[index], &vz[index],
                &Bx_center[index], &By_center[index], &Bz_center[index],
                &dummy1, &dummy2, &dummy3);
        }
      }
    }

    /* YMIN & YMAX */
    // j=jmax=outer boundary
    if(cctk_bbox[3]) {
      const int jmax=cctk_lsh[1]-cctk_nghostzones[1]+which_bdry_pt;
#pragma omp parallel for
      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int i=0; i<cctk_lsh[0]; i++) {
          const int index = CCTK_GFINDEX3D(cctkGH, i, jmax, k);
          const int indm1 = CCTK_GFINDEX3D(cctkGH, i, jmax-1, k);

          const double vtmp = vy[indm1] < 0.0 ? 0 : vy[indm1];
          ghl_primitive_quantities prims;
          ghl_initialize_primitives(
                rho_b[indm1], pressure[indm1], eps[indm1],
                vx[indm1], vtmp, vz[indm1],
                Bx_center[index], By_center[index], Bz_center[index],
                poison, poison, poison, &prims);

          GRHayLMHD_enforce_primitive_limits_and_compute_conservs(cctkGH, index, &prims);

          ghl_return_primitives(
                &prims,
                &rho_b[index], &pressure[index], &eps[index],
                &vx[index], &vy[index], &vz[index],
                &Bx_center[index], &By_center[index], &Bz_center[index],
                &dummy1, &dummy2, &dummy3);
        }
      }
    }
    // j=jmin=outer boundary
    if(cctk_bbox[2]) {
      const int jmin=cctk_nghostzones[1]-which_bdry_pt-1;
#pragma omp parallel for
      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int i=0; i<cctk_lsh[0]; i++) {
          const int index = CCTK_GFINDEX3D(cctkGH, i, jmin, k);
          const int indp1 = CCTK_GFINDEX3D(cctkGH, i, jmin+1, k);

          const double vtmp = vy[indp1] > 0.0 ? 0 : vy[indp1];
          ghl_primitive_quantities prims;
          ghl_initialize_primitives(
                rho_b[indp1], pressure[indp1], eps[indp1],
                vx[indp1], vtmp, vz[indp1],
                Bx_center[index], By_center[index], Bz_center[index],
                poison, poison, poison, &prims);

          GRHayLMHD_enforce_primitive_limits_and_compute_conservs(cctkGH, index, &prims);

          ghl_return_primitives(
                &prims,
                &rho_b[index], &pressure[index], &eps[index],
                &vx[index], &vy[index], &vz[index],
                &Bx_center[index], &By_center[index], &Bz_center[index],
                &dummy1, &dummy2, &dummy3);
        }
      }
    }

    /* ZMIN & ZMAX */
    // k=kmax=outer boundary
    if(cctk_bbox[5]) {
      const int kmax=cctk_lsh[2]-cctk_nghostzones[2]+which_bdry_pt;
#pragma omp parallel for
      for(int j=0; j<cctk_lsh[1]; j++) {
        for(int i=0; i<cctk_lsh[0]; i++) {
          const int index = CCTK_GFINDEX3D(cctkGH, i, j, kmax);
          const int indm1 = CCTK_GFINDEX3D(cctkGH, i, j, kmax-1);

          const double vtmp = vz[indm1] < 0.0 ? 0 : vz[indm1];
          ghl_primitive_quantities prims;
          ghl_initialize_primitives(
                rho_b[indm1], pressure[indm1], eps[indm1],
                vx[indm1], vy[indm1], vtmp,
                Bx_center[index], By_center[index], Bz_center[index],
                poison, poison, poison, &prims);

          GRHayLMHD_enforce_primitive_limits_and_compute_conservs(cctkGH, index, &prims);

          ghl_return_primitives(
                &prims,
                &rho_b[index], &pressure[index], &eps[index],
                &vx[index], &vy[index], &vz[index],
                &Bx_center[index], &By_center[index], &Bz_center[index],
                &dummy1, &dummy2, &dummy3);
        }
      }
    }
    // k=kmin=outer boundary
    if((cctk_bbox[4]) && Symmetry_none) {
      const int kmin=cctk_nghostzones[2]-which_bdry_pt-1;
#pragma omp parallel for
      for(int j=0; j<cctk_lsh[1]; j++) {
        for(int i=0; i<cctk_lsh[0]; i++) {
          const int index = CCTK_GFINDEX3D(cctkGH, i, j, kmin);
          const int indp1 = CCTK_GFINDEX3D(cctkGH, i, j, kmin+1);

          const double vtmp = vz[indp1] > 0.0 ? 0 : vz[indp1];
          ghl_primitive_quantities prims;
          ghl_initialize_primitives(
                rho_b[indp1], pressure[indp1], eps[indp1],
                vx[indp1], vy[indp1], vtmp,
                Bx_center[index], By_center[index], Bz_center[index],
                poison, poison, poison, &prims);

          GRHayLMHD_enforce_primitive_limits_and_compute_conservs(cctkGH, index, &prims);

          ghl_return_primitives(
                &prims,
                &rho_b[index], &pressure[index], &eps[index],
                &vx[index], &vy[index], &vz[index],
                &Bx_center[index], &By_center[index], &Bz_center[index],
                &dummy1, &dummy2, &dummy3);
        }
      }
    }
  }
}

void GRHayLMHD_enforce_primitive_limits_and_compute_conservs(const cGH* cctkGH, const int index, ghl_primitive_quantities *restrict prims) {
  // We cheat here by using the argument list of the scheduled function
  // instead of explicitly passing all these variables.
  DECLARE_CCTK_ARGUMENTS_GRHayLMHD_outer_boundaries_on_P_rho_b_vx_vy_vz;

  double dummy1, dummy2;

  ghl_metric_quantities ADM_metric;
  ghl_enforce_detgtij_and_initialize_ADM_metric(
        alp[index],
        betax[index], betay[index], betaz[index],
        gxx[index], gxy[index], gxz[index],
        gyy[index], gyz[index], gzz[index],
        &ADM_metric);

  ghl_ADM_aux_quantities metric_aux;
  ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

  ghl_conservative_quantities cons;
  const int speed_limited CCTK_ATTRIBUTE_UNUSED = ghl_enforce_primitive_limits_and_compute_u0(
        ghl_params, ghl_eos, &ADM_metric, prims);

  ghl_compute_conservs(
        &ADM_metric, &metric_aux, prims, &cons);

  ghl_return_conservatives(
        &cons,
        &rho_star[index], &tau[index],
        &Stildex[index], &Stildey[index], &Stildez[index],
        &dummy1, &dummy2);
}
