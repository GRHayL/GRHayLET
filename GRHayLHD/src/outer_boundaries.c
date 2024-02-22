/*******************************************************
 * Outer boundaries are handled as follows:
 * (-1) Update RHS quantities, leave RHS quantities zero on all outer ghostzones (including outer AMR refinement, processor, and outer boundaries)
 * ( 0) Let MoL update all evolution variables
 * ( 1) Apply outer boundary conditions (BCs) on A_{\mu}
 * ( 2) Compute B^i from A_i everywhere, synchronize B^i
 * ( 3) Call con2prim to get primitives on interior pts
 * ( 4) Apply outer BCs on {P,rho,vx,vy,vz}.
 *******************************************************/

#include "GRHayLHD.h"


/*******************************************************
 * Apply outer boundary conditions on {P,rho,vx,vy,vz}
 * It is better to apply BCs on primitives than conservs,
 * because small errors in conservs can be greatly
 * amplified in con2prim, sometimes leading to unphysical
 * primitives & unnecessary fixes.
 *******************************************************/
void GRHayLHD_outer_boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(CCTK_EQUALS(Matter_BC, "frozen")) return;

  const bool do_outflow = CCTK_EQUALS(Matter_BC,"outflow");

  const bool Symmetry_none = CCTK_EQUALS(Symmetry,"none") ? true : false;

  // Don't apply approximate outer boundary conditions on initial data, which
  // should be defined everywhere, or on levels != [coarsest level].
  if(cctk_iteration==0 || GetRefinementLevel(cctkGH)!=0) return;

  if(cctk_nghostzones[0] != cctk_nghostzones[1] || cctk_nghostzones[0] != cctk_nghostzones[2])
    CCTK_ERROR("ERROR: GRHayLHD outer BC driver does not support unequal number of ghostzones in different directions!");
  for(int which_bdry_pt=0; which_bdry_pt<cctk_nghostzones[0]; which_bdry_pt++) {

    /* XMIN & XMAX */
    // i=imax=outer boundary
    if(cctk_bbox[1]) {
      const int imax = cctk_lsh[0] - cctk_nghostzones[0] + which_bdry_pt;

      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int j=0; j<cctk_lsh[1]; j++) {
          const int index = CCTK_GFINDEX3D(cctkGH,imax, j, k);
          const int indm1 = CCTK_GFINDEX3D(cctkGH,imax-1, j, k);
          rho[index]         = rho[indm1];
          press[index]       = press[indm1];
          vx[index]          = vx[indm1];
          vy[index]          = vy[indm1];
          vz[index]          = vz[indm1];
          eps[index]         = eps[indm1];
          Y_e[index]         = Y_e[indm1];
          temperature[index] = temperature[indm1];
          if (do_outflow && vx[index] < 0.0)
            vx[index] = 0.0;
        }
      }
    }
    // i=imin=outer boundary
    if(cctk_bbox[0]) {
      const int imin = cctk_nghostzones[0] - which_bdry_pt - 1;

      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int j=0; j<cctk_lsh[1]; j++) {
          const int index = CCTK_GFINDEX3D(cctkGH, imin, j, k);
          const int indp1 = CCTK_GFINDEX3D(cctkGH, imin+1, j, k);

          rho[index]         = rho[indp1];
          press[index]       = press[indp1];
          vx[index]          = vx[indp1];
          vy[index]          = vy[indp1];
          vz[index]          = vz[indp1];
          eps[index]         = eps[indp1];
          Y_e[index]         = Y_e[indp1];
          temperature[index] = temperature[indp1];
          if (do_outflow && vx[index] > 0.0)
            vx[index] = 0.0;
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

          rho[index]         = rho[indm1];
          press[index]       = press[indm1];
          vx[index]          = vx[indm1];
          vy[index]          = vy[indm1];
          vz[index]          = vz[indm1];
          eps[index]         = eps[indm1];
          Y_e[index]         = Y_e[indm1];
          temperature[index] = temperature[indm1];
          if (do_outflow && vy[index] < 0.0)
            vy[index] = 0.0;
        }
      }
    }
    // j=jmin=outer boundary
    if(cctk_bbox[2]) {
      const int jmin = cctk_nghostzones[1] - which_bdry_pt - 1;

      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int i=0; i<cctk_lsh[0]; i++) {
          const int index = CCTK_GFINDEX3D(cctkGH, i, jmin, k);
          const int indp1 = CCTK_GFINDEX3D(cctkGH, i, jmin+1, k);

          rho[index]         = rho[indp1];
          press[index]       = press[indp1];
          vx[index]          = vx[indp1];
          vy[index]          = vy[indp1];
          vz[index]          = vz[indp1];
          eps[index]         = eps[indp1];
          Y_e[index]         = Y_e[indp1];
          temperature[index] = temperature[indp1];
          if (do_outflow && vy[index] > 0.0)
            vy[index] = 0.;
        }
      }
    }

    /* ZMIN & ZMAX */
    // k=kmax=outer boundary
    if(cctk_bbox[5]) {
      const int kmax = cctk_lsh[2] - cctk_nghostzones[2] + which_bdry_pt;

      for(int j=0; j<cctk_lsh[1]; j++) {
        for(int i=0; i<cctk_lsh[0]; i++) {
          const int index = CCTK_GFINDEX3D(cctkGH, i, j, kmax);
          const int indm1 = CCTK_GFINDEX3D(cctkGH, i, j, kmax-1);

          rho[index]         = rho[indm1];
          press[index]       = press[indm1];
          vx[index]          = vx[indm1];
          vy[index]          = vy[indm1];
          vz[index]          = vz[indm1];
          eps[index]         = eps[indm1];
          Y_e[index]         = Y_e[indm1];
          temperature[index] = temperature[indm1];
          if (do_outflow && vz[index] < 0.0)
            vz[index] = 0.;
        }
      }
    }
    // k=kmin=outer boundary
    if((cctk_bbox[4]) && Symmetry_none) {
      const int kmin = cctk_nghostzones[2] - which_bdry_pt - 1;

      for(int j=0; j<cctk_lsh[1]; j++) {
        for(int i=0; i<cctk_lsh[0]; i++) {
          const int index = CCTK_GFINDEX3D(cctkGH, i, j, kmin);
          const int indp1 = CCTK_GFINDEX3D(cctkGH, i, j, kmin+1);

          rho[index]         = rho[indp1];
          press[index]       = press[indp1];
          vx[index]          = vx[indp1];
          vy[index]          = vy[indp1];
          vz[index]          = vz[indp1];
          eps[index]         = eps[indp1];
          Y_e[index]         = Y_e[indp1];
          temperature[index] = temperature[indp1];
          if (do_outflow && vz[index] > 0.0)
            vz[index] = 0.;
        }
      }
    }
  }

#pragma omp parallel for
  for (int k = 0; k < cctk_lsh[2]; k++) {
    for (int j = 0; j < cctk_lsh[1]; j++) {
      for (int i = 0; i < cctk_lsh[0]; i++) {
        if (((cctk_bbox[0]) && i < cctk_nghostzones[0]) ||
            ((cctk_bbox[1]) && i >= cctk_lsh[0] - cctk_nghostzones[0]) ||
            ((cctk_bbox[2]) && j < cctk_nghostzones[1]) ||
            ((cctk_bbox[3]) && j >= cctk_lsh[1] - cctk_nghostzones[1]) ||
            ((cctk_bbox[4]) && k < cctk_nghostzones[2] &&
             CCTK_EQUALS(Symmetry, "none")) ||
            ((cctk_bbox[5]) && k >= cctk_lsh[2] - cctk_nghostzones[2])) {
          int index = CCTK_GFINDEX3D(cctkGH, i, j, k);

          ghl_metric_quantities ADM_metric;
          ghl_enforce_detgtij_and_initialize_ADM_metric(
              alp[index], betax[index], betay[index], betaz[index], gxx[index],
              gxy[index], gxz[index], gyy[index], gyz[index], gzz[index],
              &ADM_metric);

          ghl_ADM_aux_quantities metric_aux;
          ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

          ghl_primitive_quantities prims;
          ghl_initialize_primitives(rho[index], press[index], eps[index],
                                    vx[index], vy[index], vz[index], 0.0, 0.0,
                                    0.0, 0.0, Y_e[index], temperature[index],
                                    &prims);

          ghl_conservative_quantities cons;
          const int speed_limited CCTK_ATTRIBUTE_UNUSED =
              ghl_enforce_primitive_limits_and_compute_u0(ghl_params, ghl_eos,
                                                          &ADM_metric, &prims);

          ghl_compute_conservs(&ADM_metric, &metric_aux, &prims, &cons);

          rho[index] = prims.rho;
          press[index] = prims.press;
          eps[index] = prims.eps;
          vx[index] = prims.vU[0];
          vy[index] = prims.vU[1];
          vz[index] = prims.vU[2];
          Y_e[index] = prims.Y_e;
          temperature[index] = prims.temperature;

          rho_star[index] = cons.rho;
          tau[index] = cons.tau;
          Stildex[index] = cons.SD[0];
          Stildey[index] = cons.SD[1];
          Stildez[index] = cons.SD[2];
          Ye_star[index] = cons.Y_e;
        }
      }
    }
  }
}
