/********************************
 * CONVERT ET ID TO IllinoisGRMHD
 *
 * Written in 2014 by Zachariah B. Etienne
 *
 * Sets metric & MHD variables needed
 * by IllinoisGRMHD, converting from
 * HydroBase and ADMBase.
 ********************************/

#include "GRHayLHD.h"

void GRHayLHD_prims_to_conservs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLHD_prims_to_conservs;
  DECLARE_CCTK_PARAMETERS;

  double dummy1, dummy2, dummy3;

  const int imax = cctk_lsh[0];
  const int jmax = cctk_lsh[1];
  const int kmax = cctk_lsh[2];

  // Finally, enforce limits on primitives & compute conservative variables.
#pragma omp parallel for
  for(int k=0; k<kmax; k++) {
    for(int j=0; j<jmax; j++) {
      for(int i=0; i<imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        ghl_metric_quantities ADM_metric;
        ghl_initialize_metric(
              alp[index],
              betax[index], betay[index], betaz[index],
              gxx[index], gxy[index], gxz[index],
              gyy[index], gyz[index], gzz[index],
              &ADM_metric);

        ghl_ADM_aux_quantities metric_aux;
        ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

        ghl_primitive_quantities prims;
        ghl_initialize_primitives(
              rho_b[index], pressure[index], eps[index],
              vx[index], vy[index], vz[index],
              0.0, 0.0, 0.0,
              ent[index], Ye[index], temp[index],
              &prims);

        ghl_conservative_quantities cons;
        //This applies inequality fixes on the conservatives
        const int speed_limited CCTK_ATTRIBUTE_UNUSED = ghl_enforce_primitive_limits_and_compute_u0(
              ghl_params, ghl_eos, &ADM_metric, &prims);
        //This computes the conservatives from the new primitives
        ghl_compute_conservs(
              &ADM_metric, &metric_aux, &prims, &cons);

        ghl_return_primitives(
              &prims,
              &rho_b[index], &pressure[index], &eps[index],
              &vx[index], &vy[index], &vz[index],
              &dummy1, &dummy2, &dummy3,
              &ent[index], &Ye[index], &temp[index]);

        ghl_return_conservatives(
              &cons,
              &rho_star[index], &tau[index],
              &Stildex[index], &Stildey[index], &Stildez[index],
              &ent_star[index], &Ye_star[index]);
      }
    }
  }
}
