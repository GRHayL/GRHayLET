/********************************
 * CONVERT ET ID TO IllinoisGRMHD
 *
 * Written in 2014 by Zachariah B. Etienne
 *
 * Sets metric & MHD variables needed
 * by IllinoisGRMHD, converting from
 * HydroBase and ADMBase.
 ********************************/

#include "IGM.h"

void IllinoisGRMHD_prims_to_conservs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_IllinoisGRMHD_prims_to_conservs;
  DECLARE_CCTK_PARAMETERS;

  const double poison = 0.0/0.0;
  double dummy1, dummy2, dummy3;

  const int imax = cctk_lsh[0];
  const int jmax = cctk_lsh[1];
  const int kmax = cctk_lsh[2];

#pragma omp parallel for
  for(int k=0; k<kmax; k++) {
    for(int j=0; j<jmax; j++) {
      for(int i=0; i<imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        metric_quantities ADM_metric;
        ghl_initialize_metric(
              alp[index],
              betax[index], betay[index], betaz[index],
              gxx[index], gxy[index], gxz[index],
              gyy[index], gyz[index], gzz[index],
              &ADM_metric);

        ADM_aux_quantities metric_aux;
        ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

        primitive_quantities prims;
        ghl_initialize_primitives(
              rho_b[index], pressure[index], eps[index],
              vx[index], vy[index], vz[index],
              Bx[index], By[index], Bz[index],
              poison, poison, poison,
              &prims);

        conservative_quantities cons;
        int speed_limited = 0;
        //This applies inequality fixes on the primitives
        ghl_enforce_primitive_limits_and_compute_u0(
              ghl_params, ghl_eos, &ADM_metric,
              &prims, &speed_limited);
        //This computes the conservatives and stress-energy tensor from the new primitives
        ghl_compute_conservs(
              &ADM_metric, &metric_aux, &prims, &cons);

        ghl_return_primitives(
              &prims,
              &rho_b[index], &pressure[index], &eps[index],
              &vx[index], &vy[index], &vz[index],
              &Bx[index], &By[index], &Bz[index],
              &dummy1, &dummy2, &dummy3);

        ghl_return_conservatives(
              &cons,
              &rho_star[index], &tau[index],
              &mhd_st_x[index], &mhd_st_y[index], &mhd_st_z[index],
              &dummy1, &dummy2);
      }
    }
  }
}
