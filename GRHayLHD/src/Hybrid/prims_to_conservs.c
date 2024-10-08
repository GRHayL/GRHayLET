#include "GRHayLHD.h"

void GRHayLHD_hybrid_prims_to_conservs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLHD_hybrid_prims_to_conservs;
  DECLARE_CCTK_PARAMETERS;

  const int imax = cctk_lsh[0];
  const int jmax = cctk_lsh[1];
  const int kmax = cctk_lsh[2];

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
        prims.BU[0] = prims.BU[1] = prims.BU[2] = 0.0;
        prims.rho   = rho[index];
        prims.press = press[index];
        prims.vU[0] = vx[index];
        prims.vU[1] = vy[index];
        prims.vU[2] = vz[index];

        bool speed_limited; 
        const ghl_error_codes_t error = ghl_enforce_primitive_limits_and_compute_u0(
              ghl_params, ghl_eos, &ADM_metric, &prims, &speed_limited);
        if(error)
          ghl_read_error_codes(error);

        ghl_conservative_quantities cons;
        ghl_compute_conservs(
              &ADM_metric, &metric_aux, &prims, &cons);

        rho[index]   = prims.rho;
        press[index] = prims.press;
        eps[index]   = prims.eps;
        u0[index]    = prims.u0;
        vx[index]    = prims.vU[0];
        vy[index]    = prims.vU[1];
        vz[index]    = prims.vU[2];

        rho_star[index] = cons.rho;
        tau[index]      = cons.tau;
        Stildex[index]  = cons.SD[0];
        Stildey[index]  = cons.SD[1];
        Stildez[index]  = cons.SD[2];
      }
    }
  }
}
