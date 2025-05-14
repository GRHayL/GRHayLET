#include "GRHayLM1.h"

void GRHayLM1_compute_Tmunu(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

#pragma omp parallel // FIXME: Couldn't I just use CCTK_LOOP3_ALL?
  for (int k = 0; k < cctk_lsh[2]; k++) {
    for (int j = 0; j < cctk_lsh[1]; j++) {
      for (int i = 0; i < cctk_lsh[0]; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);

        // Read in ADM quantities from grid functions
        // set ghl quantities
        ghl_metric_quantities adm_metric = {0};
        ghl_enforce_detgtij_and_initialize_ADM_metric(
            alp[index], betax[index], betay[index], betaz[index], gxx[index],
            gxy[index], gxz[index], gyy[index], gyz[index], gzz[index],
            &adm_metric);

        ghl_ADM_aux_quantities aux_metric = {0};
        ghl_compute_ADM_auxiliaries(&adm_metric, &aux_metric);

        // Read in quantities for M1 treatment
        ghl_primitive_quantities prims = {0};
        const CCTK_REAL lapse = adm_metric.lapse;
        const CCTK_REAL shiftx = adm_metric.betaU[0];
        const CCTK_REAL shifty = adm_metric.betaU[1];
        const CCTK_REAL shiftz = adm_metric.betaU[2];
        prims.rho = rho[index];
        prims.press = press[index];
        prims.eps = eps[index];
        prims.vU[0] = lapse * vel[CCTK_GFINDEX4D(cctkGH, i, j, k, 0)] - shiftx;
        prims.vU[1] = lapse * vel[CCTK_GFINDEX4D(cctkGH, i, j, k, 1)] - shifty;
        prims.vU[2] = lapse * vel[CCTK_GFINDEX4D(cctkGH, i, j, k, 2)] - shiftz;

        bool speed_limited = false;
        ghl_limit_v_and_compute_u0(ghl_params, &adm_metric, &prims,
                                   &speed_limited);

        // TODO: Fill in to assemble Tmunu
        ghl_stress_energy rTmunu;
        ghl_radiation_flux_vector F4;
        ghl_radiation_pressure_tensor P4;

        F4.D[0] = 0.0;
        F4.D[1] = GRHayLM1_rFx[index];
        F4.D[2] = GRHayLM1_rFy[index];
        F4.D[3] = GRHayLM1_rFz[index];

        P4.DD[0][0] = 0.0;
        for (int a = 0; a < 4; a++){
          P4.DD[0][a] = 0.0;
          P4.DD[a][0] = 0.0;
        }
        P4.DD[1][1] = GRHayLM1_rPxx[index];
        P4.DD[1][2] = GRHayLM1_rPxy[index];
        P4.DD[1][3] = GRHayLM1_rPxz[index];
        P4.DD[2][2] = GRHayLM1_rPyy[index];
        P4.DD[2][3] = GRHayLM1_rPyz[index];
        P4.DD[3][3] = GRHayLM1_rPzz[index];
        P4.DD[2][1] = P4.DD[1][2];
        P4.DD[3][1] = P4.DD[1][3];
        P4.DD[3][2] = P4.DD[2][3];

        m1_root_params root_params = {0};
        root_params.prims = &prims;
        root_params.metric = &adm_metric;
        root_params.adm_aux = &aux_metric;
        root_params.E = GRHayLM1_rE[index];
        root_params.F4 = &F4;
        root_params.P4 = &P4;
        
        ghl_radiation_rootSolve_closure(&root_params);
        ghl_radiation_apply_closure(&root_params.metric, &root_params.adm_aux, &root_params.prims,
                                     root_params.E, &root_params.F4, root_params.chi,
                                     &root_params.P4);
        
        double n4D[4] = {-lapse, 0, 0, 0};
        assemble_rT_lab_frame(&n4D, root_params.E, &root_params.F4,
                              &root_params.P4, &rTmunu);

        eTtt[index] += rTmunu.T4[0][0];
        eTtx[index] += rTmunu.T4[0][1];
        eTty[index] += rTmunu.T4[0][2];
        eTtz[index] += rTmunu.T4[0][3];
        eTxx[index] += rTmunu.T4[1][1];
        eTxy[index] += rTmunu.T4[1][2];
        eTxz[index] += rTmunu.T4[1][3];
        eTyy[index] += rTmunu.T4[2][2];
        eTyz[index] += rTmunu.T4[2][3];
        eTzz[index] += rTmunu.T4[3][3];
      }
    }
  }
}
