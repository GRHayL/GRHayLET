#include "GRHayLM1.h"

void GRHayLM1_compute_Tmunu(CCTK_ARGUMENTS){
  DECLARE_CCTK_ARGUMENTS_GRHayLM1_compute_Tmunu;
  DECLARE_CCTK_PARAMETERS;

#pragma omp parallel //FIXME: Couldn't I just use CCTK_LOOP3_ALL?
  for(int k=0; k<cctk_lsh[2]; k++){ 
    for(int j=0; j<cctk_lsh[1]; j++){
      for(int i=0; i<cctk_lsh[0]; i++){
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);

        //Read in ADM quantities from grid functions
        //set ghl quantities
        //TODO: Will m1_root_params be accessible globally like other ghl quantities? That means I would have access to E/F/P for assembly.
        //I don't think I should have to be calcing the closure within this function...I think that should be done elsewhere, and then I access the E/F/P.
        m1_root_params ghl_m1_root_params;
        ghl_enforce_detgtij_and_initialize_ADM_metric(
          alp[index],
          betax[index], betay[index], betaz[index],
          gxx[index], gxy[index], gxz[index],
          gyy[index], gyz[index], gzz[index],
          &ghl_m1_root_params.metric);

        ghl_compute_ADM_auxiliaries(&ghl_m1_root_params.metric, &ghl_m1_root_params.adm_aux);

        //Read in quantities for M1 treatment don't think I need these...
        ghl_m1_root_params.prims.BU[0] = ghl_m1_root_params.prims.BU[1] = ghl_m1_root_params.prims.BU[2] = 0.0;
        ghl_m1_root_params.prims.rho = rho[index];
        ghl_m1_root_params.prims.press = press[index];
        ghl_m1_root_params.prims.eps = eps[index];
        ghl_m1_root_params.prims.vU[0] = vx[index];
        ghl_m1_root_params.prims.vU[1] = vy[index];
        ghl_m1_root_params.prims.vU[2] = vz[index];
        ghl_m1_root_params.prims.u0 = u0[index];


        //TODO: Fill in to assemble Tmunu
        ghl_stress_energy rTmunu;
        //double chi = 0.5;//chi determined by choice


        //ghl_radiation_rootSolve_closure(&fparams);
        //ghl_radiation_apply_closure(&fparams.metric, &fparams.adm_aux, &fparams.prims,
        //                            fparams.E, &fparams.F4, chi, &fparams.P4);
        assemble_rT_lab_frame(&n4D, ghl_m1_root_params.E, &ghl_m1_root_params.F4, &ghl_m1_root_params.P4, &rTmunu);

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
