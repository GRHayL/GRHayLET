#include "GRHayLM1.h"

void GRHayLM1_calc_closure(CCTK_ARGUMENTS){
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

#pragma omp parallel
  for (int k = 0; k<cctk_lsh[2]; k++){
    for (int j = 0; j<cctk_lsh[1]; j++){
      for (int i = 0; i<cctk_lsh[0]; i++){
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);
      
    
  

        m1_root_params root_params = {0};

        ghl_metric_quantities metric;
        ghl_initialize_metric(alp[index], betax[index], betay[index], betaz[index],
                                          gxx[index],   gxy[index],   gxz[index],
                                                        gyy[index],   gyz[index],
                                                                      gzz[index],
                                                                      &metric); //TODO(DRB): Confirm correct settings from ADMbase.
  
        ghl_ADM_aux_quantities adm_aux;
        ghl_compute_ADM_auxiliaries(&metric, &adm_aux);

        ghl_primitive_quantities prims = {0};
        prims.rho         = rho[index];
        prims.press       = press[index];
        prims.eps         = eps[index];
        prims.entropy     = entropy[index];
        prims.Y_e         = Y_e[index];
        prims.temperature = temperature[index];
        prims.vU[0]       = alp[index] * vel[CCTK_GFINDEX4D(cctkGH, i, j, k, 0)] - betax[index];
        prims.vU[1]       = alp[index] * vel[CCTK_GFINDEX4D(cctkGH, i, j, k, 1)] - betay[index];
        prims.vU[2]       = alp[index] * vel[CCTK_GFINDEX4D(cctkGH, i, j, k, 2)] - betaz[index];
        prims.BU[0]       = Bvec[CCTK_GFINDEX4D(cctkGH, i, j, k, 0)];
        prims.BU[1]       = Bvec[CCTK_GFINDEX4D(cctkGH, i, j, k, 1)];
        prims.BU[2]       = Bvec[CCTK_GFINDEX4D(cctkGH, i, j, k, 2)];//TODO(DRB): Confirm correct calcs and settings from Hydrobase.

        bool speed_limited = false;
        ghl_enforce_primitive_limits_and_compute_u0(ghl_params, ghl_eos, &metric, &prims, &speed_limited);

        ghl_radiation_flux_vector F4;
        ghl_radiation_pressure_tensor P4;

        F4.D[0] = 0.0;
        F4.D[1] = GRHayLM1_rFx[index];
        F4.D[2] = GRHayLM1_rFy[index];
        F4.D[3] = GRHayLM1_rFz[index];


        P4.DD[0][0] = 0.0;
        for (int a = 1; a < 4; a++){
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

        root_params.metric  = &metric;
        root_params.adm_aux = &adm_aux;
        root_params.prims   = &prims;
        root_params.E       = GRHayLM1_rE[index];
        root_params.F4      = &F4;
        root_params.P4      = &P4;

        ghl_radiation_rootSolve_closure(&root_params);
        ghl_radiation_apply_closure(&root_params.metric, &root_params.adm_aux, &root_params.prims,
                                    root_params.E, &root_params.F4, root_params.chi,
                                    &root_params.P4);

        double n4D[4] = {alp[index],0,0,0};
        double u4U[4] = {root_params.prims->u0,
                         root_params.prims->u0*root_params.prims->vU[0],
                         root_params.prims->u0*root_params.prims->vU[1],
                         root_params.prims->u0*root_params.prims->vU[2]};
        double v4U[4] = {0,root_params.prims->vU[0],root_params.prims->vU[1],root_params.prims->vU[2]};
        ghl_radiation_metric_tensor proj = {0}; 
        for (int a = 0; a < 4; a++){
          for (int b = 0; b < 4; b++){
            if(a==b){
              proj.UD[a][b] += 1;
            }
            for (int c = 0; c < 4; c++){
              proj.UD[a][b] += u4U[a]*u4U[c]*root_params.adm_aux->g4DD[b][c];
            }
          }
        }

        ghl_stress_energy rT4DD = {0};
        ghl_radiation_flux_vector H4 = {0};
        ghl_radiation_con_flux_vector fnu = {0};

        assemble_rT_lab_frame(&n4D, root_params.E, &root_params.F4, &root_params.P4, &rT4DD);

        GRHayLM1_rJ[index] = calc_J_from_rT(&u4U, &proj, &rT4DD);
        calc_H4D_from_rT(&u4U, &proj, &rT4DD, &H4);
        GRHayLM1_rHt[index] = H4.U[0];
        GRHayLM1_rHx[index] = H4.U[1];
        GRHayLM1_rHy[index] = H4.U[2];
        GRHayLM1_rHz[index] = H4.U[3];

        assemble_fnu(&root_params.adm_aux, &u4U, GRHayLM1_rJ[index], &H4, &fnu);
        const double Gamma = compute_Gamma(w_lorentz[index], &v4U, GRHayLM1_rJ[index], root_params.E, 0.0, 0.0, &root_params.F4); //TODO(DRB): Add in floor limits
        GRHayLM1_rnnu[index] = GRHayLM1_rN[index]/Gamma;


      }
    }
  }

}
