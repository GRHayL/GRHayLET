#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "NRPyLeakageET.h"

#define CHECK_POINTER(pointer,name) \
  if( !pointer ) CCTK_VERROR("Failed to get pointer for gridfunction '%s'",name);

void NRPyLeakageET_compute_neutrino_opacities_and_add_source_terms_to_MHD_rhss(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS_NRPyLeakageET_compute_neutrino_opacities_and_add_source_terms_to_MHD_rhss;
  DECLARE_CCTK_PARAMETERS;

  if(verbosity_level>1) CCTK_VINFO("Inside NRPyLeakageET_compute_opacities_and_add_source_terms_to_MHD_rhss");

  const int timelevel = 0;

  const char *rhs_names[5] = {GFstring_Y_e_star_rhs,GFstring_tau_rhs,
                             GFstring_Stildex_rhs,GFstring_Stildey_rhs,GFstring_Stildez_rhs};
  CCTK_INT rhs_vars[5];
  const CCTK_INT rhs_tls[5] = {0,0,0,0,0};
  const CCTK_INT rhs_where[5] = {CCTK_VALID_INTERIOR,CCTK_VALID_INTERIOR,CCTK_VALID_INTERIOR,
                                CCTK_VALID_INTERIOR,CCTK_VALID_INTERIOR};
  for(int n=0;n<5;n++) {
    rhs_vars[n] = CCTK_VarIndex(rhs_names[n]);
    if(rhs_vars[n] < 0) CCTK_VERROR("Unknown RHS gridfunction '%s'",rhs_names[n]);
  }
  if(Driver_RequireValidData(cctkGH,rhs_vars,rhs_tls,5,rhs_where) != 0)
    CCTK_ERROR("Could not require valid hydro RHS data");

  // Step 1: Get pointers to opacity and optical depth gridfunctions
  CCTK_REAL *Y_e_star_rhs = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH, timelevel, GFstring_Y_e_star_rhs));
  CCTK_REAL *tau_rhs      = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH, timelevel, GFstring_tau_rhs));
  CCTK_REAL *Stildex_rhs  = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH, timelevel, GFstring_Stildex_rhs));
  CCTK_REAL *Stildey_rhs  = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH, timelevel, GFstring_Stildey_rhs));
  CCTK_REAL *Stildez_rhs  = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH, timelevel, GFstring_Stildez_rhs));

  // Step 2: Check pointers are ok
  CHECK_POINTER(Y_e_star_rhs, GFstring_Y_e_star_rhs);
  CHECK_POINTER(tau_rhs     , GFstring_tau_rhs     );
  CHECK_POINTER(Stildex_rhs , GFstring_Stildex_rhs );
  CHECK_POINTER(Stildey_rhs , GFstring_Stildey_rhs );
  CHECK_POINTER(Stildez_rhs , GFstring_Stildez_rhs );

  // Step 3: Ghostzones begin and end index
  const int imin = cctk_nghostzones[0];
  const int imax = cctk_lsh[0] - cctk_nghostzones[0];
  const int jmin = cctk_nghostzones[1];
  const int jmax = cctk_lsh[1] - cctk_nghostzones[1];
  const int kmin = cctk_nghostzones[2];
  const int kmax = cctk_lsh[2] - cctk_nghostzones[2];

  // Step 4: Compute opacities and leakage source terms
  int num_points=0;
  CCTK_REAL Y_e_star_rhs_avg=0,tau_rhs_avg=0,Stildex_rhs_avg=0,Stildey_rhs_avg=0,Stildez_rhs_avg=0;
#pragma omp parallel for reduction(+:num_points,Y_e_star_rhs_avg,tau_rhs_avg,Stildex_rhs_avg,Stildey_rhs_avg,Stildez_rhs_avg)
  for(int k=kmin;k<kmax;k++) {
    for(int j=jmin;j<jmax;j++) {
      for(int i=imin;i<imax;i++) {

        // Step 4.a: Set the index of the current gridpoint
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        // Step 4.b: Check if we are within the threshold
        const CCTK_REAL rhoL = rho[index];
        if( rhoL < rho_min_threshold || rhoL > rho_max_threshold ) {
          // Step 4.b.i: Below density threshold.
          //             Set opacities to zero; don't add anything to the RHSs
          kappa_0_nue [index] = 0.0;
          kappa_1_nue [index] = 0.0;
          kappa_0_anue[index] = 0.0;
          kappa_1_anue[index] = 0.0;
          kappa_0_nux [index] = 0.0;
          kappa_1_nux [index] = 0.0;
        }
        else {
          CCTK_REAL gxxL        = gxx[index];
          CCTK_REAL gxyL        = gxy[index];
          CCTK_REAL gxzL        = gxz[index];
          CCTK_REAL gyyL        = gyy[index];
          CCTK_REAL gyzL        = gyz[index];
          CCTK_REAL gzzL        = gzz[index];
          const CCTK_REAL gdet  = (gxxL * gyyL * gzzL
                                     + gxyL * gyzL * gxzL
                                     + gxzL * gxyL * gyzL
                                     - gxzL * gyyL * gxzL
                                     - gxyL * gxyL * gzzL
                                     - gxxL * gyzL * gyzL);
          if(!robust_isfinite(gxxL) || !robust_isfinite(gxyL) || !robust_isfinite(gxzL) ||
             !robust_isfinite(gyyL) || !robust_isfinite(gyzL) || !robust_isfinite(gzzL) ||
             !robust_isfinite(gdet) || gxxL <= 0 || gxxL*gyyL-gxyL*gxyL <= 0 || gdet <= 0) {
            CCTK_VERROR("Invalid spatial metric at (%d,%d,%d), level %d",i,j,k,GetRefinementLevel(cctkGH));
            continue;
          }
          const CCTK_REAL phiL  = (1.0/12.0) * log(gdet);
          const CCTK_REAL psiL  = exp(phiL);
          const CCTK_REAL psi2L = psiL *psiL;
          const CCTK_REAL psi4L = psi2L*psi2L;
          const CCTK_REAL psi6L = psi4L*psi2L;
          if( psi6L > psi6_threshold ) {
            kappa_0_nue [index] = 0.0;
            kappa_1_nue [index] = 0.0;
            kappa_0_anue[index] = 0.0;
            kappa_1_anue[index] = 0.0;
            kappa_0_nux [index] = 0.0;
            kappa_1_nux [index] = 0.0;
          }
          else {
            // Step 4.c: Read from main memory
            const CCTK_REAL alpL         = alp[index];
            const CCTK_REAL betaxL       = betax[index];
            const CCTK_REAL betayL       = betay[index];
            const CCTK_REAL betazL       = betaz[index];
            const CCTK_REAL Y_eL         = Y_e[index];
            const CCTK_REAL temperatureL = temperature[index];
            ghl_neutrino_optical_depths tauL;
            tauL.nue [0] = tau_0_nue_p [index];
            tauL.nue [1] = tau_1_nue_p [index];
            tauL.anue[0] = tau_0_anue_p[index];
            tauL.anue[1] = tau_1_anue_p[index];
            tauL.nux [0] = tau_0_nux_p [index];
            tauL.nux [1] = tau_1_nux_p [index];

            // HydroBase supplies Eulerian velocity V^i. Validate before sqrt,
            // then scale its metric norm to 1-W_max^-2 if the leakage cap is reached.
            CCTK_REAL Vx = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)];
            CCTK_REAL Vy = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)];
            CCTK_REAL Vz = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)];
            const CCTK_REAL speed2 = gxxL*Vx*Vx + gyyL*Vy*Vy + gzzL*Vz*Vz
                                  + 2*(gxyL*Vx*Vy + gxzL*Vx*Vz + gyzL*Vy*Vz);
            if(!robust_isfinite(alpL) || alpL <= 0 || !robust_isfinite(betaxL) ||
               !robust_isfinite(betayL) || !robust_isfinite(betazL) ||
               !robust_isfinite(Vx) || !robust_isfinite(Vy) || !robust_isfinite(Vz) ||
               !robust_isfinite(speed2) || speed2 < 0 || speed2 >= 1 ||
               !robust_isfinite(W_max) || W_max < 1) {
              CCTK_VERROR("Invalid leakage velocity/geometry at (%d,%d,%d), level %d: alpha=%g speed2=%g W_max=%g",
                          i,j,k,GetRefinementLevel(cctkGH),alpL,speed2,W_max);
              continue;
            }
            CCTK_REAL W = 1/sqrt(1-speed2);
            if(W > W_max) {
              const CCTK_REAL scale = sqrt((1-1/(W_max*W_max))/speed2);
              Vx *= scale; Vy *= scale; Vz *= scale;
              W = W_max;
            }
            const CCTK_REAL vxL = alpL*Vx - betaxL;
            const CCTK_REAL vyL = alpL*Vy - betayL;
            const CCTK_REAL vzL = alpL*Vz - betazL;

            // Step 4.f: Compute u^{mu} using:
            //  - W = alpha u^{0}     => u^{0} = W / alpha
            //  - v^{i} = u^{i}/u^{0} => u^{i} = v^{i}u^{0}
            const CCTK_REAL u0L = W / alpL;
            const CCTK_REAL uxL = vxL * u0L;
            const CCTK_REAL uyL = vyL * u0L;
            const CCTK_REAL uzL = vzL * u0L;

            // Step 4.g: Compute u_{mu} = g_{mu nu}u^{nu}
            // Step 4.g.i: Set g_{mu nu}
            // Step 4.g.i.1: Set gamma_{ij}
            CCTK_REAL gammaDD[3][3];
            gammaDD[0][0] = gxxL;
            gammaDD[0][1] = gammaDD[1][0] = gxyL;
            gammaDD[0][2] = gammaDD[2][0] = gxzL;
            gammaDD[1][1] = gyyL;
            gammaDD[1][2] = gammaDD[2][1] = gyzL;
            gammaDD[2][2] = gzzL;

            // Step 4.g.i.2: Compute beta_{i}
            CCTK_REAL betaU[3] = {betaxL,betayL,betazL};
            CCTK_REAL betaD[3] = {0.0,0.0,0.0};
            for(int ii=0;ii<3;ii++)
              for(int jj=0;jj<3;jj++)
                betaD[ii] += gammaDD[ii][jj] * betaU[jj];

            // Step 4.g.i.3: Compute beta^{2} = beta_{i}beta^{i}
            CCTK_REAL betasqr = 0.0;
            for(int ii=0;ii<3;ii++) betasqr += betaD[ii]*betaU[ii];

            // Step 4.g.i.4: Set g_{mu nu}
            CCTK_REAL g4DD[4][4];
            g4DD[0][0] = -alpL*alpL + betasqr;
            for(int ii=0;ii<3;ii++) {
              g4DD[0][ii+1] = g4DD[ii+1][0] = betaD[ii];
              for(int jj=ii;jj<3;jj++) {
                g4DD[ii+1][jj+1] = g4DD[jj+1][ii+1] = gammaDD[ii][jj];
              }
            }

            // Step 4.g.ii: Compute u_{mu} = g_{mu nu}u^{nu}
            CCTK_REAL u4U[4] = {u0L,uxL,uyL,uzL};
            CCTK_REAL u4D[4] = {0.0,0.0,0.0,0.0};
            for(int mu=0;mu<4;mu++)
              for(int nu=0;nu<4;nu++)
                u4D[mu] += g4DD[mu][nu] * u4U[nu];

            // Step 4.h: Compute R, Q, and the neutrino opacities
            ghl_neutrino_opacities kappaL;
            CCTK_REAL R_sourceL, Q_sourceL;
            const ghl_error_codes_t status = NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms(ghl_eos,
                                                                          rhoL, Y_eL, temperatureL,
                                                                          &tauL, &kappaL, &R_sourceL, &Q_sourceL);

            if(status != ghl_success) {
              CCTK_VERROR("Leakage source evaluation failed (status %d) at (%d,%d,%d), level %d: rho=%g Ye=%g T=%g",
                          (int)status,i,j,k,GetRefinementLevel(cctkGH),rhoL,Y_eL,temperatureL);
              continue;
            }

            // Step 4.i: Compute MHD right-hand sides
            const CCTK_REAL sqrtmgL       = alpL * psi6L;
            const CCTK_REAL sqrtmgR       = sqrtmgL * R_sourceL;
            const CCTK_REAL sqrtmgQ       = sqrtmgL * Q_sourceL;
            const CCTK_REAL Y_e_star_rhsL = sqrtmgR;
            const CCTK_REAL tau_rhsL      = alpL * sqrtmgQ * u4U[0];
            const CCTK_REAL Stildex_rhsL  = sqrtmgQ * u4D[1];
            const CCTK_REAL Stildey_rhsL  = sqrtmgQ * u4D[2];
            const CCTK_REAL Stildez_rhsL  = sqrtmgQ * u4D[3];

            // Check each result and proposed sum before committing any cell output.
            const CCTK_REAL increments[5] = {Y_e_star_rhsL,tau_rhsL,Stildex_rhsL,Stildey_rhsL,Stildez_rhsL};
            CCTK_REAL *rhss[5] = {Y_e_star_rhs,tau_rhs,Stildex_rhs,Stildey_rhs,Stildez_rhs};
            int valid = NRPyLeakageET_opacities_finite(&kappaL) && robust_isfinite(R_sourceL) && robust_isfinite(Q_sourceL);
            for(int mu=0;mu<4;mu++)
              valid = valid && robust_isfinite(u4U[mu]) && robust_isfinite(u4D[mu]);
            for(int n=0;n<5;n++)
              valid = valid && robust_isfinite(increments[n]) && robust_isfinite(rhss[n][index])
                            && robust_isfinite(rhss[n][index]+increments[n]);
            if(!valid) {
              CCTK_VERROR("Nonfinite leakage output/RHS at (%d,%d,%d), level %d: rho=%g Ye=%g T=%g W=%g R=%g Q=%g",
                          i,j,k,GetRefinementLevel(cctkGH),rhoL,Y_eL,temperatureL,W,R_sourceL,Q_sourceL);
              continue;
            }

            // Step 4.j: Write to main memory
            kappa_0_nue [index]  = kappaL.nue [0];
            kappa_1_nue [index]  = kappaL.nue [1];
            kappa_0_anue[index]  = kappaL.anue[0];
            kappa_1_anue[index]  = kappaL.anue[1];
            kappa_0_nux [index]  = kappaL.nux [0];
            kappa_1_nux [index]  = kappaL.nux [1];

            // Step 4.k: Update right-hand sides only in the grid interior
            Y_e_star_rhs[index] += Y_e_star_rhsL;
            tau_rhs     [index] += tau_rhsL;
            Stildex_rhs [index] += Stildex_rhsL;
            Stildey_rhs [index] += Stildey_rhsL;
            Stildez_rhs [index] += Stildez_rhsL;

            Y_e_star_rhs_avg += Y_e_star_rhsL;
            tau_rhs_avg      += tau_rhsL;
            Stildex_rhs_avg  += Stildex_rhsL;
            Stildey_rhs_avg  += Stildey_rhsL;
            Stildez_rhs_avg  += Stildez_rhsL;

            num_points++;
          }
        }
      }
    }
  }
  CCTK_REAL inv_numpts = num_points > 0 ? 1.0/((CCTK_REAL)num_points) : 1.0;
  if(verbosity_level>0) {
    CCTK_VINFO("***** Iter. # %d, Lev: %d, Averages -- Ye_rhs: %e | tau_rhs: %e | st_i_rhs: %e,%e,%e *****",cctk_iteration,GetRefinementLevel(cctkGH),
               Y_e_star_rhs_avg*inv_numpts,
               tau_rhs_avg*inv_numpts,
               Stildex_rhs_avg*inv_numpts,
               Stildey_rhs_avg*inv_numpts,
               Stildez_rhs_avg*inv_numpts);
    if(verbosity_level>1) CCTK_INFO("Finished NRPyLeakageET_compute_opacities_and_add_source_terms_to_MHD_rhss");
  }
  if(Driver_NotifyDataModified(cctkGH,rhs_vars,rhs_tls,5,rhs_where) != 0)
    CCTK_ERROR("Could not notify hydro RHS modification");
}
