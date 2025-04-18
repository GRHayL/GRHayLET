#include "IllinoisGRMHD.h"

void IllinoisGRMHD_evaluate_phitilde_and_A_gauge_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_IllinoisGRMHD_evaluate_phitilde_and_A_gauge_rhs;
  DECLARE_CCTK_PARAMETERS;

  // Note that in this function, we don't bother with reconstruction, instead interpolating.
  // We need A^i, but only have A_i. So we use the BSSN metric gtupij.
  // The reconstruction variables are temporary variables and the data in them can be safely overwritten,
  // saving some memory.
  CCTK_REAL *betax_interp = vxr;
  CCTK_REAL *betay_interp = vyr;
  CCTK_REAL *betaz_interp = vzr;
  CCTK_REAL *alpha_interp = vxrr;
  CCTK_REAL *alpha_Phi_minus_betaj_A_j_interp = vxll;
  CCTK_REAL *sqrtg_Ax_interp = vxl;
  CCTK_REAL *sqrtg_Ay_interp = vyl;
  CCTK_REAL *sqrtg_Az_interp = vzl;

  /* Compute \partial_t psi6phi = -\partial_i (  \alpha psi^6 A^i - psi6phi \beta^i)
   *    (Eq 13 of http://arxiv.org/pdf/1110.4633.pdf), using Lorenz gauge.
   * Note that the RHS consists of a shift advection term on psi6phi and
   *    a term depending on the vector potential.
   * psi6phi is defined at (i+1/2,j+1/2,k+1/2), but instead of reconstructing
   *    to compute the RHS of \partial_t psi6phi, we instead use standard
   *    interpolations.
   */

  // We declare these values to be over the interior so the setting of ghostzone points is
  // more transparent.
  const int imin = cctkGH->cctk_nghostzones[0];
  const int jmin = cctkGH->cctk_nghostzones[1];
  const int kmin = cctkGH->cctk_nghostzones[2];
  const int imax = cctkGH->cctk_lsh[0]-cctkGH->cctk_nghostzones[0];
  const int jmax = cctkGH->cctk_lsh[1]-cctkGH->cctk_nghostzones[1];
  const int kmax = cctkGH->cctk_lsh[2]-cctkGH->cctk_nghostzones[2];

  // The RHS loop requires 2 ghostzones, so this loop must set those values, hence this loop
  // goes into the ghostzones. This loop requires a stencil of {-1,1},{-1,1},{-1,1},
  // which uses the last of the 3 ghostzones required by the simulation.
  // Note that ALL input variables are defined at ALL gridpoints, so no
  // worries about ghostzones.
#pragma omp parallel for
  for(int k=kmin-2; k<kmax+2; k++) {
    for(int j=jmin-2; j<jmax+2; j++) {
      for(int i=imin-2; i<imax+2; i++) {
        const int index=CCTK_GFINDEX3D(cctkGH,i,j,k);

        // First compute \partial_j \alpha \sqrt{\gamma} A^j (RHS of \partial_i psi6phi)
        // FIXME: Would be much cheaper & easier to unstagger A_i, raise, then interpolate A^i.
        //        However, we keep it this way to be completely compatible with the original
        //        Illinois GRMHD thorn, called mhd_evolve.
        //
        //Step 1) j=x: Need to raise A_i, but to do that, we must have all variables at the same gridpoints:
        // The goal is to compute \partial_j (\alpha \sqrt{\gamma} A^j) at (i+1/2,j+1/2,k+1/2)
        //    We do this by first interpolating (RHS1x) = (\alpha \sqrt{\gamma} A^x) at
        //    (i,j+1/2,k+1/2)and (i+1,j+1/2,k+1/2), then taking \partial_x (RHS1x) =
        //    [ RHS1x(i+1,j+1/2,k+1/2) - RHS1x(i,j+1/2,k+1/2) ]/dX.
        // First bring gtup's, psi, and alpha to (i,j+1/2,k+1/2):
        ghl_metric_quantities metric_stencil[2][2][2];
        CCTK_REAL Ax_stencil[3][3][3];
        CCTK_REAL Ay_stencil[3][3][3];
        CCTK_REAL Az_stencil[3][3][3];
        ghl_induction_interp_vars interp_vars;

        // Read in variable at interpolation stencil points from main memory.
        for(int iterz=0; iterz<2; iterz++)
          for(int itery=0; itery<2; itery++)
            for(int iterx=0; iterx<2; iterx++) {
              const int ind = CCTK_GFINDEX3D(cctkGH,i+iterx,j+itery,k+iterz);

              ghl_initialize_metric(
                    alp[ind], betax[ind], betay[ind], betaz[ind],
                    gxx[ind], gxy[ind], gxz[ind],
                    gyy[ind], gyz[ind], gzz[ind],
                    &metric_stencil[iterz][itery][iterx]);
        }
        // A_x needs a stencil s.t. interp_limits={ 0,1,-1,1,-1,1}.
        // A_y needs a stencil s.t. interp_limits={-1,1, 0,1,-1,1}.
        // A_z needs a stencil s.t. interp_limits={-1,1,-1,1, 0,1}.
        // We could fill only the needed elements, but it is cleaner
        // to fill in the whole 3x3x3 array.
        // TODO: the old code explicitly only filled in the necessary
        // elements. If we want to remove ~15 memcopies, do that here.
        for(int iterz=-1; iterz<2; iterz++)
          for(int itery=-1; itery<2; itery++)
            for(int iterx=-1; iterx<2; iterx++) {
              const int ind = CCTK_GFINDEX3D(cctkGH,i+iterx,j+itery,k+iterz);
              Ax_stencil[iterz+1][itery+1][iterx+1] = Ax[ind];
              Ay_stencil[iterz+1][itery+1][iterx+1] = Ay[ind];
              Az_stencil[iterz+1][itery+1][iterx+1] = Az[ind];
        }
// This code should only copy the needed data that isn't copied in the loop for other variables, but it is untested.
//        for(int iter2=0; iter2<2; iter2++)
//        for(int iter1=0; iter1<2; iter1++) {
//          gauge_vars.A_x[iter2+1][0][iter1+1] = in_vars[A_XI][CCTK_GFINDEX3D(cctkGH, i+iter1,     j-1, k+iter2)]; // { (0,1),    -1, (0,1)}
//          gauge_vars.A_x[0][iter2+1][iter1+1] = in_vars[A_XI][CCTK_GFINDEX3D(cctkGH, i+iter1, j+iter2, k-1    )]; // { (0,1), (0,1), -1}
//          gauge_vars.A_y[iter2+1][iter1+1][0] = in_vars[A_YI][CCTK_GFINDEX3D(cctkGH,     i-1, j+iter1, k+iter2)]; // { -1,    (0,1), (0,1)}
//          gauge_vars.A_y[0][iter1+1][iter2+1] = in_vars[A_YI][CCTK_GFINDEX3D(cctkGH, i+iter2, j+iter1, k-1    )]; // { (0,1), (0,1), -1}
//          gauge_vars.A_z[iter1+1][iter2+1][0] = in_vars[A_ZI][CCTK_GFINDEX3D(cctkGH,     i-1, j+iter2, k+iter1)]; // { -1,    (0,1), (0,1)}
//          gauge_vars.A_z[iter1+1][0][iter2+1] = in_vars[A_ZI][CCTK_GFINDEX3D(cctkGH, i+iter2,     j-1, k+iter1)]; // { (0,1),    -1, (0,1)}
//        }

        ghl_interpolate_with_cell_centered_ADM(metric_stencil, Ax_stencil, Ay_stencil, Az_stencil, phitilde[index], &interp_vars);

        alpha_interp[index] = interp_vars.alpha;
        sqrtg_Ax_interp[index] = interp_vars.sqrtg_Ai[0];
        sqrtg_Ay_interp[index] = interp_vars.sqrtg_Ai[1];
        sqrtg_Az_interp[index] = interp_vars.sqrtg_Ai[2];
        alpha_Phi_minus_betaj_A_j_interp[index] = interp_vars.alpha_Phi_minus_betaj_A_j;
        betax_interp[index] = interp_vars.betai[0];
        betay_interp[index] = interp_vars.betai[1];
        betaz_interp[index] = interp_vars.betai[2];
      }
    }
  }

  const CCTK_REAL dxi[3] = { 1.0/CCTK_DELTA_SPACE(0), 1.0/CCTK_DELTA_SPACE(1), 1.0/CCTK_DELTA_SPACE(2) };

#pragma omp parallel for
  for(int k=kmin; k<kmax; k++) {
    for(int j=jmin; j<jmax; j++) {
      for(int i=imin; i<imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        // \partial_t A_i = [reconstructed stuff] + [gauge stuff],
        //    where [gauge stuff] = -\partial_i (\alpha \Phi - \beta^j A_j)
        Ax_rhs[index] += dxi[0]*(alpha_Phi_minus_betaj_A_j_interp[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] - alpha_Phi_minus_betaj_A_j_interp[index]);
        Ay_rhs[index] += dxi[1]*(alpha_Phi_minus_betaj_A_j_interp[CCTK_GFINDEX3D(cctkGH,i,j-1,k)] - alpha_Phi_minus_betaj_A_j_interp[index]);
        Az_rhs[index] += dxi[2]*(alpha_Phi_minus_betaj_A_j_interp[CCTK_GFINDEX3D(cctkGH,i,j,k-1)] - alpha_Phi_minus_betaj_A_j_interp[index]);

        CCTK_REAL betax_stencil[5], betay_stencil[5], betaz_stencil[5];
        CCTK_REAL phitilde_stencil[3][5], sqrtg_Ai_stencil[3][2];

        sqrtg_Ai_stencil[0][0] = sqrtg_Ax_interp[index];
        sqrtg_Ai_stencil[1][0] = sqrtg_Ay_interp[index];
        sqrtg_Ai_stencil[2][0] = sqrtg_Az_interp[index];

        sqrtg_Ai_stencil[0][1] = sqrtg_Ax_interp[CCTK_GFINDEX3D(cctkGH,i+1,j,k)];
        sqrtg_Ai_stencil[1][1] = sqrtg_Ay_interp[CCTK_GFINDEX3D(cctkGH,i,j+1,k)];
        sqrtg_Ai_stencil[2][1] = sqrtg_Az_interp[CCTK_GFINDEX3D(cctkGH,i,j,k+1)];

        for(int iter=-2; iter<3; iter++) {
          const int indexx = CCTK_GFINDEX3D(cctkGH,i+iter,j,     k     );
          const int indexy = CCTK_GFINDEX3D(cctkGH,i,     j+iter,k     );
          const int indexz = CCTK_GFINDEX3D(cctkGH,i,     j,     k+iter);
          betax_stencil[iter+2] = betax_interp[indexx];
          betay_stencil[iter+2] = betay_interp[indexy];
          betaz_stencil[iter+2] = betaz_interp[indexz];
          phitilde_stencil[0][iter+2] = phitilde[indexx];
          phitilde_stencil[1][iter+2] = phitilde[indexy];
          phitilde_stencil[2][iter+2] = phitilde[indexz];
        }
        phitilde_rhs[index] = ghl_calculate_phitilde_rhs(dxi, ghl_params->Lorenz_damping_factor, alpha_interp[index], betax_stencil, betay_stencil, betaz_stencil, sqrtg_Ai_stencil, phitilde_stencil);
      }
    }
  }
}
