#include "GRHayLMHD.h"

static inline double get_Gamma_eff_hybrid(
      const double rho_in,
      const double press_in) {
  double K, Gamma;
  ghl_hybrid_get_K_and_Gamma(ghl_eos, rho_in, &K, &Gamma);
  const double P_cold = K*pow(rho_in, Gamma);
  return ghl_eos->Gamma_th + (Gamma - ghl_eos->Gamma_th)*P_cold/press_in;
}

static inline double get_Gamma_eff_tabulated(
      const double rho_in,
      const double press_in) {
  return 1.0;
}

static double (*get_Gamma_eff)(const double, const double) = &get_Gamma_eff_hybrid;

void GRHayLMHD_reconstruction_loop(const cGH *restrict cctkGH, const int flux_dir, const int num_vars,
                         const int *restrict var_indices,
                         const double *rho_b,
                         const double *pressure,
                         const double *v_flux,
                         const double **in_prims,
                         double **out_prims_r,
                         double **out_prims_l) {

  const int xdir = (flux_dir == 0);
  const int ydir = (flux_dir == 1);
  const int zdir = (flux_dir == 2);

  const int imin = cctkGH->cctk_nghostzones[0];
  const int imax = (cctkGH->cctk_lsh[0] - cctkGH->cctk_nghostzones[0]) + 1;
  const int jmin = cctkGH->cctk_nghostzones[1];
  const int jmax = (cctkGH->cctk_lsh[1] - cctkGH->cctk_nghostzones[1]) + 1;
  const int kmin = cctkGH->cctk_nghostzones[2];
  const int kmax = (cctkGH->cctk_lsh[2] - cctkGH->cctk_nghostzones[2]) + 1;

  if( ghl_eos->eos_type == ghl_eos_tabulated )
    get_Gamma_eff = &get_Gamma_eff_tabulated;

#pragma omp parallel for
  for(int k=kmin; k<kmax; k++) {
    for(int j=jmin; j<jmax; j++) {
      for(int i=imin; i<imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);
        double press_stencil[6], v_flux_stencil[6];
        double var_data[num_vars][6], vars_r[num_vars], vars_l[num_vars];

        for(int ind=0; ind<6; ind++) {
          // Stencil from -3 to +2 reconstructs to e.g. i-1/2
          const int stencil = CCTK_GFINDEX3D(cctkGH, i+xdir*(ind-3), j+ydir*(ind-3), k+zdir*(ind-3));
          v_flux_stencil[ind] = v_flux[stencil]; // Could be smaller; doesn't use full stencil
          press_stencil[ind] = pressure[stencil];
          for(int var=0; var<num_vars; var++) {
            var_data[var][ind] = in_prims[var_indices[var]][stencil];
          }
        }

        // Compute Gamma
        const double Gamma = get_Gamma_eff(rho_b[index], pressure[index]);

        ghl_ppm_no_rho_P(
              press_stencil, var_data,
              num_vars, v_flux_stencil, Gamma,
              vars_r, vars_l);

        for(int var=0; var<num_vars; var++) {
          out_prims_r[var_indices[var]][index] = vars_r[var];
          out_prims_l[var_indices[var]][index] = vars_l[var];
        }
      }
    }
  }
}
