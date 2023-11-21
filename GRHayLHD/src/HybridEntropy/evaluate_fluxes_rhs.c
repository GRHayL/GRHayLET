#include "GRHayLHD.h"

static inline CCTK_REAL get_Gamma_eff(
      const CCTK_REAL rho_in,
      const CCTK_REAL press_in) {
  CCTK_REAL K, Gamma;
  ghl_hybrid_get_K_and_Gamma(ghl_eos, rho_in, &K, &Gamma);
  const CCTK_REAL P_cold = K*pow(rho_in, Gamma);
  return ghl_eos->Gamma_th + (Gamma - ghl_eos->Gamma_th)*P_cold/press_in;
}

void GRHayLHD_hybrid_entropy_evaluate_fluxes_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLHD_hybrid_entropy_evaluate_fluxes_rhs;
  DECLARE_CCTK_PARAMETERS;

  const int imin = cctkGH->cctk_nghostzones[0];
  const int jmin = cctkGH->cctk_nghostzones[1];
  const int kmin = cctkGH->cctk_nghostzones[2];
  const int imax = cctkGH->cctk_lsh[0] - cctkGH->cctk_nghostzones[0];
  const int jmax = cctkGH->cctk_lsh[1] - cctkGH->cctk_nghostzones[1];
  const int kmax = cctkGH->cctk_lsh[2] - cctkGH->cctk_nghostzones[2];

  void (*calculate_characteristic_speed)(
        ghl_primitive_quantities *restrict prims_r,
        ghl_primitive_quantities *restrict prims_l,
        const ghl_eos_parameters *restrict eos,
        const ghl_metric_quantities *restrict ADM_metric_face,
        double *cmin, double *cmax);

  void (*calculate_HLLE_fluxes)(
        ghl_primitive_quantities *restrict prims_r,
        ghl_primitive_quantities *restrict prims_l,
        const ghl_eos_parameters *restrict eos,
        const ghl_metric_quantities *restrict ADM_metric_face,
        const double cmin,
        const double cmax,
        ghl_conservative_quantities *restrict cons_fluxes);

  for(int flux_dir=0; flux_dir<3; flux_dir++) {
    const int xdir = (flux_dir == 0);
    const int ydir = (flux_dir == 1);
    const int zdir = (flux_dir == 2);
    const double *v_flux_dir;

    // Set function pointer to specific function for a given direction
    switch(flux_dir) {
      case 0:
        v_flux_dir = vx;
        calculate_characteristic_speed = ghl_calculate_characteristic_speed_dirn0;
        calculate_HLLE_fluxes = ghl_calculate_HLLE_fluxes_dirn0_hybrid_entropy;
        break;
      case 1:
        v_flux_dir = vy;
        calculate_characteristic_speed = ghl_calculate_characteristic_speed_dirn1;
        calculate_HLLE_fluxes = ghl_calculate_HLLE_fluxes_dirn1_hybrid_entropy;
        break;
      case 2:
        v_flux_dir = vz;
        calculate_characteristic_speed = ghl_calculate_characteristic_speed_dirn2;
        calculate_HLLE_fluxes = ghl_calculate_HLLE_fluxes_dirn2_hybrid_entropy;
        break;
      default:
        CCTK_ERROR("Invalid flux_dir value (not 0, 1, or 2) has been passed to calculate_MHD_rhs.");
    }

    // This loop includes 1 ghostzone because the RHS calculation for e.g. the x direction
    // requires (i,j,k) and (i+1,j,k)
#pragma omp parallel for
    for(int k=kmin; k<kmax+1; k++) {
      for(int j=jmin; j<jmax+1; j++) {
        for(int i=imin; i<imax+1; i++) {
          const int index = CCTK_GFINDEX3D(cctkGH, i, j ,k);

          ghl_metric_quantities ADM_metric_face;
          GRHayLHD_interpolate_metric_to_face(
                cctkGH, i, j, k,
                flux_dir, alp,
                betax, betay, betaz,
                gxx, gxy, gxz,
                gyy, gyz, gzz,
                &ADM_metric_face);

          CCTK_REAL rho_stencil[6], press_stencil[6], v_flux[6];
          CCTK_REAL vx_stencil[6], vy_stencil[6], vz_stencil[6];
          CCTK_REAL ent_stencil[6];
          ghl_primitive_quantities prims_r, prims_l;

          for(int ind=0; ind<6; ind++) {
            // Stencil from -3 to +2 reconstructs to e.g. i-1/2
            const int stencil  = CCTK_GFINDEX3D(cctkGH, i+xdir*(ind-3), j+ydir*(ind-3), k+zdir*(ind-3));
            v_flux[ind]        = v_flux_dir[stencil]; // Could be smaller; doesn't use full stencil
            rho_stencil[ind]   = rho[stencil];
            press_stencil[ind] = press[stencil];
            vx_stencil[ind]    = vx[stencil];
            vy_stencil[ind]    = vy[stencil];
            vz_stencil[ind]    = vz[stencil];
            ent_stencil[ind]   = entropy[stencil];
          }

          CCTK_REAL ftilde[2];
          ghl_compute_ftilde(ghl_params, press_stencil, v_flux, ftilde);

          const CCTK_REAL Gamma = get_Gamma_eff(rho[index], press[index]);
          ghl_ppm_reconstruction_with_steepening(ghl_params, press_stencil, Gamma, ftilde, rho_stencil, &prims_r.rho, &prims_l.rho);

          ghl_ppm_reconstruction(ftilde, press_stencil, &prims_r.press, &prims_l.press);
          ghl_ppm_reconstruction(ftilde, vx_stencil, &prims_r.vU[0], &prims_l.vU[0]);
          ghl_ppm_reconstruction(ftilde, vy_stencil, &prims_r.vU[1], &prims_l.vU[1]);
          ghl_ppm_reconstruction(ftilde, vz_stencil, &prims_r.vU[2], &prims_l.vU[2]);
          ghl_ppm_reconstruction(ftilde, ent_stencil, &prims_r.entropy, &prims_l.entropy);

          prims_r.BU[0] = prims_r.BU[1] = prims_r.BU[2] = 0.0;
          prims_l.BU[0] = prims_l.BU[1] = prims_l.BU[2] = 0.0;

          int speed_limited CCTK_ATTRIBUTE_UNUSED = ghl_limit_v_and_compute_u0(ghl_params, &ADM_metric_face, &prims_r);
          speed_limited = ghl_limit_v_and_compute_u0(ghl_params, &ADM_metric_face, &prims_l);

          CCTK_REAL cmin, cmax;
          ghl_conservative_quantities cons_fluxes;
          calculate_characteristic_speed(&prims_r, &prims_l, ghl_eos, &ADM_metric_face, &cmin, &cmax);
          calculate_HLLE_fluxes(&prims_r, &prims_l, ghl_eos, &ADM_metric_face, cmin, cmax, &cons_fluxes);

          rho_star_flux[index] = cons_fluxes.rho;
          tau_flux     [index] = cons_fluxes.tau;
          Stildex_flux [index] = cons_fluxes.SD[0];
          Stildey_flux [index] = cons_fluxes.SD[1];
          Stildez_flux [index] = cons_fluxes.SD[2];
          ent_star_flux[index] = cons_fluxes.entropy;
        }
      }
    }

    CCTK_REAL dxi = 1.0/CCTK_DELTA_SPACE(flux_dir);

#pragma omp parallel for
    for(int k=kmin; k<kmax; k++) {
      for(int j=jmin; j<jmax; j++) {
        for(int i=imin; i<imax; i++) {
          const int index = CCTK_GFINDEX3D(cctkGH, i, j ,k);
          const int indp1 = CCTK_GFINDEX3D(cctkGH, i+xdir, j+ydir, k+zdir);

          rho_star_rhs[index] += dxi*(rho_star_flux[index] - rho_star_flux[indp1]);
          tau_rhs     [index] += dxi*(tau_flux     [index] - tau_flux     [indp1]);
          Stildex_rhs [index] += dxi*(Stildex_flux [index] - Stildex_flux [indp1]);
          Stildey_rhs [index] += dxi*(Stildey_flux [index] - Stildey_flux [indp1]);
          Stildez_rhs [index] += dxi*(Stildez_flux [index] - Stildez_flux [indp1]);
          ent_star_rhs[index] += dxi*(ent_star_flux[index] - ent_star_flux[indp1]);
        }
      }
    }
  }
}
