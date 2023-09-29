#include "GRHayLMHD.h"

void GRHayLMHD_tabulated_calculate_flux_dir_rhs(
      const cGH *restrict cctkGH,
      const int flux_dir,
      const CCTK_REAL **B_center,
      const CCTK_REAL *restrict B_stagger,
      CCTK_REAL **vel_r,
      CCTK_REAL **vel_l,
      CCTK_REAL *restrict cmin,
      CCTK_REAL *restrict cmax) {
  DECLARE_CCTK_ARGUMENTS_GRHayLMHD_tabulated_evaluate_fluxes_rhs;

  const int imin = cctkGH->cctk_nghostzones[0];
  const int jmin = cctkGH->cctk_nghostzones[1];
  const int kmin = cctkGH->cctk_nghostzones[2];
  const int imax = cctkGH->cctk_lsh[0] - cctkGH->cctk_nghostzones[0];
  const int jmax = cctkGH->cctk_lsh[1] - cctkGH->cctk_nghostzones[1];
  const int kmax = cctkGH->cctk_lsh[2] - cctkGH->cctk_nghostzones[2];

  // Function pointer to allow for loop over fluxes and sources
  void (*calculate_characteristic_speed)(ghl_primitive_quantities *restrict prims_r,
                                         ghl_primitive_quantities *restrict prims_l,
                                         const ghl_eos_parameters *restrict eos,
                                         const ghl_metric_quantities *restrict ADM_metric_face,
                                         double *cmin, double *cmax);

  void (*calculate_HLLE_fluxes)(ghl_primitive_quantities *restrict prims_r,
                                ghl_primitive_quantities *restrict prims_l,
                                const ghl_eos_parameters *restrict eos,
                                const ghl_metric_quantities *restrict ADM_metric_face,
                                const double cmin,
                                const double cmax,
                                ghl_conservative_quantities *restrict cons_fluxes);

  const int xdir = (flux_dir == 0);
  const int ydir = (flux_dir == 1);
  const int zdir = (flux_dir == 2);

  const CCTK_REAL *v_flux;
  int B_recon[3];
  // Set function pointer to specific function for a given direction
  switch(flux_dir) {
    case 0:
      v_flux = vx;
      B_recon[0] = 0;
      B_recon[1] = 1;
      B_recon[2] = 2;
      calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn0;
      calculate_HLLE_fluxes = ghl_calculate_HLLE_fluxes_dirn0_tabulated;
      break;
    case 1:
      v_flux = vy;
      B_recon[0] = 1;
      B_recon[1] = 2;
      B_recon[2] = 0;
      calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn1;
      calculate_HLLE_fluxes = ghl_calculate_HLLE_fluxes_dirn1_tabulated;
      break;
    case 2:
      v_flux = vz;
      B_recon[0] = 2;
      B_recon[1] = 0;
      B_recon[2] = 1;
      calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn2;
      calculate_HLLE_fluxes = ghl_calculate_HLLE_fluxes_dirn2_tabulated;
      break;
    default:
      CCTK_ERROR("Invalid flux_dir value (not 0, 1, or 2) has been passed to calculate_MHD_rhs.");
  }

  // This loop fills in all the data for reconstructed velocities. This loop is larger
  // than the others because we will need to reconstruct a second time.
  const int vimin = xdir*cctkGH->cctk_nghostzones[0];
  const int vjmin = ydir*cctkGH->cctk_nghostzones[1];
  const int vkmin = zdir*cctkGH->cctk_nghostzones[2];
  const int vimax = cctkGH->cctk_lsh[0] - xdir*(cctkGH->cctk_nghostzones[0] - 1);
  const int vjmax = cctkGH->cctk_lsh[1] - ydir*(cctkGH->cctk_nghostzones[1] - 1);
  const int vkmax = cctkGH->cctk_lsh[2] - zdir*(cctkGH->cctk_nghostzones[2] - 1);

#pragma omp parallel for
  for(int k=vkmin; k<vkmax; k++) {
    for(int j=vjmin; j<vjmax; j++) {
      for(int i=vimin; i<vimax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);

        CCTK_REAL press_stencil[6], v_flux_dir[6];
        CCTK_REAL var_data[3][6], vars_r[3], vars_l[3];

        for(int ind=0; ind<6; ind++) {
          // Stencil from -3 to +2 reconstructs to e.g. i-1/2
          const int stencil = CCTK_GFINDEX3D(cctkGH, i+xdir*(ind-3), j+ydir*(ind-3), k+zdir*(ind-3));
          v_flux_dir[ind] = v_flux[stencil]; // Could be smaller; doesn't use full stencil
          press_stencil[ind] = pressure[stencil];
          var_data[0][ind] = vx[stencil];
          var_data[1][ind] = vy[stencil];
          var_data[2][ind] = vz[stencil];
        }

        // Compute Gamma
        const CCTK_REAL Gamma = get_Gamma_eff_tabulated(rho_b[index], pressure[index]);

        ghl_ppm_no_rho_P(
              ghl_params, press_stencil, var_data,
              3, v_flux_dir, Gamma,
              vars_r, vars_l);

        vel_r[0][index] = vars_r[0];
        vel_r[1][index] = vars_r[1];
        vel_r[2][index] = vars_r[2];

        vel_l[0][index] = vars_l[0];
        vel_l[1][index] = vars_l[1];
        vel_l[2][index] = vars_l[2];
      }
    }
  }

  // This loop includes 1 ghostzone because the RHS calculation for e.g. the x direction
  // requires (i,j,k) and (i+1,j,k); if cmin/max weren't also needed for A_i, we could
  // technically have the loop only go 1 extra point in the flux_dir direction
#pragma omp parallel for
  for(int k=kmin; k<kmax+1; k++) {
    for(int j=jmin; j<jmax+1; j++) {
      for(int i=imin; i<imax+1; i++) {
        const int indm1 = CCTK_GFINDEX3D(cctkGH, i-xdir, j-ydir, k-zdir);
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);

        CCTK_REAL rho_stencil[6], press_stencil[6], v_flux_dir[6];
        CCTK_REAL rhor, rhol, pressr, pressl, B_r[3], B_l[3];
        CCTK_REAL others_stencil[3][6], others_r[3], others_l[3];

        for(int ind=0; ind<6; ind++) {
          // Stencil from -3 to +2 reconstructs to e.g. i-1/2
          const int stencil      = CCTK_GFINDEX3D(cctkGH, i+xdir*(ind-3), j+ydir*(ind-3), k+zdir*(ind-3));
          v_flux_dir[ind]        = v_flux[stencil]; // Could be smaller; doesn't use full stencil
          rho_stencil[ind]       = rho_b[stencil];
          press_stencil[ind]     = pressure[stencil];
          others_stencil[0][ind] = B_center[B_recon[1]][stencil];
          others_stencil[1][ind] = B_center[B_recon[2]][stencil];
          others_stencil[2][ind] = Y_e[stencil];
        }

        // Compute Gamma
        const CCTK_REAL Gamma = get_Gamma_eff_tabulated(rho_b[index], pressure[index]);

        ghl_ppm(
              ghl_params, rho_stencil,
              press_stencil, others_stencil,
              3, v_flux_dir, Gamma,
              &rhor, &rhol,
              &pressr, &pressl,
              others_r, others_l);

        B_r[B_recon[1]] = others_r[0];
        B_r[B_recon[2]] = others_r[1];

        B_l[B_recon[1]] = others_l[0];
        B_l[B_recon[2]] = others_l[1];

        ghl_metric_quantities ADM_metric_face;
        GRHayLMHD_interpolate_metric_to_face(
              cctkGH, i, j, k,
              flux_dir, alp,
              betax, betay, betaz,
              gxx, gxy, gxz,
              gyy, gyz, gzz,
              &ADM_metric_face);

        // B_stagger is densitized, but B_center is not.
        B_r[B_recon[0]] = B_stagger[indm1]/ADM_metric_face.sqrt_detgamma;
        B_l[B_recon[0]] = B_stagger[indm1]/ADM_metric_face.sqrt_detgamma;

        ghl_primitive_quantities prims_r, prims_l;
        prims_r.rho   = rhor;
        prims_r.press = pressr;
        prims_r.vU[0] = vel_r[0][index];
        prims_r.vU[1] = vel_r[1][index];
        prims_r.vU[2] = vel_r[2][index];
        prims_r.BU[0] = B_r[0];
        prims_r.BU[1] = B_r[1];
        prims_r.BU[2] = B_r[2];
        prims_r.Y_e = others_r[2];
        prims_r.temperature = temperature[index];
        int speed_limited CCTK_ATTRIBUTE_UNUSED = ghl_limit_v_and_compute_u0(ghl_params, &ADM_metric_face, &prims_r);

        prims_l.rho   = rhol;
        prims_l.press = pressl;
        prims_l.vU[0] = vel_l[0][index];
        prims_l.vU[1] = vel_l[1][index];
        prims_l.vU[2] = vel_l[2][index];
        prims_l.BU[0] = B_l[0];
        prims_l.BU[1] = B_l[1];
        prims_l.BU[2] = B_l[2];
        prims_l.Y_e = others_l[2];
        prims_l.temperature = temperature[index];
        speed_limited = ghl_limit_v_and_compute_u0(ghl_params, &ADM_metric_face, &prims_l);

        ghl_conservative_quantities cons_fluxes;
        calculate_characteristic_speed(&prims_r, &prims_l, ghl_eos, &ADM_metric_face, &cmin[index], &cmax[index]);
        calculate_HLLE_fluxes(&prims_r, &prims_l, ghl_eos, &ADM_metric_face, cmin[index], cmax[index], &cons_fluxes);

        rho_star_flux[index] = cons_fluxes.rho;
        tau_flux     [index] = cons_fluxes.tau;
        Stildex_flux [index] = cons_fluxes.SD[0];
        Stildey_flux [index] = cons_fluxes.SD[1];
        Stildez_flux [index] = cons_fluxes.SD[2];
        Ye_star_flux [index] = cons_fluxes.Y_e;
      }
    }
  }

  const CCTK_REAL dxi = 1.0/CCTK_DELTA_SPACE(flux_dir);

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
        Ye_star_rhs [index] += dxi*(Ye_star_flux [index] - Ye_star_flux [indp1]);
      }
    }
  }
}
