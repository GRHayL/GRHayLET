#include "GRHayLMHD.h"

void GRHayLMHD_hybrid_calculate_flux_dir_rhs(
      const cGH *restrict cctkGH,
      const int flux_dir,
      const CCTK_REAL **B_center,
      const CCTK_REAL *restrict B_stagger,
      CCTK_REAL **vel_r,
      CCTK_REAL **vel_l,
      CCTK_REAL *restrict cmin,
      CCTK_REAL *restrict cmax) {
  DECLARE_CCTK_ARGUMENTS_GRHayLMHD_hybrid_evaluate_fluxes_rhs;

  const int imin = cctkGH->cctk_nghostzones[0];
  const int jmin = cctkGH->cctk_nghostzones[1];
  const int kmin = cctkGH->cctk_nghostzones[2];
  const int imax = cctkGH->cctk_lsh[0] - cctkGH->cctk_nghostzones[0];
  const int jmax = cctkGH->cctk_lsh[1] - cctkGH->cctk_nghostzones[1];
  const int kmax = cctkGH->cctk_lsh[2] - cctkGH->cctk_nghostzones[2];

  // Function pointer to allow for loop over fluxes and sources
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
      calculate_HLLE_fluxes = ghl_calculate_HLLE_fluxes_dirn0_hybrid;
      break;
    case 1:
      v_flux = vy;
      B_recon[0] = 1;
      B_recon[1] = 2;
      B_recon[2] = 0;
      calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn1;
      calculate_HLLE_fluxes = ghl_calculate_HLLE_fluxes_dirn1_hybrid;
      break;
    case 2:
      v_flux = vz;
      B_recon[0] = 2;
      B_recon[1] = 0;
      B_recon[2] = 1;
      calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn2;
      calculate_HLLE_fluxes = ghl_calculate_HLLE_fluxes_dirn2_hybrid;
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
        CCTK_REAL vx_data[6], vy_data[6], vz_data[6];
        CCTK_REAL vxr, vxl, vyr, vyl, vzr, vzl;

        for(int ind=0; ind<6; ind++) {
          // Stencil from -3 to +2 reconstructs to e.g. i-1/2
          const int stencil = CCTK_GFINDEX3D(cctkGH, i+xdir*(ind-3), j+ydir*(ind-3), k+zdir*(ind-3));
          v_flux_dir[ind] = v_flux[stencil]; // Could be smaller; doesn't use full stencil
          press_stencil[ind] = pressure[stencil];
          vx_data[ind] = vx[stencil];
          vy_data[ind] = vy[stencil];
          vz_data[ind] = vz[stencil];
        }

        double ftilde[2];
        ghl_compute_ftilde(ghl_params, press_stencil, v_flux_dir, ftilde);

        ghl_ppm_reconstruction(ftilde, vx_data, &vxr, &vxl);
        ghl_ppm_reconstruction(ftilde, vy_data, &vyr, &vyl);
        ghl_ppm_reconstruction(ftilde, vz_data, &vzr, &vzl);

        vel_r[0][index] = vxr;
        vel_r[1][index] = vyr;
        vel_r[2][index] = vzr;

        vel_l[0][index] = vxl;
        vel_l[1][index] = vyl;
        vel_l[2][index] = vzl;
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

        ghl_metric_quantities ADM_metric_face;
        GRHayLMHD_interpolate_metric_to_face(
              cctkGH, i, j, k,
              flux_dir, alp,
              betax, betay, betaz,
              gxx, gxy, gxz,
              gyy, gyz, gzz,
              &ADM_metric_face);

        CCTK_REAL rho_stencil[6], press_stencil[6], v_flux_dir[6];
        CCTK_REAL B1_stencil[6], B2_stencil[6];
        ghl_primitive_quantities prims_r, prims_l;

        for(int ind=0; ind<6; ind++) {
          // Stencil from -3 to +2 reconstructs to e.g. i-1/2
          const int stencil  = CCTK_GFINDEX3D(cctkGH, i+xdir*(ind-3), j+ydir*(ind-3), k+zdir*(ind-3));
          v_flux_dir[ind]    = v_flux[stencil]; // Could be smaller; doesn't use full stencil
          rho_stencil[ind]   = rho_b[stencil];
          press_stencil[ind] = pressure[stencil];
          B1_stencil[ind]    = B_center[B_recon[1]][stencil];
          B2_stencil[ind]    = B_center[B_recon[2]][stencil];
        }

        double ftilde[2];
        ghl_compute_ftilde(ghl_params, press_stencil, v_flux_dir, ftilde);

        const CCTK_REAL Gamma = get_Gamma_eff_hybrid(rho_b[index], pressure[index]);
        ghl_ppm_reconstruction_with_steepening(ghl_params, press_stencil, Gamma, ftilde, rho_stencil, &prims_r.rho, &prims_l.rho);

        ghl_ppm_reconstruction(ftilde, press_stencil, &prims_r.press, &prims_l.press);
        ghl_ppm_reconstruction(ftilde, B1_stencil, &prims_r.BU[B_recon[1]], &prims_l.BU[B_recon[1]]);
        ghl_ppm_reconstruction(ftilde, B2_stencil, &prims_r.BU[B_recon[2]], &prims_l.BU[B_recon[2]]);

        // B_stagger is densitized, but B_center is not.
        prims_r.BU[B_recon[0]] = prims_l.BU[B_recon[0]] = B_stagger[indm1]/ADM_metric_face.sqrt_detgamma;

        prims_r.vU[0] = vel_r[0][index];
        prims_r.vU[1] = vel_r[1][index];
        prims_r.vU[2] = vel_r[2][index];

        prims_l.vU[0] = vel_l[0][index];
        prims_l.vU[1] = vel_l[1][index];
        prims_l.vU[2] = vel_l[2][index];

        int speed_limited CCTK_ATTRIBUTE_UNUSED = ghl_limit_v_and_compute_u0(ghl_params, &ADM_metric_face, &prims_r);
        speed_limited = ghl_limit_v_and_compute_u0(ghl_params, &ADM_metric_face, &prims_l);

        ghl_conservative_quantities cons_fluxes;
        calculate_characteristic_speed(&prims_r, &prims_l, ghl_eos, &ADM_metric_face, &cmin[index], &cmax[index]);
        calculate_HLLE_fluxes(&prims_r, &prims_l, ghl_eos, &ADM_metric_face, cmin[index], cmax[index], &cons_fluxes);

        rho_star_flux[index] = cons_fluxes.rho;
        tau_flux     [index] = cons_fluxes.tau;
        Stildex_flux [index] = cons_fluxes.SD[0];
        Stildey_flux [index] = cons_fluxes.SD[1];
        Stildez_flux [index] = cons_fluxes.SD[2];
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
      }
    }
  }
}
