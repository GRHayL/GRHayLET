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

void reconstruct_v_edges(
      const cGH *restrict cctkGH,
      const int flux_dir,
      const int imin,
      const int imax,
      const int jmin,
      const int jmax,
      const int kmin,
      const int kmax,
      const CCTK_REAL *restrict rho_b,
      const CCTK_REAL *restrict pressure,
      const CCTK_REAL *v_flux,
      const CCTK_REAL *vx,
      const CCTK_REAL *vy,
      const CCTK_REAL *vz,
      CCTK_REAL **vel_r,
      CCTK_REAL **vel_l);

// We could reduce computational cost by introducing grid functions for the
// metric face interpolations (used in both loops
void GRHayLMHD_calculate_MHD_dirn_rhs(
      const cGH *restrict cctkGH,
      const int flux_dir,
      const CCTK_REAL *restrict dX,
      const CCTK_REAL *restrict lapse,
      const CCTK_REAL *restrict betax,
      const CCTK_REAL *restrict betay,
      const CCTK_REAL *restrict betaz,
      const CCTK_REAL *restrict gxx,
      const CCTK_REAL *restrict gxy,
      const CCTK_REAL *restrict gxz,
      const CCTK_REAL *restrict gyy,
      const CCTK_REAL *restrict gyz,
      const CCTK_REAL *restrict gzz,
      const CCTK_REAL *restrict rho_b,
      const CCTK_REAL *restrict pressure,
      const CCTK_REAL *restrict ent,
      const CCTK_REAL *restrict Ye,
      const CCTK_REAL *restrict temp,
      const CCTK_REAL *vx,
      const CCTK_REAL *vy,
      const CCTK_REAL *vz,
      const CCTK_REAL **B_center,
      const CCTK_REAL *restrict B_stagger,
      CCTK_REAL **vel_r,
      CCTK_REAL **vel_l,
      CCTK_REAL *restrict cmin,
      CCTK_REAL *restrict cmax,
      CCTK_REAL *restrict rho_star_flux,
      CCTK_REAL *restrict tau_flux,
      CCTK_REAL *restrict Stildex_flux,
      CCTK_REAL *restrict Stildey_flux,
      CCTK_REAL *restrict Stildez_flux,
      CCTK_REAL *restrict ent_star_flux,
      CCTK_REAL *restrict Ye_star_flux,
      CCTK_REAL *restrict rho_star_rhs,
      CCTK_REAL *restrict tau_rhs,
      CCTK_REAL *restrict Stildex_rhs,
      CCTK_REAL *restrict Stildey_rhs,
      CCTK_REAL *restrict Stildez_rhs,
      CCTK_REAL *restrict ent_star_rhs,
      CCTK_REAL *restrict Ye_star_rhs) {

  const CCTK_REAL dxi = 1.0/dX[flux_dir];
  const CCTK_REAL poison = 0.0/0.0;

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

  void (*calculate_source_terms)(ghl_primitive_quantities *restrict prims,
                                 const ghl_eos_parameters *restrict eos,
                                 const ghl_metric_quantities *restrict ADM_metric,
                                 const ghl_metric_quantities *restrict metric_derivs,
                                 ghl_conservative_quantities *restrict cons_sources);

  void (*calculate_HLLE_fluxes)(ghl_primitive_quantities *restrict prims_r,
                                ghl_primitive_quantities *restrict prims_l,
                                const ghl_eos_parameters *restrict eos,
                                const ghl_metric_quantities *restrict ADM_metric_face,
                                const double cmin,
                                const double cmax,
                                ghl_conservative_quantities *restrict cons_fluxes);

  void (*calculate_HLLE_fluxes_dirn0)(ghl_primitive_quantities *restrict prims_r,
                                      ghl_primitive_quantities *restrict prims_l,
                                      const ghl_eos_parameters *restrict eos,
                                      const ghl_metric_quantities *restrict ADM_metric_face,
                                      const double cmin,
                                      const double cmax,
                                      ghl_conservative_quantities *restrict cons_fluxes);

  void (*calculate_HLLE_fluxes_dirn1)(ghl_primitive_quantities *restrict prims_r,
                                      ghl_primitive_quantities *restrict prims_l,
                                      const ghl_eos_parameters *restrict eos,
                                      const ghl_metric_quantities *restrict ADM_metric_face,
                                      const double cmin,
                                      const double cmax,
                                      ghl_conservative_quantities *restrict cons_fluxes);

  void (*calculate_HLLE_fluxes_dirn2)(ghl_primitive_quantities *restrict prims_r,
                                      ghl_primitive_quantities *restrict prims_l,
                                      const ghl_eos_parameters *restrict eos,
                                      const ghl_metric_quantities *restrict ADM_metric_face,
                                      const double cmin,
                                      const double cmax,
                                      ghl_conservative_quantities *restrict cons_fluxes);

  if( ghl_eos->eos_type == ghl_eos_hybrid ) {
    if( ghl_params->evolve_entropy ) {
      calculate_HLLE_fluxes_dirn0 = &ghl_calculate_HLLE_fluxes_dirn0_hybrid_entropy;
      calculate_HLLE_fluxes_dirn1 = &ghl_calculate_HLLE_fluxes_dirn1_hybrid_entropy;
      calculate_HLLE_fluxes_dirn2 = &ghl_calculate_HLLE_fluxes_dirn2_hybrid_entropy;
    }
    else {
      calculate_HLLE_fluxes_dirn0 = &ghl_calculate_HLLE_fluxes_dirn0_hybrid;
      calculate_HLLE_fluxes_dirn1 = &ghl_calculate_HLLE_fluxes_dirn1_hybrid;
      calculate_HLLE_fluxes_dirn2 = &ghl_calculate_HLLE_fluxes_dirn2_hybrid;
    }
  }
  else {
    if( ghl_params->evolve_entropy ) {
      calculate_HLLE_fluxes_dirn0 = &ghl_calculate_HLLE_fluxes_dirn0_tabulated_entropy;
      calculate_HLLE_fluxes_dirn1 = &ghl_calculate_HLLE_fluxes_dirn1_tabulated_entropy;
      calculate_HLLE_fluxes_dirn2 = &ghl_calculate_HLLE_fluxes_dirn2_tabulated_entropy;
    }
    else {
      calculate_HLLE_fluxes_dirn0 = &ghl_calculate_HLLE_fluxes_dirn0_tabulated;
      calculate_HLLE_fluxes_dirn1 = &ghl_calculate_HLLE_fluxes_dirn1_tabulated;
      calculate_HLLE_fluxes_dirn2 = &ghl_calculate_HLLE_fluxes_dirn2_tabulated;
    }
  }

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
      calculate_HLLE_fluxes = calculate_HLLE_fluxes_dirn0;
      calculate_source_terms = &ghl_calculate_source_terms_dirn0;
      break;
    case 1:
      v_flux = vy;
      B_recon[0] = 1;
      B_recon[1] = 2;
      B_recon[2] = 0;
      calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn1;
      calculate_HLLE_fluxes = calculate_HLLE_fluxes_dirn1;
      calculate_source_terms = &ghl_calculate_source_terms_dirn1;
      break;
    case 2:
      v_flux = vz;
      B_recon[0] = 2;
      B_recon[1] = 0;
      B_recon[2] = 1;
      calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn2;
      calculate_HLLE_fluxes = calculate_HLLE_fluxes_dirn2;
      calculate_source_terms = &ghl_calculate_source_terms_dirn2;
      break;
    default:
      CCTK_ERROR("Invalid flux_dir value (not 0, 1, or 2) has been passed to calculate_MHD_rhs.");
  }

  // Count number of additional reconstructed variables
  const int num_others = 5 + ghl_params->evolve_entropy + (ghl_eos->eos_type == ghl_eos_tabulated);

  // If using the entropy, it should be the first reconstructed variable
  // after the three velocities
  const int ent_index = 5 + !ghl_params->evolve_entropy;

  // If not using the entropy, then Ye should be the first reconstructed
  // variable after the three velocities
  const int Ye_index = 6 - !ghl_params->evolve_entropy;

  if( ghl_eos->eos_type == ghl_eos_tabulated )
    get_Gamma_eff = &get_Gamma_eff_tabulated;

#pragma omp parallel for
  for(int k=0; k<cctkGH->cctk_lsh[2]; k++) {
    for(int j=0; j<cctkGH->cctk_lsh[1]; j++) {
      for(int i=0; i<cctkGH->cctk_lsh[0]; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);
        vel_r[0][index] = poison;
        vel_r[1][index] = poison;
        vel_r[2][index] = poison;

        vel_l[0][index] = poison;
        vel_l[1][index] = poison;
        vel_l[2][index] = poison;
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
        const int indm1 = CCTK_GFINDEX3D(cctkGH, i-xdir, j-ydir, k-zdir); /* indexim1=0 when i=0 */
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);

        CCTK_REAL rho_stencil[6], press_stencil[6], v_flux_dir[6];
        CCTK_REAL rhor, rhol, pressr, pressl, B_r[3], B_l[3];
        CCTK_REAL others_stencil[7][6], others_r[7], others_l[7];

        for(int ind=0; ind<6; ind++) {
          // Stencil from -3 to +2 reconstructs to e.g. i-1/2
          const int stencil              = CCTK_GFINDEX3D(cctkGH, i+xdir*(ind-3), j+ydir*(ind-3), k+zdir*(ind-3));
          v_flux_dir[ind]                = v_flux[stencil]; // Could be smaller; doesn't use full stencil
          rho_stencil[ind]               = rho_b[stencil];
          press_stencil[ind]             = pressure[stencil];
          others_stencil[0][ind]         = vx[stencil];
          others_stencil[1][ind]         = vy[stencil];
          others_stencil[2][ind]         = vz[stencil];
          others_stencil[3][ind]         = B_center[B_recon[1]][stencil];
          others_stencil[4][ind]         = B_center[B_recon[2]][stencil];
          others_stencil[ent_index][ind] = ent[stencil];
          others_stencil[Ye_index ][ind] = Ye[stencil];
        }

        // Compute Gamma
        const CCTK_REAL Gamma = get_Gamma_eff(rho_b[index], pressure[index]);

        ghl_ppm(
              rho_stencil, press_stencil, others_stencil,
              num_others, v_flux_dir, Gamma,
              &rhor, &rhol, &pressr, &pressl, others_r, others_l);

        vel_r[0][index] = others_r[0];
        vel_r[1][index] = others_r[1];
        vel_r[2][index] = others_r[2];

        vel_l[0][index] = others_l[0];
        vel_l[1][index] = others_l[1];
        vel_l[2][index] = others_l[2];

        B_r[B_recon[1]] = others_r[3];
        B_r[B_recon[2]] = others_r[4];

        B_l[B_recon[1]] = others_l[3];
        B_l[B_recon[2]] = others_l[4];

        ghl_metric_quantities ADM_metric_face;
        GRHayLMHD_interpolate_metric_to_face(
              cctkGH, i, j, k,
              flux_dir, lapse,
              betax, betay, betaz,
              gxx, gxy, gxz,
              gyy, gyz, gzz,
              &ADM_metric_face);

        // B_stagger is densitized, but B_center is not.
        B_r[B_recon[0]] = B_stagger[indm1]/ADM_metric_face.sqrt_detgamma;
        B_l[B_recon[0]] = B_stagger[indm1]/ADM_metric_face.sqrt_detgamma;

        ghl_primitive_quantities prims_r, prims_l;
        ghl_initialize_primitives(
              rhor, pressr, poison,
              vel_r[0][index], vel_r[1][index], vel_r[2][index],
              B_r[0], B_r[1], B_r[2],
              others_r[ent_index], others_r[Ye_index], temp[index],
              &prims_r);

        ghl_initialize_primitives(
              rhol, pressl, poison,
              vel_l[0][index], vel_l[1][index], vel_l[2][index],
              B_l[0], B_l[1], B_l[2],
              others_l[ent_index], others_l[Ye_index], temp[index],
              &prims_l);

        int speed_limited CCTK_ATTRIBUTE_UNUSED = ghl_limit_v_and_compute_u0(ghl_eos, &ADM_metric_face, &prims_r);
        speed_limited = ghl_limit_v_and_compute_u0(ghl_eos, &ADM_metric_face, &prims_l);

        ghl_conservative_quantities cons_fluxes;
        calculate_characteristic_speed(&prims_r, &prims_l, ghl_eos, &ADM_metric_face, &cmin[index], &cmax[index]);
        calculate_HLLE_fluxes(&prims_r, &prims_l, ghl_eos, &ADM_metric_face, cmin[index], cmax[index], &cons_fluxes);

        rho_star_flux[index] = cons_fluxes.rho;
        tau_flux     [index] = cons_fluxes.tau;
        Stildex_flux [index] = cons_fluxes.SD[0];
        Stildey_flux [index] = cons_fluxes.SD[1];
        Stildez_flux [index] = cons_fluxes.SD[2];
        ent_star_flux[index] = cons_fluxes.entropy;
        Ye_star_flux [index] = cons_fluxes.Y_e;
      }
    }
  }

  /*
     The following loops fill in the outer ghost zones for the
     reconstructed velocities. To compute the flux terms for A_RHS, we
     need to reconstruct the velocities a second time. This means we
     will need to do another reconstruction stencil in the interior.
     However, the previous loop only provides the minimum loop size.
     We could needlessly compute lots of extra math in the previous
     loop to simplify the code, but instead we loop over the extra
     ghost zone points here and just reconstruct the velocities.

     Example: consider the x direction. If we reconstruct in x first,
     then we may want to reconstruct in y or z. Then, the full grid
     we want is
        i: loop over interior + 1 one upper ghost zone (due to different centerings)
        j: loop over all points
        k: loop over all points
     We can visualize this with a box

              ------------
              |          |
              | -------- |
              | |xxxxxx| |
              | |xxxxxx| |
              | |xxxxxx| |
           ^  | -------- |
           |  |          |
           k  ------------
           j-->

     where the 'x' positions have already been filled.
     We can do several loops to fill the remainder:
     1) loop k @ (0, kmin), j @ (0,lsh[1]), i @ (imin,imax+1)
        loop shown with w
     2) loop k @ (kmax+1, lsh[2]), j @ (0,lsh[1]), i @ (imin,imax+1)
        loop shown with y
     3) loop k @ (kmin, kmax+1), j @ (0,jmin), i @ (imin,imax+1)
        loop shown with z
     4) loop k @ (kmin, kmax+1), j @ (jmax+1,lsh[1]), i @ (imin,imax+1)
        loop shown with q
              ------------
              |yyyyyyyyyy|
              |z--------q|
              |z|xxxxxx|q|
              |z|xxxxxx|q|
              |z|xxxxxx|q|
              |z--------q|
              |wwwwwwwwww|
              ------------
  */
  if(xdir) {
    reconstruct_v_edges(cctkGH, flux_dir, // w zone
          imin, imax+1,
          0, cctkGH->cctk_lsh[1],
          0, kmin,
          rho_b, pressure, v_flux, vx, vy, vz, vel_r, vel_l);

    reconstruct_v_edges(cctkGH, flux_dir, // y zone
          imin, imax+1,
          0, cctkGH->cctk_lsh[1],
          kmax+1, cctkGH->cctk_lsh[2],
          rho_b, pressure, v_flux, vx, vy, vz, vel_r, vel_l);

    reconstruct_v_edges(cctkGH, flux_dir, // z zone
          imin, imax+1,
          0, jmin,
          kmin, kmax+1,
          rho_b, pressure, v_flux, vx, vy, vz, vel_r, vel_l);

    reconstruct_v_edges(cctkGH, flux_dir,  // q zone
          imin, imax+1,
          jmax+1, cctkGH->cctk_lsh[1],
          kmin, kmax+1,
          rho_b, pressure, v_flux, vx, vy, vz, vel_r, vel_l);
  } else if(ydir) {
    reconstruct_v_edges(cctkGH, flux_dir, // w zone
          0, cctkGH->cctk_lsh[0],
          jmin, jmax+1,
          0, kmin,
          rho_b, pressure, v_flux, vx, vy, vz, vel_r, vel_l);

    reconstruct_v_edges(cctkGH, flux_dir, // y zone
          0, cctkGH->cctk_lsh[0],
          jmin, jmax+1,
          kmax+1, cctkGH->cctk_lsh[2],
          rho_b, pressure, v_flux, vx, vy, vz, vel_r, vel_l);

    reconstruct_v_edges(cctkGH, flux_dir, // z zone
          0, imin,
          jmin, jmax+1,
          kmin, kmax+1,
          rho_b, pressure, v_flux, vx, vy, vz, vel_r, vel_l);

    reconstruct_v_edges(cctkGH, flux_dir,  // q zone
          imax+1, cctkGH->cctk_lsh[0],
          jmin, jmax+1,
          kmin, kmax+1,
          rho_b, pressure, v_flux, vx, vy, vz, vel_r, vel_l);
  } else if(zdir) {
    reconstruct_v_edges(cctkGH, flux_dir, // w zone
          0, cctkGH->cctk_lsh[0],
          0, jmin,
          kmin, kmax+1,
          rho_b, pressure, v_flux, vx, vy, vz, vel_r, vel_l);

    reconstruct_v_edges(cctkGH, flux_dir, // y zone
          0, cctkGH->cctk_lsh[0],
          jmax+1, cctkGH->cctk_lsh[1],
          kmin, kmax+1,
          rho_b, pressure, v_flux, vx, vy, vz, vel_r, vel_l);

    reconstruct_v_edges(cctkGH, flux_dir, // z zone
          0, imin,
          jmin, jmax+1,
          kmin, kmax+1,
          rho_b, pressure, v_flux, vx, vy, vz, vel_r, vel_l);

    reconstruct_v_edges(cctkGH, flux_dir,  // q zone
          imax+1, cctkGH->cctk_lsh[0],
          jmin, jmax+1,
          kmin, kmax+1,
          rho_b, pressure, v_flux, vx, vy, vz, vel_r, vel_l);
  }

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
        Ye_star_rhs [index] += dxi*(Ye_star_flux [index] - Ye_star_flux [indp1]);

        ghl_metric_quantities ADM_metric;
        ghl_initialize_metric(
              lapse[index],
              betax[index], betay[index], betaz[index],
              gxx[index], gxy[index], gxz[index],
              gyy[index], gyz[index], gzz[index],
              &ADM_metric);

        ghl_primitive_quantities prims;
        ghl_initialize_primitives(
              rho_b[index], pressure[index], poison,
              vx[index], vy[index], vz[index],
              B_center[0][index], B_center[1][index], B_center[2][index],
              ent[index], Ye[index], temp[index],
              &prims);

        const int speed_limited CCTK_ATTRIBUTE_UNUSED = ghl_limit_v_and_compute_u0(
              ghl_eos, &ADM_metric, &prims);

        ghl_metric_quantities ADM_metric_derivs;
        GRHayLMHD_compute_metric_derivs(
              cctkGH, i, j, k,
              flux_dir, dxi, lapse,
              betax, betay, betaz,
              gxx, gxy, gxz,
              gyy, gyz, gzz,
              &ADM_metric_derivs);

        ghl_conservative_quantities cons_source;
        cons_source.tau = 0.0;
        cons_source.SD[0] = 0.0;
        cons_source.SD[1] = 0.0;
        cons_source.SD[2] = 0.0;

        calculate_source_terms(&prims, ghl_eos, &ADM_metric, &ADM_metric_derivs, &cons_source);
        tau_rhs[index]     += cons_source.tau;
        Stildex_rhs[index] += cons_source.SD[0];
        Stildey_rhs[index] += cons_source.SD[1];
        Stildez_rhs[index] += cons_source.SD[2];
      }
    }
  }
}

void reconstruct_v_edges(
      const cGH *restrict cctkGH,
      const int flux_dir,
      const int imin,
      const int imax,
      const int jmin,
      const int jmax,
      const int kmin,
      const int kmax,
      const CCTK_REAL *restrict rho_b,
      const CCTK_REAL *restrict pressure,
      const CCTK_REAL *v_flux,
      const CCTK_REAL *vx,
      const CCTK_REAL *vy,
      const CCTK_REAL *vz,
      CCTK_REAL **vel_r,
      CCTK_REAL **vel_l) {

  const int xdir = (flux_dir == 0);
  const int ydir = (flux_dir == 1);
  const int zdir = (flux_dir == 2);

#pragma omp parallel for
  for(int k=kmin; k<kmax; k++) {
    for(int j=jmin; j<jmax; j++) {
      for(int i=imin; i<imax; i++) {
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
        const CCTK_REAL Gamma = get_Gamma_eff(rho_b[index], pressure[index]);

        ghl_ppm_no_rho_P(
              press_stencil, var_data,
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
}
