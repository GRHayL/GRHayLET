#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "GRHayLib.h"

// The inner two points of the interpolation function use
// the value of A_in, and the outer two points use A_out.
#define A_out -0.0625
#define A_in  0.5625
// Interpolates to the +1/2 face of point Var
#define COMPUTE_FCVAL(Varm1,Var,Varp1,Varp2) (A_out*(Varm1) + A_in*(Var) + A_in*(Varp1) + A_out*(Varp2))


static void GRHayLHD_interpolate_metric_to_face(
      const cGH *cctkGH,
      const int i, const int j, const int k,
      const int flux_dir,
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
      ghl_metric_quantities *restrict metric) {

  const int xdir = (flux_dir == 0);
  const int ydir = (flux_dir == 1);
  const int zdir = (flux_dir == 2);

  const int indm2  = CCTK_GFINDEX3D(cctkGH, i-2*xdir, j-2*ydir, k-2*zdir);
  const int indm1  = CCTK_GFINDEX3D(cctkGH, i-xdir,   j-ydir,   k-zdir);
  const int index  = CCTK_GFINDEX3D(cctkGH, i,        j ,       k);
  const int indp1  = CCTK_GFINDEX3D(cctkGH, i+xdir,   j+ydir,   k+zdir);

  const CCTK_REAL face_lapse = COMPUTE_FCVAL(lapse[indm2], lapse[indm1], lapse[index], lapse[indp1]);
  const CCTK_REAL face_betax = COMPUTE_FCVAL(betax[indm2], betax[indm1], betax[index], betax[indp1]);
  const CCTK_REAL face_betay = COMPUTE_FCVAL(betay[indm2], betay[indm1], betay[index], betay[indp1]);
  const CCTK_REAL face_betaz = COMPUTE_FCVAL(betaz[indm2], betaz[indm1], betaz[index], betaz[indp1]);

  const CCTK_REAL face_gxx = COMPUTE_FCVAL(gxx[indm2], gxx[indm1], gxx[index], gxx[indp1]);
  const CCTK_REAL face_gxy = COMPUTE_FCVAL(gxy[indm2], gxy[indm1], gxy[index], gxy[indp1]);
  const CCTK_REAL face_gxz = COMPUTE_FCVAL(gxz[indm2], gxz[indm1], gxz[index], gxz[indp1]);
  const CCTK_REAL face_gyy = COMPUTE_FCVAL(gyy[indm2], gyy[indm1], gyy[index], gyy[indp1]);
  const CCTK_REAL face_gyz = COMPUTE_FCVAL(gyz[indm2], gyz[indm1], gyz[index], gyz[indp1]);
  const CCTK_REAL face_gzz = COMPUTE_FCVAL(gzz[indm2], gzz[indm1], gzz[index], gzz[indp1]);

  ghl_initialize_metric(
        face_lapse,
        face_betax, face_betay, face_betaz,
        face_gxx, face_gxy, face_gxz,
        face_gyy, face_gyz, face_gzz,
        metric);
}


void GRHayLHD_tabulated_evaluate_fluxes_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
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
        calculate_HLLE_fluxes = ghl_calculate_HLLE_fluxes_dirn0_tabulated_entropy;
        break;
      case 1:
        v_flux_dir = vy;
        calculate_characteristic_speed = ghl_calculate_characteristic_speed_dirn1;
        calculate_HLLE_fluxes = ghl_calculate_HLLE_fluxes_dirn1_tabulated_entropy;
        break;
      case 2:
        v_flux_dir = vz;
        calculate_characteristic_speed = ghl_calculate_characteristic_speed_dirn2;
        calculate_HLLE_fluxes = ghl_calculate_HLLE_fluxes_dirn2_tabulated_entropy;
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
          CCTK_REAL Ye_stencil[6];
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
            Ye_stencil[ind]    = Y_e[stencil];
          }

          CCTK_REAL ftilde[2];
          ghl_compute_ftilde(ghl_params, press_stencil, v_flux, ftilde);

          const CCTK_REAL Gamma = 1.0;
          ghl_ppm_reconstruction_with_steepening(ghl_params, press_stencil, Gamma, ftilde, rho_stencil, &prims_r.rho, &prims_l.rho);

          ghl_ppm_reconstruction(ftilde, press_stencil, &prims_r.press, &prims_l.press);
          ghl_ppm_reconstruction(ftilde, vx_stencil, &prims_r.vU[0], &prims_l.vU[0]);
          ghl_ppm_reconstruction(ftilde, vy_stencil, &prims_r.vU[1], &prims_l.vU[1]);
          ghl_ppm_reconstruction(ftilde, vz_stencil, &prims_r.vU[2], &prims_l.vU[2]);
          ghl_ppm_reconstruction(ftilde, Ye_stencil, &prims_r.Y_e, &prims_l.Y_e);

          prims_r.BU[0] = prims_r.BU[1] = prims_r.BU[2] = 0.0;
          prims_l.BU[0] = prims_l.BU[1] = prims_l.BU[2] = 0.0;

          prims_r.temperature = prims_l.temperature = temperature[index];

          int speed_limited CCTK_ATTRIBUTE_UNUSED = ghl_limit_v_and_compute_u0(ghl_params, &ADM_metric_face, &prims_r);
          speed_limited = ghl_limit_v_and_compute_u0(ghl_params, &ADM_metric_face, &prims_l);

          // We must now compute eps and T
          ghl_tabulated_enforce_bounds_rho_Ye_P(ghl_eos, &prims_r.rho, &prims_r.Y_e, &prims_r.press);
          ghl_tabulated_compute_eps_T_from_P(ghl_eos, prims_r.rho, prims_r.Y_e, prims_r.press,
                                             &prims_r.eps, &prims_r.temperature);
  
          ghl_tabulated_enforce_bounds_rho_Ye_P(ghl_eos, &prims_l.rho, &prims_l.Y_e, &prims_l.press);
          ghl_tabulated_compute_eps_T_from_P(ghl_eos, prims_l.rho, prims_l.Y_e, prims_l.press,
                                             &prims_l.eps, &prims_l.temperature);

          CCTK_REAL cmin, cmax;
          ghl_conservative_quantities cons_fluxes;
          calculate_characteristic_speed(&prims_r, &prims_l, ghl_eos, &ADM_metric_face, &cmin, &cmax);
          calculate_HLLE_fluxes(&prims_r, &prims_l, ghl_eos, &ADM_metric_face, cmin, cmax, &cons_fluxes);

          rho_star_flux[index] = cons_fluxes.rho;
          tau_flux     [index] = cons_fluxes.tau;
          st_x_flux [index] = cons_fluxes.SD[0];
          st_y_flux [index] = cons_fluxes.SD[1];
          st_z_flux [index] = cons_fluxes.SD[2];
          Ye_star_flux [index] = cons_fluxes.Y_e;
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
          st_x_rhs [index] += dxi*(st_x_flux [index] - st_x_flux [indp1]);
          st_y_rhs [index] += dxi*(st_y_flux [index] - st_y_flux [indp1]);
          st_z_rhs [index] += dxi*(st_z_flux [index] - st_z_flux [indp1]);
          Ye_star_rhs [index] += dxi*(Ye_star_flux [index] - Ye_star_flux [indp1]);
        }
      }
    }
  }
}


