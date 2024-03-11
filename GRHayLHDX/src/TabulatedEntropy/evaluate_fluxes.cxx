#include "GRHayLHDX.h"

template <int flux_dir>
void GRHayLHDX_tabulated_entropy_evaluate_fluxes_dir(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_tabulated_entropy_evaluate_fluxes;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D2layout ccc_layout(cctkGH, {1, 1, 1});

  constexpr std::array<int, Loop::dim> facetype = {flux_dir!=0, flux_dir!=1, flux_dir!=2};
  const Loop::GF3D2layout flux_layout(cctkGH, facetype);

  constexpr void (*calculate_characteristic_speed)(
        ghl_primitive_quantities *restrict prims_r,
        ghl_primitive_quantities *restrict prims_l,
        const ghl_eos_parameters *restrict eos,
        const ghl_metric_quantities *restrict ADM_metric_face,
        CCTK_REAL *cmin, CCTK_REAL *cmax)
    = flux_dir==0 ? &ghl_calculate_characteristic_speed_dirn0 :
      flux_dir==1 ? &ghl_calculate_characteristic_speed_dirn1 :
                    &ghl_calculate_characteristic_speed_dirn2 ;

  constexpr void (*calculate_HLLE_fluxes)(
        ghl_primitive_quantities *restrict prims_r,
        ghl_primitive_quantities *restrict prims_l,
        const ghl_eos_parameters *restrict eos,
        const ghl_metric_quantities *restrict ADM_metric_face,
        const CCTK_REAL cmin,
        const CCTK_REAL cmax,
        ghl_conservative_quantities *restrict cons_fluxes)
    = flux_dir==0 ? &ghl_calculate_HLLE_fluxes_dirn0_tabulated_entropy :
      flux_dir==1 ? &ghl_calculate_HLLE_fluxes_dirn1_tabulated_entropy :
                    &ghl_calculate_HLLE_fluxes_dirn2_tabulated_entropy ;

  Loop::GF3D2<const CCTK_REAL> v_flux_dir = flux_dir==0 ? vx :
                                            flux_dir==1 ? vy : vz;

  Loop::GF3D2<CCTK_REAL> rho_star_flux = flux_dir==0 ? rho_star_flux_x :
                                         flux_dir==1 ? rho_star_flux_y : rho_star_flux_z;

  Loop::GF3D2<CCTK_REAL> tau_flux = flux_dir==0 ? tau_flux_x :
                                    flux_dir==1 ? tau_flux_y : tau_flux_z;

  Loop::GF3D2<CCTK_REAL> Sx_flux = flux_dir==0 ? Sx_flux_x :
                                   flux_dir==1 ? Sx_flux_y : Sx_flux_z;

  Loop::GF3D2<CCTK_REAL> Sy_flux = flux_dir==0 ? Sy_flux_x :
                                   flux_dir==1 ? Sy_flux_y : Sy_flux_z;

  Loop::GF3D2<CCTK_REAL> Sz_flux = flux_dir==0 ? Sz_flux_x :
                                   flux_dir==1 ? Sz_flux_y : Sz_flux_z;

  Loop::GF3D2<CCTK_REAL> ent_flux = flux_dir==0 ? ent_flux_x :
                                    flux_dir==1 ? ent_flux_y : ent_flux_z;

  Loop::GF3D2<CCTK_REAL> Ye_flux = flux_dir==0 ? Ye_flux_x :
                                   flux_dir==1 ? Ye_flux_y : Ye_flux_z;

  grid.loop_int<facetype[0], facetype[1], facetype[2]>(
      grid.nghostzones,
      [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index indm2(ccc_layout, p.I - 2*p.DI[flux_dir]);
    const Loop::GF3D2index indm1(ccc_layout, p.I - p.DI[flux_dir]);
    const Loop::GF3D2index index(ccc_layout, p.I);
    const Loop::GF3D2index indp1(ccc_layout, p.I + p.DI[flux_dir]);

    const Loop::GF3D2index ind_flux(flux_layout, p.I);

    ghl_metric_quantities ADM_metric_face;
    GRHayLHDX_interpolate_metric_to_face(
          indm2, indm1, index, indp1,
          ccc_lapse, ccc_betax, ccc_betay, ccc_betaz,
          ccc_gxx, ccc_gxy, ccc_gxz,
          ccc_gyy, ccc_gyz, ccc_gzz,
          &ADM_metric_face);

    CCTK_REAL rho_stencil[6], press_stencil[6], v_flux[6];
    CCTK_REAL vx_stencil[6], vy_stencil[6], vz_stencil[6];
    CCTK_REAL ent_stencil[6], Ye_stencil[6];
    ghl_primitive_quantities prims_r, prims_l;

    for(int ind=0; ind<6; ind++) {
      // Stencil from -3 to +2 reconstructs to e.g. i-1/2
      const Loop::GF3D2index stencil(ccc_layout, p.I + (ind-3)*p.DI[flux_dir]);

      v_flux[ind]        = v_flux_dir(stencil);
      rho_stencil[ind]   = rho(stencil);
      press_stencil[ind] = press(stencil);
      vx_stencil[ind]    = vx(stencil);
      vy_stencil[ind]    = vy(stencil);
      vz_stencil[ind]    = vz(stencil);
      ent_stencil[ind]   = entropy(stencil);
      Ye_stencil[ind]    = Ye(stencil);
    }

    CCTK_REAL ftilde[2];
    ghl_compute_ftilde(ghl_params, press_stencil, v_flux, ftilde);

    ghl_ppm_reconstruction_with_steepening(ghl_params, press_stencil, 1.0, ftilde, rho_stencil, &prims_r.rho, &prims_l.rho);

    ghl_ppm_reconstruction(ftilde, press_stencil, &prims_r.press, &prims_l.press);
    ghl_ppm_reconstruction(ftilde, vx_stencil, &prims_r.vU[0], &prims_l.vU[0]);
    ghl_ppm_reconstruction(ftilde, vy_stencil, &prims_r.vU[1], &prims_l.vU[1]);
    ghl_ppm_reconstruction(ftilde, vz_stencil, &prims_r.vU[2], &prims_l.vU[2]);
    ghl_ppm_reconstruction(ftilde, ent_stencil, &prims_r.entropy, &prims_l.entropy);
    ghl_ppm_reconstruction(ftilde, Ye_stencil, &prims_r.Y_e, &prims_l.Y_e);

    prims_r.BU[0] = prims_r.BU[1] = prims_r.BU[2] = 0.0;
    prims_l.BU[0] = prims_l.BU[1] = prims_l.BU[2] = 0.0;

    prims_r.temperature = prims_l.temperature = temperature(index);

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

    rho_star_flux(ind_flux) = cons_fluxes.rho;
    tau_flux(ind_flux)      = cons_fluxes.tau;
    Sx_flux(ind_flux)       = cons_fluxes.SD[0];
    Sy_flux(ind_flux)       = cons_fluxes.SD[1];
    Sz_flux(ind_flux)       = cons_fluxes.SD[2];
    ent_flux(ind_flux)      = cons_fluxes.entropy;
    Ye_flux (ind_flux)      = cons_fluxes.Y_e;
  }); // staggered loop interior (e.g. flux_dir=0 gives vcc)
}

extern "C" void GRHayLHDX_tabulated_entropy_evaluate_fluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_tabulated_entropy_evaluate_fluxes;
  DECLARE_CCTK_PARAMETERS;

  GRHayLHDX_tabulated_entropy_evaluate_fluxes_dir<0>(cctkGH);
  GRHayLHDX_tabulated_entropy_evaluate_fluxes_dir<1>(cctkGH);
  GRHayLHDX_tabulated_entropy_evaluate_fluxes_dir<2>(cctkGH);
}
