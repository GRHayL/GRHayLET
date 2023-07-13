#include "GRHayLHDX.h"

/*
  This file gives an example of how we could get rid of individual
  flux variables. Note that we have to loop over a bunch of extra
  points. This is because the vvv-centered variables *have* to be
  completely filled in the interior for PreSync to be happy.
  It also requires multiple scheduled functions (see bottom of file)
  to repeatedly schedule this.
*/
template <int flux_dir>
void GRHayLHDX_evaluate_flux_dir(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_evaluate_flux_x;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL poison = 0.0/0.0;
  const Loop::GF3D2layout ccc_layout(cctkGH, {1, 1, 1});

  constexpr std::array<int, Loop::dim> facetype = {flux_dir!=0, flux_dir!=1, flux_dir!=2};
  const Loop::GF3D2layout flux_layout(cctkGH, {0, 0, 0});

  // These nested condtional ternary operators let us tell the compiler that
  // these pointers can be set at compiler time while making the templated functions
  // instead of using a switch statement at runtime.
  constexpr void (*calculate_characteristic_speed)(const ghl_primitive_quantities *restrict prims_r,
                                         const ghl_primitive_quantities *restrict prims_l,
                                         const ghl_eos_parameters *restrict eos,
                                         const ghl_metric_quantities *restrict ADM_metric_face,
                                         CCTK_REAL *cmin, CCTK_REAL *cmax)
    = flux_dir==0 ? &ghl_calculate_characteristic_speed_dirn0 :
      flux_dir==1 ? &ghl_calculate_characteristic_speed_dirn1 :
                    &ghl_calculate_characteristic_speed_dirn2 ;

  constexpr void (*calculate_HLLE_fluxes)(const ghl_primitive_quantities *restrict prims_r,
                                const ghl_primitive_quantities *restrict prims_l,
                                const ghl_eos_parameters *restrict eos,
                                const ghl_metric_quantities *restrict ADM_metric_face,
                                const CCTK_REAL cmin,
                                const CCTK_REAL cmax,
                                ghl_conservative_quantities *restrict cons_fluxes)
    = flux_dir==0 ? &ghl_calculate_HLLE_fluxes_dirn0 :
      flux_dir==1 ? &ghl_calculate_HLLE_fluxes_dirn1 :
                    &ghl_calculate_HLLE_fluxes_dirn2 ;

  Loop::GF3D2<const CCTK_REAL> v_flux_dir = flux_dir==0 ? vx :
                                            flux_dir==1 ? vy : vz;

  //Loop::GF3D2<CCTK_REAL> rho_star_flux = flux_dir==0 ? rho_star_flux_x :
  //                                       flux_dir==1 ? rho_star_flux_y : rho_star_flux_z;

  //Loop::GF3D2<CCTK_REAL> tau_flux = flux_dir==0 ? tau_flux_x :
  //                                  flux_dir==1 ? tau_flux_y : tau_flux_z;

  //Loop::GF3D2<CCTK_REAL> Sx_flux = flux_dir==0 ? Sx_flux_x :
  //                                 flux_dir==1 ? Sx_flux_y : Sx_flux_z;

  //Loop::GF3D2<CCTK_REAL> Sy_flux = flux_dir==0 ? Sy_flux_x :
  //                                 flux_dir==1 ? Sy_flux_y : Sy_flux_z;

  //Loop::GF3D2<CCTK_REAL> Sz_flux = flux_dir==0 ? Sz_flux_x :
  //                                 flux_dir==1 ? Sz_flux_y : Sz_flux_z;

  //grid.loop_int_device<facetype[0], facetype[1], facetype[2]>(
  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index indm2(ccc_layout, p.I - 2*p.DI[flux_dir]);
    const Loop::GF3D2index indm1(ccc_layout, p.I - p.DI[flux_dir]);
    const Loop::GF3D2index index(ccc_layout, p.I);
    const Loop::GF3D2index indp1(ccc_layout, p.I + p.DI[flux_dir]);

    const Loop::GF3D2index ind_flux(flux_layout, p.I);

    CCTK_REAL rho_stencil[6], press_stencil[6], v_flux[6];
    CCTK_REAL rhor, rhol, pressr, pressl;
    CCTK_REAL vel_stencil[3][6], vel_r[3], vel_l[3];

    for(int ind=0; ind<6; ind++) {
      // Stencil from -3 to +2 reconstructs to e.g. i-1/2
      const Loop::GF3D2index stencil(ccc_layout, p.I + (ind-3)*p.DI[flux_dir]);
      v_flux[ind] = v_flux_dir(stencil); // Could be smaller; doesn't use full stencil
      rho_stencil[ind] = rho_b(stencil);
      press_stencil[ind] = pressure(stencil);
      vel_stencil[0][ind] = vx(stencil);
      vel_stencil[1][ind] = vy(stencil);
      vel_stencil[2][ind] = vz(stencil);
    }

    // Compute Gamma_eff
    CCTK_REAL K, Gamma;
    ghl_hybrid_get_K_and_Gamma(ghl_eos, rho_b(index), &K, &Gamma);
    const CCTK_REAL P_cold = K*pow(rho_b(index), Gamma);
    const CCTK_REAL Gamma_eff = ghl_eos->Gamma_th + (Gamma - ghl_eos->Gamma_th)*P_cold/pressure(index);

    ghl_ppm(
          rho_stencil, press_stencil, vel_stencil,
          3, v_flux, Gamma_eff,
          &rhor, &rhol, &pressr, &pressl, vel_r, vel_l);

    ghl_metric_quantities ADM_metric_face;
    GRHayLHDX_interpolate_metric_to_face(
          indm2, indm1, index, indp1,
          ccc_lapse, ccc_betax, ccc_betay, ccc_betaz,
          ccc_gxx, ccc_gxy, ccc_gxz,
          ccc_gyy, ccc_gyz, ccc_gzz,
          &ADM_metric_face);

    ghl_primitive_quantities prims_r, prims_l;
    ghl_initialize_primitives(
          rhor, pressr, poison,
          vel_r[0], vel_r[1], vel_r[2],
          0.0, 0.0, 0.0,
          poison, poison, poison, // entropy, Y_e, temp
          &prims_r);

    ghl_initialize_primitives(
          rhol, pressl, poison,
          vel_l[0], vel_l[1], vel_l[2],
          0.0, 0.0, 0.0,
          poison, poison, poison, // entropy, Y_e, temp
          &prims_l);

    int speed_limited = 0;
    ghl_limit_v_and_compute_u0(
          ghl_eos, &ADM_metric_face, &prims_r, &speed_limited);
    ghl_limit_v_and_compute_u0(
          ghl_eos, &ADM_metric_face, &prims_l, &speed_limited);

    CCTK_REAL cmin, cmax;
    ghl_conservative_quantities cons_fluxes;
    calculate_characteristic_speed(&prims_r, &prims_l, ghl_eos, &ADM_metric_face, &cmin, &cmax);
    calculate_HLLE_fluxes(&prims_r, &prims_l, ghl_eos, &ADM_metric_face, cmin, cmax, &cons_fluxes);

    rho_star_flux(ind_flux) = cons_fluxes.rho;
    tau_flux(ind_flux) = cons_fluxes.tau;
    Sx_flux(ind_flux)  = cons_fluxes.SD[0];
    Sy_flux(ind_flux)  = cons_fluxes.SD[1];
    Sz_flux(ind_flux)  = cons_fluxes.SD[2];
  }); // staggered loop interior (e.g. flux_dir=0 gives vcc)
}

extern "C" void GRHayLHDX_evaluate_flux_x(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_evaluate_flux_x;
  DECLARE_CCTK_PARAMETERS;

  GRHayLHDX_evaluate_flux_dir<0>(cctkGH);
}

extern "C" void GRHayLHDX_evaluate_flux_y(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_evaluate_flux_y;
  DECLARE_CCTK_PARAMETERS;

  GRHayLHDX_evaluate_flux_dir<1>(cctkGH);
}

extern "C" void GRHayLHDX_evaluate_flux_z(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_evaluate_flux_z;
  DECLARE_CCTK_PARAMETERS;

  GRHayLHDX_evaluate_flux_dir<2>(cctkGH);
}
