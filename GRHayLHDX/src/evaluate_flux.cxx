#include "GRHayLHDX.h"

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

template <int flux_dir, int tabulated_eos, int evolve_entropy>
void GRHayLHDX_evaluate_flux_dir(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_evaluate_flux;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL poison = 0.0/0.0;
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

  // I used constexpr with ternary operators elsewhere, but here there would be way too
  // many operators here, so I resort to if statements.
  void (*calculate_HLLE_fluxes)(
        ghl_primitive_quantities *restrict prims_r,
        ghl_primitive_quantities *restrict prims_l,
        const ghl_eos_parameters *restrict eos,
        const ghl_metric_quantities *restrict ADM_metric_face,
        const CCTK_REAL cmin,
        const CCTK_REAL cmax,
        ghl_conservative_quantities *restrict cons_fluxes);

  if(tabulated_eos) {
    if(evolve_entropy) {
      calculate_HLLE_fluxes = flux_dir==0 ? &ghl_calculate_HLLE_fluxes_dirn0_tabulated_entropy :
                              flux_dir==1 ? &ghl_calculate_HLLE_fluxes_dirn1_tabulated_entropy :
                                            &ghl_calculate_HLLE_fluxes_dirn2_tabulated_entropy ;
    } else {
      calculate_HLLE_fluxes = flux_dir==0 ? &ghl_calculate_HLLE_fluxes_dirn0_tabulated :
                              flux_dir==1 ? &ghl_calculate_HLLE_fluxes_dirn1_tabulated :
                                            &ghl_calculate_HLLE_fluxes_dirn2_tabulated ;
    }
  } else {
    if(evolve_entropy) {
      calculate_HLLE_fluxes = flux_dir==0 ? &ghl_calculate_HLLE_fluxes_dirn0_hybrid_entropy :
                              flux_dir==1 ? &ghl_calculate_HLLE_fluxes_dirn1_hybrid_entropy :
                                            &ghl_calculate_HLLE_fluxes_dirn2_hybrid_entropy ;
    } else {
      calculate_HLLE_fluxes = flux_dir==0 ? &ghl_calculate_HLLE_fluxes_dirn0_hybrid :
                              flux_dir==1 ? &ghl_calculate_HLLE_fluxes_dirn1_hybrid :
                                            &ghl_calculate_HLLE_fluxes_dirn2_hybrid ;
    }
  }

  constexpr double (*get_Gamma_eff)(
        const double,
        const double)
    = tabulated_eos ? &get_Gamma_eff_tabulated : &get_Gamma_eff_hybrid;

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

  //Loop::GF3D2<CCTK_REAL> ent_flux = flux_dir==0 ? ent_flux_x :
  //                                  flux_dir==1 ? ent_flux_y : ent_flux_z;

  //Loop::GF3D2<CCTK_REAL> Ye_flux = flux_dir==0 ? Ye_flux_x :
  //                                 flux_dir==1 ? Ye_flux_y : Ye_flux_z;

  // Count number of additional reconstructed variables
  constexpr int num_others = 3 + evolve_entropy + tabulated_eos;

  // If using the entropy, it should be the first reconstructed variable
  // after the three velocities
  constexpr int ent_index = 3 + !evolve_entropy;

  // If not using the entropy, then Ye should be the first reconstructed
  // variable after the three velocities
  constexpr int Ye_index = 4 - !evolve_entropy;

  grid.loop_int_device<facetype[0], facetype[1], facetype[2]>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index indm2(ccc_layout, p.I - 2*p.DI[flux_dir]);
    const Loop::GF3D2index indm1(ccc_layout, p.I - p.DI[flux_dir]);
    const Loop::GF3D2index index(ccc_layout, p.I);
    const Loop::GF3D2index indp1(ccc_layout, p.I + p.DI[flux_dir]);

    const Loop::GF3D2index ind_flux(flux_layout, p.I);

    CCTK_REAL rho_stencil[6], press_stencil[6], v_flux[6];
    CCTK_REAL rhor, rhol, pressr, pressl;
    CCTK_REAL others_stencil[num_others][6], others_r[num_others], others_l[num_others];

    for(int ind=0; ind<6; ind++) {
      // Stencil from -3 to +2 reconstructs to e.g. i-1/2
      const Loop::GF3D2index stencil(ccc_layout, p.I + (ind-3)*p.DI[flux_dir]);
      v_flux[ind]                    = v_flux_dir(stencil); // Could be smaller; doesn't use full stencil
      rho_stencil[ind]               = rho(stencil);
      press_stencil[ind]             = press(stencil);
      others_stencil[0][ind]         = vx(stencil);
      others_stencil[1][ind]         = vy(stencil);
      others_stencil[2][ind]         = vz(stencil);
      //others_stencil[ent_index][ind] = entropy(stencil);
      //others_stencil[Ye_index][ind] = Ye(stencil);
    }

    // Compute Gamma
    const double Gamma = get_Gamma_eff(rho(index), press(index));

    ghl_ppm(
          ghl_params, rho_stencil,
          press_stencil, others_stencil,
          3, v_flux, Gamma,
          &rhor, &rhol,
          &pressr, &pressl,
          others_r, others_l);

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
          others_r[0], others_r[1], others_r[2],
          0.0, 0.0, 0.0,
          0.0, 0.0, 0.0,
          //others_r[ent_index], others_r[Ye_index], poison,
          &prims_r);

    ghl_initialize_primitives(
          rhol, pressl, poison,
          others_l[0], others_l[1], others_l[2],
          0.0, 0.0, 0.0,
          0.0, 0.0, 0.0,
          //others_l[ent_index], others_l[Ye_index], poison,
          &prims_l);

    int speed_limited CCTK_ATTRIBUTE_UNUSED = ghl_limit_v_and_compute_u0(ghl_params, &ADM_metric_face, &prims_r);
    speed_limited = ghl_limit_v_and_compute_u0(ghl_params, &ADM_metric_face, &prims_l);

    CCTK_REAL cmin, cmax;
    ghl_conservative_quantities cons_fluxes;
    calculate_characteristic_speed(&prims_r, &prims_l, ghl_eos, &ADM_metric_face, &cmin, &cmax);
    calculate_HLLE_fluxes(&prims_r, &prims_l, ghl_eos, &ADM_metric_face, cmin, cmax, &cons_fluxes);

    rho_star_flux(ind_flux) = cons_fluxes.rho;
    tau_flux(ind_flux) = cons_fluxes.tau;
    Sx_flux(ind_flux)  = cons_fluxes.SD[0];
    Sy_flux(ind_flux)  = cons_fluxes.SD[1];
    Sz_flux(ind_flux)  = cons_fluxes.SD[2];
    //ent_flux(ind_flux) = cons_fluxes.entropy;
    //Ye_flux (ind_flux) = cons_fluxes.Y_e;
  }); // staggered loop interior (e.g. flux_dir=0 gives vcc)
}

extern "C" void GRHayLHDX_evaluate_flux(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_evaluate_flux;
  DECLARE_CCTK_PARAMETERS;

  if(ghl_eos->eos_type==ghl_eos_hybrid) {
    if(ghl_params->evolve_entropy) {
      GRHayLHDX_evaluate_flux_dir<0, 0, 1>(cctkGH);
      GRHayLHDX_evaluate_flux_dir<1, 0, 1>(cctkGH);
      GRHayLHDX_evaluate_flux_dir<2, 0, 1>(cctkGH);
    } else {
      GRHayLHDX_evaluate_flux_dir<0, 0, 0>(cctkGH);
      GRHayLHDX_evaluate_flux_dir<1, 0, 0>(cctkGH);
      GRHayLHDX_evaluate_flux_dir<2, 0, 0>(cctkGH);
    }
  } else if(ghl_eos->eos_type==ghl_eos_tabulated) {
    if(ghl_params->evolve_entropy) {
      GRHayLHDX_evaluate_flux_dir<0, 1, 1>(cctkGH);
      GRHayLHDX_evaluate_flux_dir<1, 1, 1>(cctkGH);
      GRHayLHDX_evaluate_flux_dir<2, 1, 1>(cctkGH);
    } else {
      GRHayLHDX_evaluate_flux_dir<0, 1, 0>(cctkGH);
      GRHayLHDX_evaluate_flux_dir<1, 1, 0>(cctkGH);
      GRHayLHDX_evaluate_flux_dir<2, 1, 0>(cctkGH);
    }
  }
}
