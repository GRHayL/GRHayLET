#include "GRHayLMHDX.hxx"

template <int flux_dir>
void GRHayLMHDX_evaluate_flux_dir(
      CCTK_ARGUMENTS,
      const Loop::GF3D5layout local_layout,
      Loop::GF3D5<CCTK_REAL> v1_r,
      Loop::GF3D5<CCTK_REAL> v1_l,
      Loop::GF3D5<CCTK_REAL> v2_r,
      Loop::GF3D5<CCTK_REAL> v2_l,
      Loop::GF3D5<CCTK_REAL> cmin,
      Loop::GF3D5<CCTK_REAL> cmax,
      Loop::GF3D5<CCTK_REAL> rho_flux,
      Loop::GF3D5<CCTK_REAL> tau_flux,
      Loop::GF3D5<CCTK_REAL> S_x_flux,
      Loop::GF3D5<CCTK_REAL> S_y_flux,
      Loop::GF3D5<CCTK_REAL> S_z_flux) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLMHDX_evaluate_MHD_rhs;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL poison = 0.0/0.0;
  const Loop::GF3D2layout ccc_layout(cctkGH, {1, 1, 1});

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

  // Here, we're setting up an index for Bvec to avoid reconstructing by using
  // B_stagger.
  constexpr int B0 = flux_dir==0 ? 0 :
                     flux_dir==1 ? 1 :
                                   2 ;

  constexpr int B1 = flux_dir==0 ? 1 :
                     flux_dir==1 ? 2 :
                                   0 ;

  constexpr int B2 = flux_dir==0 ? 2 :
                     flux_dir==1 ? 0 :
                                   1 ;

  Loop::GF3D2<const CCTK_REAL> Bvec[3] = {Bx_center, By_center, Bz_center};

  Loop::GF3D2<const CCTK_REAL> B_stagger = flux_dir==0 ? Bx_stagger :
                                           flux_dir==1 ? By_stagger : Bz_stagger;

  Loop::GF3D2<const CCTK_REAL> v_flux_dir = flux_dir==0 ? vx :
                                            flux_dir==1 ? vy : vz;

  constexpr std::array<int, Loop::dim> facetype = {flux_dir!=0, flux_dir!=1, flux_dir!=2};

  Arith::vect<int, 3> bnd_min, bnd_max;
  grid.boundary_box<facetype[0], facetype[1], facetype[2]>(grid.nghostzones, bnd_min, bnd_max);
  Arith::vect<int, 3> imin, imax;
  grid.box_int<facetype[0], facetype[1], facetype[2]>(grid.nghostzones, imin, imax);

  // Extend the loop box. You must make sure that you extend by at most the
  // number of ghosts.
  for(int dir=0; dir<3; dir++) {
    if(dir==flux_dir) continue;
    imin[dir] -= 3;
    imax[dir] += 3;
  }
  //if(flux_dir != 0)
  //  imin[0] -= 3; imax[0] += 3;
  //if(flux_dir != 1)
  //  imin[1] -= 3; imax[1] += 3;
  //if(flux_dir != 2)
  //  imin[2] -= 3; imax[2] += 3;

  grid.loop_box_device<facetype[0], facetype[1], facetype[2]>(
      bnd_min, bnd_max, imin, imax,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

    const Loop::GF3D2index indm2(ccc_layout, p.I - 2*p.DI[flux_dir]);
    const Loop::GF3D2index indm1(ccc_layout, p.I - p.DI[flux_dir]);
    const Loop::GF3D2index index(ccc_layout, p.I);
    const Loop::GF3D2index indp1(ccc_layout, p.I + p.DI[flux_dir]);

    CCTK_REAL rho_stencil[6], press_stencil[6], v_flux[6];
    CCTK_REAL rhor, rhol, pressr, pressl, B_r[3], B_l[3];
    CCTK_REAL vars_stencil[5][6], vars_r[5], vars_l[5];

    for(int ind=0; ind<6; ind++) {
      // Stencil from -3 to +2 reconstructs to e.g. i-1/2
      const Loop::GF3D2index stencil(ccc_layout, p.I + (ind-3)*p.DI[flux_dir]);
      v_flux[ind] = v_flux_dir(stencil); // Could be smaller; doesn't use full stencil
      rho_stencil[ind] = rho_b(stencil);
      press_stencil[ind] = pressure(stencil);
      vars_stencil[0][ind] = vx(stencil);
      vars_stencil[1][ind] = vy(stencil);
      vars_stencil[2][ind] = vz(stencil);
      vars_stencil[3][ind] = Bvec[B1](stencil);
      vars_stencil[4][ind] = Bvec[B2](stencil);
    }

    // Compute Gamma_eff
    CCTK_REAL K, Gamma;
    ghl_hybrid_get_K_and_Gamma(ghl_eos, rho_b(index), &K, &Gamma);
    const CCTK_REAL P_cold = K*pow(rho_b(index), Gamma);
    const CCTK_REAL Gamma_eff = ghl_eos->Gamma_th + (Gamma - ghl_eos->Gamma_th)*P_cold/pressure(index);

    ghl_ppm(
          rho_stencil, press_stencil, vars_stencil,
          5, v_flux, Gamma_eff,
          &rhor, &rhol, &pressr, &pressl, vars_r, vars_l);

    ghl_metric_quantities ADM_metric_face;
    GRHayLMHDX_interpolate_metric_to_face(
          indm2, indm1, index, indp1,
          ccc_lapse, ccc_betax, ccc_betay, ccc_betaz,
          ccc_gxx, ccc_gxy, ccc_gxz,
          ccc_gyy, ccc_gyz, ccc_gzz,
          &ADM_metric_face);

    const Loop::GF3D5index ind_flux(local_layout, p.I);
    // These if's should vanish when the compile
    // makes the actual templated functions.
    if(flux_dir==0) {
      v1_r(ind_flux) = vars_r[0];
      v1_l(ind_flux) = vars_l[0];
      v2_r(ind_flux) = vars_r[1];
      v2_l(ind_flux) = vars_l[1];
    }
    if(flux_dir==1) {
      v1_r(ind_flux) = vars_r[1];
      v1_l(ind_flux) = vars_l[1];
      v2_r(ind_flux) = vars_r[2];
      v2_l(ind_flux) = vars_l[2];
    }
    if(flux_dir==2) {
      v1_r(ind_flux) = vars_r[2];
      v1_l(ind_flux) = vars_l[2];
      v2_r(ind_flux) = vars_r[0];
      v2_l(ind_flux) = vars_l[0];
    }

    B_r[B1] = vars_r[3];
    B_r[B2] = vars_r[4];

    B_l[B1] = vars_l[3];
    B_l[B2] = vars_l[4];

    // B_stagger is densitized, but B_center is not.
    B_r[B0] = B_l[B0] = B_stagger(p.I)/ADM_metric_face.sqrt_detgamma;

    ghl_primitive_quantities prims_r, prims_l;
    ghl_initialize_primitives(
          rhor, pressr, poison,
          vars_r[0], vars_r[1], vars_r[2],
          B_r[0], B_r[1], B_r[2],
          poison, poison, poison, // entropy, Y_e, temp
          &prims_r);

    ghl_initialize_primitives(
          rhol, pressl, poison,
          vars_l[0], vars_l[1], vars_l[2],
          B_l[0], B_l[1], B_l[2],
          poison, poison, poison, // entropy, Y_e, temp
          &prims_l);

    int speed_limited = 0;
    ghl_limit_v_and_compute_u0(
          ghl_eos, &ADM_metric_face, &prims_r, &speed_limited);
    ghl_limit_v_and_compute_u0(
          ghl_eos, &ADM_metric_face, &prims_l, &speed_limited);

    ghl_conservative_quantities cons_fluxes;
    calculate_characteristic_speed(&prims_r, &prims_l, ghl_eos, &ADM_metric_face, &cmin(ind_flux), &cmax(ind_flux));
    calculate_HLLE_fluxes(&prims_r, &prims_l, ghl_eos, &ADM_metric_face, cmin(ind_flux), cmax(ind_flux), &cons_fluxes);

    rho_flux(ind_flux) = cons_fluxes.rho;
    tau_flux(ind_flux) = cons_fluxes.tau;
    S_x_flux(ind_flux) = cons_fluxes.SD[0];
    S_y_flux(ind_flux) = cons_fluxes.SD[1];
    S_z_flux(ind_flux) = cons_fluxes.SD[2];
  }); // staggered loop box
}

template <int flux_dir>
void GRHayLMHDX_evaluate_flux_source_rhs_dir(
      CCTK_ARGUMENTS,
      const Loop::GF3D5layout local_layout,
      const Loop::GF3D5<CCTK_REAL> rho_flux,
      const Loop::GF3D5<CCTK_REAL> tau_flux,
      const Loop::GF3D5<CCTK_REAL> S_x_flux,
      const Loop::GF3D5<CCTK_REAL> S_y_flux,
      const Loop::GF3D5<CCTK_REAL> S_z_flux) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLMHDX_evaluate_MHD_rhs;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL poison = 0.0/0.0;
  const Loop::GF3D2layout ccc_layout(cctkGH, {1, 1, 1});

  // These nested condtional ternary operators let us tell the compiler that
  // these pointers can be set at compiler time while making the templated functions
  // instead of using a switch statement at runtime.
  constexpr void (*calculate_source_terms)(const ghl_primitive_quantities *restrict prims,
                                 const ghl_eos_parameters *restrict eos,
                                 const ghl_metric_quantities *restrict ADM_metric,
                                 const ghl_metric_quantities *restrict metric_derivs,
                                 ghl_conservative_quantities *restrict cons_sources)
    = flux_dir==0 ? &ghl_calculate_source_terms_dirn0 :
      flux_dir==1 ? &ghl_calculate_source_terms_dirn1 :
                    &ghl_calculate_source_terms_dirn2 ;

  const CCTK_REAL dxi = 1.0/CCTK_DELTA_SPACE(flux_dir);

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index indm2(ccc_layout, p.I - 2*p.DI[flux_dir]);
    const Loop::GF3D2index indm1(ccc_layout, p.I - p.DI[flux_dir]);
    const Loop::GF3D2index index(ccc_layout, p.I);
    const Loop::GF3D2index indp1(ccc_layout, p.I + p.DI[flux_dir]);
    const Loop::GF3D2index indp2(ccc_layout, p.I + 2*p.DI[flux_dir]);

    const Loop::GF3D5index ind_flux(local_layout, p.I);
    const Loop::GF3D5index ind_flp1(local_layout, p.I + p.DI[flux_dir]);

    rho_star_rhs(index) += dxi*(rho_flux(ind_flux) - rho_flux(ind_flp1));
    tau_rhs(index)      += dxi*(tau_flux(ind_flux) - tau_flux(ind_flp1));
    Stildex_rhs(index)  += dxi*(S_x_flux(ind_flux) - S_x_flux(ind_flp1));
    Stildey_rhs(index)  += dxi*(S_y_flux(ind_flux) - S_y_flux(ind_flp1));
    Stildez_rhs(index)  += dxi*(S_z_flux(ind_flux) - S_z_flux(ind_flp1));

    ghl_metric_quantities ADM_metric;
    ghl_initialize_metric(ccc_lapse(index),
          ccc_betax(index), ccc_betay(index), ccc_betaz(index),
          ccc_gxx(index), ccc_gxy(index), ccc_gxz(index),
          ccc_gyy(index), ccc_gyz(index), ccc_gzz(index),
          &ADM_metric);

    ghl_primitive_quantities prims;
    ghl_initialize_primitives(
          rho_b(index), pressure(index), eps(index),
          vx(index), vy(index), vz(index),
          Bx_center(index), By_center(index), Bz_center(index),
          poison, poison, poison, // entropy, Y_e, temp
          &prims);

    int speed_limited = 0;
    ghl_limit_v_and_compute_u0(
          ghl_eos, &ADM_metric, &prims, &speed_limited);

    ghl_metric_quantities ADM_metric_derivs;

    GRHayLMHDX_compute_metric_derivs(
          dxi, indm2, indm1, indp1, indp2,
          ccc_lapse, ccc_betax, ccc_betay, ccc_betaz,
          ccc_gxx, ccc_gxy, ccc_gxz,
          ccc_gyy, ccc_gyz, ccc_gzz,
          &ADM_metric_derivs);

    ghl_conservative_quantities cons_source;
    cons_source.tau = 0.0;
    cons_source.SD[0] = 0.0;
    cons_source.SD[1] = 0.0;
    cons_source.SD[2] = 0.0;

    calculate_source_terms(&prims, ghl_eos, &ADM_metric, &ADM_metric_derivs, &cons_source);
    tau_rhs(index)     += cons_source.tau;
    Stildex_rhs(index) += cons_source.SD[0];
    Stildey_rhs(index) += cons_source.SD[1];
    Stildez_rhs(index) += cons_source.SD[2];
  }); // ccc loop interior
}

extern "C" void GRHayLMHDX_evaluate_MHD_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLMHDX_evaluate_MHD_rhs;
  DECLARE_CCTK_PARAMETERS;

  Arith::vect<int, Loop::dim> imin, imax;

  // Set up temporary variables with vcc layout including
  // y and z ghostzones
  grid.box_int<0, 1, 1>(grid.nghostzones, imin, imax);
  imin[1] -= 3; imax[1] += 3;
  imin[2] -= 3; imax[2] += 3;
  const Loop::GF3D5layout vcc_layout(imin, imax);
  Loop::GF3D5vector<CCTK_REAL> tmpsx(vcc_layout, 11);
  int itmp = 0; const auto make_gfx = [&]() { return Loop::GF3D5<CCTK_REAL>(tmpsx(itmp++)); };

  Loop::GF3D5<CCTK_REAL> vx_xr(make_gfx());
  Loop::GF3D5<CCTK_REAL> vx_xl(make_gfx());
  Loop::GF3D5<CCTK_REAL> vy_xr(make_gfx());
  Loop::GF3D5<CCTK_REAL> vy_xl(make_gfx());
  Loop::GF3D5<CCTK_REAL> cmin_x(make_gfx());
  Loop::GF3D5<CCTK_REAL> cmax_x(make_gfx());
  Loop::GF3D5<CCTK_REAL> rho_flux_x(make_gfx());
  Loop::GF3D5<CCTK_REAL> tau_flux_x(make_gfx());
  Loop::GF3D5<CCTK_REAL> S_x_flux_x(make_gfx());
  Loop::GF3D5<CCTK_REAL> S_y_flux_x(make_gfx());
  Loop::GF3D5<CCTK_REAL> S_z_flux_x(make_gfx());

  GRHayLMHDX_evaluate_flux_dir<0>(cctkGH, vcc_layout, vx_xr, vx_xl, vy_xr, vy_xl, cmin_x, cmax_x, rho_flux_x, tau_flux_x, S_x_flux_x, S_y_flux_x, S_z_flux_x);
  GRHayLMHDX_evaluate_flux_source_rhs_dir<0>(cctkGH, vcc_layout, rho_flux_x, tau_flux_x, S_x_flux_x, S_y_flux_x, S_z_flux_x);

  // Set up temporary variables with cvc layout including
  // x and z ghostzones
  grid.box_int<1, 0, 1>(grid.nghostzones, imin, imax);
  imin[0] -= 3; imax[0] += 3;
  imin[2] -= 3; imax[2] += 3;
  const Loop::GF3D5layout cvc_layout(imin, imax);
  Loop::GF3D5vector<CCTK_REAL> tmpsy(cvc_layout, 11);
  itmp = 0; const auto make_gfy = [&]() { return Loop::GF3D5<CCTK_REAL>(tmpsy(itmp++)); };

  Loop::GF3D5<CCTK_REAL> vy_yr(make_gfy());
  Loop::GF3D5<CCTK_REAL> vy_yl(make_gfy());
  Loop::GF3D5<CCTK_REAL> vz_yr(make_gfy());
  Loop::GF3D5<CCTK_REAL> vz_yl(make_gfy());
  Loop::GF3D5<CCTK_REAL> cmin_y(make_gfy());
  Loop::GF3D5<CCTK_REAL> cmax_y(make_gfy());
  Loop::GF3D5<CCTK_REAL> rho_flux_y(make_gfy());
  Loop::GF3D5<CCTK_REAL> tau_flux_y(make_gfy());
  Loop::GF3D5<CCTK_REAL> S_x_flux_y(make_gfy());
  Loop::GF3D5<CCTK_REAL> S_y_flux_y(make_gfy());
  Loop::GF3D5<CCTK_REAL> S_z_flux_y(make_gfy());

  GRHayLMHDX_evaluate_flux_dir<1>(cctkGH, cvc_layout, vy_yr, vy_yl, vz_yr, vz_yl, cmin_y, cmax_y, rho_flux_y, tau_flux_y, S_x_flux_y, S_y_flux_y, S_z_flux_y);
  GRHayLMHDX_evaluate_flux_source_rhs_dir<1>(cctkGH, cvc_layout, rho_flux_y, tau_flux_y, S_x_flux_y, S_y_flux_y, S_z_flux_y);

  // Set up temporary variables with ccv layout including
  // x and y ghostzones
  grid.box_int<1, 1, 0>(grid.nghostzones, imin, imax);
  imin[0] -= 3; imax[0] += 3;
  imin[1] -= 3; imax[1] += 3;
  const Loop::GF3D5layout ccv_layout(imin, imax);
  Loop::GF3D5vector<CCTK_REAL> tmpsz(ccv_layout, 11);
  itmp = 0; const auto make_gfz = [&]() { return Loop::GF3D5<CCTK_REAL>(tmpsz(itmp++)); };

  Loop::GF3D5<CCTK_REAL> vz_zr(make_gfz());
  Loop::GF3D5<CCTK_REAL> vz_zl(make_gfz());
  Loop::GF3D5<CCTK_REAL> vx_zr(make_gfz());
  Loop::GF3D5<CCTK_REAL> vx_zl(make_gfz());
  Loop::GF3D5<CCTK_REAL> cmin_z(make_gfz());
  Loop::GF3D5<CCTK_REAL> cmax_z(make_gfz());
  Loop::GF3D5<CCTK_REAL> rho_flux_z(make_gfz());
  Loop::GF3D5<CCTK_REAL> tau_flux_z(make_gfz());
  Loop::GF3D5<CCTK_REAL> S_x_flux_z(make_gfz());
  Loop::GF3D5<CCTK_REAL> S_y_flux_z(make_gfz());
  Loop::GF3D5<CCTK_REAL> S_z_flux_z(make_gfz());

  GRHayLMHDX_evaluate_flux_dir<2>(cctkGH, ccv_layout, vz_yr, vz_zl, vx_zr, vx_zl, cmin_z, cmax_z, rho_flux_z, tau_flux_z, S_x_flux_z, S_y_flux_z, S_z_flux_z);
  GRHayLMHDX_evaluate_flux_source_rhs_dir<2>(cctkGH, ccv_layout, rho_flux_z, tau_flux_z, S_x_flux_z, S_y_flux_z, S_z_flux_z);

  grid.loop_int_device<1, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    Ax_rhs(p.I) = 0.0;
  });

  grid.loop_int_device<0, 1, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    Ay_rhs(p.I) = 0.0;
  });

  grid.loop_int_device<0, 0, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    Az_rhs(p.I) = 0.0;
  });
}
