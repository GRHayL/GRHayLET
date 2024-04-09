#include "IllinoisGRMHDX.hxx"

static inline CCTK_REAL get_Gamma_eff(
      const CCTK_REAL rho_in,
      const CCTK_REAL press_in) {
  CCTK_REAL K, Gamma;
  ghl_hybrid_get_K_and_Gamma(ghl_eos, rho_in, &K, &Gamma);
  const CCTK_REAL P_cold = K*pow(rho_in, Gamma);
  return ghl_eos->Gamma_th + (Gamma - ghl_eos->Gamma_th)*P_cold/press_in;
}

template <int flux_dir>
void IllinoisGRMHDX_evaluate_flux_dir(
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
  DECLARE_CCTK_ARGUMENTSX_IllinoisGRMHDX_hybrid_evaluate_fluxes_rhs;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL poison = 0.0/0.0;
  const Loop::GF3D2layout ccc_layout(cctkGH, {1, 1, 1});

  // These nested condtional ternary operators let us tell the compiler that
  // these pointers can be set at compiler time while making the templated functions
  // instead of using a switch statement at runtime.
  constexpr void (*calculate_characteristic_speed)(
                                         ghl_primitive_quantities *restrict prims_r,
                                         ghl_primitive_quantities *restrict prims_l,
                                         const ghl_eos_parameters *restrict eos,
                                         const ghl_metric_quantities *restrict ADM_metric_face,
                                         CCTK_REAL *cmin, CCTK_REAL *cmax)
    = flux_dir==0 ? ghl_calculate_characteristic_speed_dirn0 :
      flux_dir==1 ? ghl_calculate_characteristic_speed_dirn1 :
                    ghl_calculate_characteristic_speed_dirn2 ;

  constexpr void (*calculate_HLLE_fluxes)(
                                ghl_primitive_quantities *restrict prims_r,
                                ghl_primitive_quantities *restrict prims_l,
                                const ghl_eos_parameters *restrict eos,
                                const ghl_metric_quantities *restrict ADM_metric_face,
                                const CCTK_REAL cmin,
                                const CCTK_REAL cmax,
                                ghl_conservative_quantities *restrict cons_fluxes)
    = flux_dir==0 ? ghl_calculate_HLLE_fluxes_dirn0_hybrid :
      flux_dir==1 ? ghl_calculate_HLLE_fluxes_dirn1_hybrid :
                    ghl_calculate_HLLE_fluxes_dirn2_hybrid ;

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

  Loop::GF3D2<const CCTK_REAL> Bvec[3] = {Bvecx, Bvecy, Bvecz};

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

    ghl_metric_quantities ADM_metric_face;
    IllinoisGRMHDX_interpolate_metric_to_face(
          indm2, indm1, index, indp1,
          ccc_lapse, ccc_betax, ccc_betay, ccc_betaz,
          ccc_gxx, ccc_gxy, ccc_gxz,
          ccc_gyy, ccc_gyz, ccc_gzz,
          &ADM_metric_face);

    CCTK_REAL rho_stencil[6], press_stencil[6], v_flux[6];
    CCTK_REAL B1_stencil[6], B2_stencil[6];
    CCTK_REAL vx_stencil[6], vy_stencil[6], vz_stencil[6];
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
      B1_stencil[ind]    = Bvec[B1](stencil);
      B2_stencil[ind]    = Bvec[B2](stencil);
    }

    // Compute Gamma_eff
    CCTK_REAL ftilde[2];
    ghl_compute_ftilde(ghl_params, press_stencil, v_flux, ftilde);

    const CCTK_REAL Gamma = get_Gamma_eff(rho(index), press(index));
    ghl_ppm_reconstruction_with_steepening(ghl_params, press_stencil, Gamma, ftilde, rho_stencil, &prims_r.rho, &prims_l.rho);

    ghl_ppm_reconstruction(ftilde, press_stencil, &prims_r.press, &prims_l.press);
    ghl_ppm_reconstruction(ftilde, vx_stencil, &prims_r.vU[0], &prims_l.vU[0]);
    ghl_ppm_reconstruction(ftilde, vy_stencil, &prims_r.vU[1], &prims_l.vU[1]);
    ghl_ppm_reconstruction(ftilde, vz_stencil, &prims_r.vU[2], &prims_l.vU[2]);
    ghl_ppm_reconstruction(ftilde, B1_stencil, &prims_r.BU[B1], &prims_l.BU[B1]);
    ghl_ppm_reconstruction(ftilde, B2_stencil, &prims_r.BU[B2], &prims_l.BU[B2]);

    // B_stagger is densitized, but B_center is not.
    prims_r.BU[B0] = prims_l.BU[B0] = B_stagger(p.I)/ADM_metric_face.sqrt_detgamma;

    const Loop::GF3D5index ind_flux(local_layout, p.I);
    // These if's should vanish when the compile
    // makes the actual templated functions.
    if(flux_dir==0) {
      v1_r(ind_flux) = prims_r.vU[0];
      v1_l(ind_flux) = prims_l.vU[0];
      v2_r(ind_flux) = prims_r.vU[1];
      v2_l(ind_flux) = prims_l.vU[1];
    }
    if(flux_dir==1) {
      v1_r(ind_flux) = prims_r.vU[1];
      v1_l(ind_flux) = prims_l.vU[1];
      v2_r(ind_flux) = prims_r.vU[2];
      v2_l(ind_flux) = prims_l.vU[2];
    }
    if(flux_dir==2) {
      v1_r(ind_flux) = prims_r.vU[2];
      v1_l(ind_flux) = prims_l.vU[2];
      v2_r(ind_flux) = prims_r.vU[0];
      v2_l(ind_flux) = prims_l.vU[0];
    }

    int speed_limited CCTK_ATTRIBUTE_UNUSED = ghl_limit_v_and_compute_u0(
          ghl_params, &ADM_metric_face, &prims_r);
    speed_limited = ghl_limit_v_and_compute_u0(
          ghl_params, &ADM_metric_face, &prims_l);

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

//template <int flux_dir>
//void IllinoisGRMHDX_reconstruct_for_A_dir(
//      CCTK_ARGUMENTS,
//      const Loop::GF3D5layout local_layout,
//      Loop::GF3D5<CCTK_REAL> v1_r,
//      Loop::GF3D5<CCTK_REAL> v1_l,
//      Loop::GF3D5<CCTK_REAL> v2_r,
//      Loop::GF3D5<CCTK_REAL> v2_l,
//      Loop::GF3D5<CCTK_REAL> v1_r,
//      Loop::GF3D5<CCTK_REAL> v1_l,
//      Loop::GF3D5<CCTK_REAL> v2_r,
//      Loop::GF3D5<CCTK_REAL> v2_l,
//      Loop::GF3D5<CCTK_REAL> cmin,
//      Loop::GF3D5<CCTK_REAL> cmax,
//      Loop::GF3D5<CCTK_REAL> rho_flux,
//      Loop::GF3D5<CCTK_REAL> tau_flux,
//      Loop::GF3D5<CCTK_REAL> S_x_flux,
//      Loop::GF3D5<CCTK_REAL> S_y_flux,
//      Loop::GF3D5<CCTK_REAL> S_z_flux) {
//  DECLARE_CCTK_ARGUMENTSX_IllinoisGRMHDX_hybrid_evaluate_fluxes_rhs;
//  DECLARE_CCTK_PARAMETERS;
//
//  const CCTK_REAL poison = 0.0/0.0;
//  const Loop::GF3D2layout ccc_layout(cctkGH, {1, 1, 1});
//
//}
//
extern "C" void IllinoisGRMHDX_hybrid_evaluate_fluxes_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_IllinoisGRMHDX_hybrid_evaluate_fluxes_rhs;
  DECLARE_CCTK_PARAMETERS;

  Arith::vect<int, Loop::dim> imin, imax;

  // Technically could reduce memory by only including gzs for
  // velocities, as these are the only ones that need them.

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

  IllinoisGRMHDX_evaluate_flux_dir<0>(cctkGH, vcc_layout, vx_xr, vx_xl, vy_xr, vy_xl, cmin_x, cmax_x, rho_flux_x, tau_flux_x, S_x_flux_x, S_y_flux_x, S_z_flux_x);
  IllinoisGRMHDX_evaluate_flux_dir<1>(cctkGH, cvc_layout, vy_yr, vy_yl, vz_yr, vz_yl, cmin_y, cmax_y, rho_flux_y, tau_flux_y, S_x_flux_y, S_y_flux_y, S_z_flux_y);
  IllinoisGRMHDX_evaluate_flux_dir<2>(cctkGH, ccv_layout, vz_yr, vz_zl, vx_zr, vx_zl, cmin_z, cmax_z, rho_flux_z, tau_flux_z, S_x_flux_z, S_y_flux_z, S_z_flux_z);

  // Set up temporary variables with vvc layout
  grid.box_int<0, 0, 1>(grid.nghostzones, imin, imax);
  const Loop::GF3D5layout vvc_layout(imin, imax);
  Loop::GF3D5vector<CCTK_REAL> tmpsvvc(vvc_layout, 4);
  itmp = 0; const auto make_gfvvc = [&]() { return Loop::GF3D5<CCTK_REAL>(tmpsvvc(itmp++)); };

  Loop::GF3D5<CCTK_REAL> vx_xyl(make_gfvvc());
  Loop::GF3D5<CCTK_REAL> vy_xyl(make_gfvvc());
  Loop::GF3D5<CCTK_REAL> Bx_xyl(make_gfvvc());
  Loop::GF3D5<CCTK_REAL> By_xyl(make_gfvvc());

  // Set up temporary variables with cvv layout
  grid.box_int<1, 0, 0>(grid.nghostzones, imin, imax);
  const Loop::GF3D5layout cvv_layout(imin, imax);
  Loop::GF3D5vector<CCTK_REAL> tmpscvv(cvv_layout, 4);
  itmp = 0; const auto make_gfcvv = [&]() { return Loop::GF3D5<CCTK_REAL>(tmpscvv(itmp++)); };

  Loop::GF3D5<CCTK_REAL> vy_yzl(make_gfcvv());
  Loop::GF3D5<CCTK_REAL> vz_yzl(make_gfcvv());
  Loop::GF3D5<CCTK_REAL> By_yzl(make_gfcvv());
  Loop::GF3D5<CCTK_REAL> Bz_yzl(make_gfcvv());

  // Set up temporary variables with vcv layout
  grid.box_int<0, 1, 0>(grid.nghostzones, imin, imax);
  const Loop::GF3D5layout vcv_layout(imin, imax);
  Loop::GF3D5vector<CCTK_REAL> tmpsvcv(vcv_layout, 4);
  itmp = 0; const auto make_gfvcv = [&]() { return Loop::GF3D5<CCTK_REAL>(tmpsvcv(itmp++)); };

  Loop::GF3D5<CCTK_REAL> vz_zxl(make_gfvcv());
  Loop::GF3D5<CCTK_REAL> vx_zxl(make_gfvcv());
  Loop::GF3D5<CCTK_REAL> Bz_zxl(make_gfvcv());
  Loop::GF3D5<CCTK_REAL> Bx_zxl(make_gfvcv());


}
