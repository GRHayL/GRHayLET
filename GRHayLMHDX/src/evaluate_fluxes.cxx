#include "GRHayLMHDX.hxx"

template <int flux_dir>
void GRHayLMHDX_evaluate_flux_dir(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLMHDX_evaluate_fluxes;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL poison = 0.0/0.0;
  const Loop::GF3D2layout ccc_layout(cctkGH, {1, 1, 1});

  constexpr std::array<int, Loop::dim> facetype = {flux_dir!=0, flux_dir!=1, flux_dir!=2};
  const Loop::GF3D2layout flux_layout(cctkGH, facetype);

  // These nested condtional ternary operators let us tell the compiler that
  // these pointers can be set at compiler time while making the templated functions
  // instead of using a switch statement at runtime.
  constexpr void (*calculate_characteristic_speed)(const primitive_quantities *restrict prims_r,
                                         const primitive_quantities *restrict prims_l,
                                         struct eos_parameters const *restrict eos,
                                         const metric_quantities *restrict ADM_metric_face,
                                         CCTK_REAL *cmin, CCTK_REAL *cmax)
    = flux_dir==0 ? &ghl_calculate_characteristic_speed_dirn0 :
      flux_dir==1 ? &ghl_calculate_characteristic_speed_dirn1 :
                    &ghl_calculate_characteristic_speed_dirn2 ;

  constexpr void (*calculate_HLLE_fluxes)(const primitive_quantities *restrict prims_r,
                                const primitive_quantities *restrict prims_l,
                                const eos_parameters *restrict eos,
                                const metric_quantities *restrict ADM_metric_face,
                                const CCTK_REAL cmin,
                                const CCTK_REAL cmax,
                                conservative_quantities *restrict cons_fluxes)
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

  Loop::GF3D2<CCTK_REAL> cmin = flux_dir==0 ? cmin_x :
                                flux_dir==1 ? cmin_y : cmin_z;

  Loop::GF3D2<CCTK_REAL> cmax = flux_dir==0 ? cmax_x :
                                flux_dir==1 ? cmax_y : cmax_z;

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
/*
  template <int CI, int CJ, int CK, int VS = 1, typename F>
  void loop_box_device(const F &f, const vect<int, dim> &restrict bnd_min,
                       const vect<int, dim> &restrict bnd_max,
                       const vect<int, dim> &restrict loop_min,
                       const vect<int, dim> &restrict loop_max) const {

  template <int CI, int CJ, int CK, int VS = 1, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_int_device(const vect<int, dim> &group_nghostzones, const F &f) const {
*/
  Arith::vect<int, 3> bnd_min, bnd_max;
  grid.boundary_box<facetype[0], facetype[1], facetype[2]>(grid.nghostzones, bnd_min, bnd_max);
  Arith::vect<int, 3> imin, imax;
  grid.box_int<facetype[0], facetype[1], facetype[2]>(grid.nghostzones, imin, imax);

  // Extend the loop box. You must make sure that you extend by at most the
  // number of ghosts.
  for (int d = 0; d < 3; ++d) {
    if(d==flux_dir) continue;
    imin[d] -= 3;
    imax[d] += 3;
  }

//  grid.loop_box_device<0, 0, 0>(
//      bnd_min, bnd_max, imin, imax,
//      [=] CCTK_DEVICE(const Loop::PointDesc &p)
//          CCTK_ATTRIBUTE_ALWAYS_INLINE { printf("a"); });
  // These variables let us compute all the velocities we need for A_RHS
  // without stencils going out of bounds. This does compute the fluxes
  // on some ghost zones where we don't need them, but this is mostly
  // unavoidable and faster than including more if statements.
  const int loop_min = cctkGH->cctk_nghostzones[flux_dir];
  const int loop_max = cctkGH->cctk_lsh[flux_dir] - cctkGH->cctk_nghostzones[flux_dir];

  grid.loop_all_device<facetype[0], facetype[1], facetype[2]>(
      grid.nghostzones,
  //grid.loop_box_device<facetype[0], facetype[1], facetype[2]>(
  //    bnd_min, bnd_max, imin, imax,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

    // Honestly, this is really hacky. I'd much prefer a way to loop over just
    // some ghost zones.
    const int loop_indices[3] = {p.i, p.j, p.k};
    if(loop_indices[flux_dir] >= loop_min && loop_indices[flux_dir] < loop_max) {
      const Loop::GF3D2index indm2(ccc_layout, p.I - 2*p.DI[flux_dir]);
      const Loop::GF3D2index indm1(ccc_layout, p.I - p.DI[flux_dir]);
      const Loop::GF3D2index index(ccc_layout, p.I);
      const Loop::GF3D2index indp1(ccc_layout, p.I + p.DI[flux_dir]);

      const Loop::GF3D2index ind_flux(flux_layout, p.I);

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
  
      metric_quantities ADM_metric_face;
      GRHayLMHDX_interpolate_metric_to_face(
            indm2, indm1, index, indp1,
            ccc_lapse, ccc_betax, ccc_betay, ccc_betaz,
            ccc_gxx, ccc_gxy, ccc_gxz,
            ccc_gyy, ccc_gyz, ccc_gzz,
            &ADM_metric_face);

      // These if's should vanish when the compile
      // makes the actual templated functions.
      if(flux_dir==0) {
        vx_xr(ind_flux) = vars_r[0];
        vy_xr(ind_flux) = vars_r[1];
        vz_xr(ind_flux) = vars_r[2];
        vx_xl(ind_flux) = vars_l[0];
        vy_xl(ind_flux) = vars_l[1];
        vz_xl(ind_flux) = vars_l[2];
      } else if(flux_dir==1) {
        vy_yr(ind_flux) = vars_r[1];
        vz_yr(ind_flux) = vars_r[2];
        vy_yl(ind_flux) = vars_l[1];
        vz_yl(ind_flux) = vars_l[2];
      }
  
      B_r[B1] = vars_r[3];
      B_r[B2] = vars_r[4];

      B_l[B1] = vars_l[3];
      B_l[B2] = vars_l[4];

      // B_stagger is densitized, but B_center is not.
      B_r[B0] = B_l[B0] = B_stagger(p.I)/ADM_metric_face.sqrt_detgamma;

      primitive_quantities prims_r, prims_l;
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
  
      conservative_quantities cons_fluxes;
      calculate_characteristic_speed(&prims_r, &prims_l, ghl_eos, &ADM_metric_face, &cmin(ind_flux), &cmax(ind_flux));
//if(isnan(cmin(ind_flux)*cmax(ind_flux)))
//CCTK_VINFO("c %e %e\n"
//"r %e %e v %e %e %e B %e %e %e u %e\n"
//"l %e %e v %e %e %e B %e %e %e u %e\n",
//cmin(ind_flux), cmax(ind_flux),
//rhor, pressr, vars_r[0], vars_r[1], vars_r[2], B_r[0], B_r[1], B_r[2], prims_r.u0,
//rhol, pressl, vars_l[0], vars_l[1], vars_l[2], B_l[0], B_l[1], B_l[2], prims_l.u0);
      calculate_HLLE_fluxes(&prims_r, &prims_l, ghl_eos, &ADM_metric_face, cmin(ind_flux), cmax(ind_flux), &cons_fluxes);

      rho_star_flux(ind_flux) = cons_fluxes.rho;
      tau_flux(ind_flux) = cons_fluxes.tau;
      Sx_flux(ind_flux)  = cons_fluxes.SD[0];
      Sy_flux(ind_flux)  = cons_fluxes.SD[1];
      Sz_flux(ind_flux)  = cons_fluxes.SD[2];
    } else {
      const Loop::GF3D2index ind_flux(flux_layout, p.I);
      rho_star_flux(ind_flux) = 1e300;
      tau_flux(ind_flux) = 1e300;
      Sx_flux(ind_flux)  = 1e300;
      Sy_flux(ind_flux)  = 1e300;
      Sz_flux(ind_flux)  = 1e300;
      cmin(ind_flux) = 1e300;
      cmax(ind_flux) = 1e300;
      if(flux_dir==0) {
        vx_xr(ind_flux) = 1e300;
        vy_xr(ind_flux) = 1e300;
        vz_xr(ind_flux) = 1e300;
        vx_xl(ind_flux) = 1e300;
        vy_xl(ind_flux) = 1e300;
        vz_xl(ind_flux) = 1e300;
      } else if(flux_dir==1) {
        vy_yr(ind_flux) = 1e300;
        vz_yr(ind_flux) = 1e300;
        vy_yl(ind_flux) = 1e300;
        vz_yl(ind_flux) = 1e300;
      }
    }
  }); // staggered loop interior (e.g. flux_dir=0 gives vcc)
}

extern "C" void GRHayLMHDX_evaluate_fluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLMHDX_evaluate_fluxes;
  DECLARE_CCTK_PARAMETERS;

  GRHayLMHDX_evaluate_flux_dir<0>(cctkGH);
  GRHayLMHDX_evaluate_flux_dir<1>(cctkGH);
  GRHayLMHDX_evaluate_flux_dir<2>(cctkGH);
}
