#include "GRHayLHDX.h"

  /*
   *  Computation of \partial_i on RHS of \partial_t {rho_star,tau,Stilde{x,y,z}},
   *  via PPM reconstruction onto e.g. (i+1/2,j,k), so that
   *  \partial_x F = [ F(i+1/2,j,k) - F(i-1/2,j,k) ] / dx
   */
template <int flux_dir>
void GRHayLHDX_evaluate_flux_source_rhs_dir(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_evaluate_flux_source_rhs;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D2layout ccc_layout(cctkGH, {1, 1, 1});

  constexpr std::array<int, Loop::dim> facetype = {flux_dir!=0, flux_dir!=1, flux_dir!=2};
  const Loop::GF3D2layout flux_layout(cctkGH, facetype);

  // These nested condtional ternary operators let us tell the compiler that
  // these pointers can be set at compiler time while making the templated functions
  // instead of using a switch statement at runtime.
  constexpr void (*calculate_source_terms)(
        ghl_primitive_quantities *restrict prims,
        const ghl_eos_parameters *restrict eos,
        const ghl_metric_quantities *restrict ADM_metric,
        const ghl_metric_quantities *restrict metric_derivs,
        ghl_conservative_quantities *restrict cons_sources)
    = flux_dir==0 ? &ghl_calculate_source_terms_dirn0 :
      flux_dir==1 ? &ghl_calculate_source_terms_dirn1 :
                    &ghl_calculate_source_terms_dirn2 ;

  Loop::GF3D2<const CCTK_REAL> rho_star_flux = flux_dir==0 ? rho_star_flux_x :
                                               flux_dir==1 ? rho_star_flux_y : rho_star_flux_z;

  Loop::GF3D2<const CCTK_REAL> tau_flux = flux_dir==0 ? tau_flux_x :
                                          flux_dir==1 ? tau_flux_y : tau_flux_z;

  Loop::GF3D2<const CCTK_REAL> Sx_flux = flux_dir==0 ? Sx_flux_x :
                                         flux_dir==1 ? Sx_flux_y : Sx_flux_z;

  Loop::GF3D2<const CCTK_REAL> Sy_flux = flux_dir==0 ? Sy_flux_x :
                                         flux_dir==1 ? Sy_flux_y : Sy_flux_z;

  Loop::GF3D2<const CCTK_REAL> Sz_flux = flux_dir==0 ? Sz_flux_x :
                                         flux_dir==1 ? Sz_flux_y : Sz_flux_z;

  //Loop::GF3D2<const CCTK_REAL> ent_flux = flux_dir==0 ? ent_flux_x :
  //                                        flux_dir==1 ? ent_flux_y : ent_flux_z;

  //Loop::GF3D2<const CCTK_REAL> Ye_flux = flux_dir==0 ? Ye_flux_x :
  //                                       flux_dir==1 ? Ye_flux_y : Ye_flux_z;

  const CCTK_REAL dxi = 1.0/CCTK_DELTA_SPACE(flux_dir);

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index indm2(ccc_layout, p.I - 2*p.DI[flux_dir]);
    const Loop::GF3D2index indm1(ccc_layout, p.I - p.DI[flux_dir]);
    const Loop::GF3D2index index(ccc_layout, p.I);
    const Loop::GF3D2index indp1(ccc_layout, p.I + p.DI[flux_dir]);
    const Loop::GF3D2index indp2(ccc_layout, p.I + 2*p.DI[flux_dir]);

    const Loop::GF3D2index ind_flux(flux_layout, p.I);
    const Loop::GF3D2index ind_flp1(flux_layout, p.I + p.DI[flux_dir]);

    rho_star_rhs(index) += dxi*(rho_star_flux(ind_flux) - rho_star_flux(ind_flp1));
    tau_rhs(index)      += dxi*(tau_flux(ind_flux)      - tau_flux(ind_flp1));
    Stildex_rhs(index)  += dxi*(Sx_flux(ind_flux)       - Sx_flux(ind_flp1));
    Stildey_rhs(index)  += dxi*(Sy_flux(ind_flux)       - Sy_flux(ind_flp1));
    Stildez_rhs(index)  += dxi*(Sz_flux(ind_flux)       - Sz_flux(ind_flp1));
    //ent_star_rhs(index) += dxi*(ent_flux(ind_flux)      - ent_flux(ind_flp1));
    //Ye_star_rhs (index) += dxi*(Ye_flux (ind_flux)      - Ye_flux(ind_flp1));

    ghl_metric_quantities ADM_metric;
    ghl_initialize_metric(ccc_lapse(index),
          ccc_betax(index), ccc_betay(index), ccc_betaz(index),
          ccc_gxx(index), ccc_gxy(index), ccc_gxz(index),
          ccc_gyy(index), ccc_gyz(index), ccc_gzz(index),
          &ADM_metric);

    ghl_primitive_quantities prims;
    ghl_initialize_primitives(
          rho(index), press(index), eps(index),
          vx(index), vy(index), vz(index),
          0.0, 0.0, 0.0,
          0.0, 0.0, 0.0,
    //      entropy(index), Ye(index), temperature(index),
          &prims);

    const int speed_limited CCTK_ATTRIBUTE_UNUSED = ghl_limit_v_and_compute_u0(ghl_eos, &ADM_metric, &prims);

    ghl_metric_quantities ADM_metric_derivs;

    GRHayLHDX_compute_metric_derivs(
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

extern "C" void GRHayLHDX_evaluate_flux_source_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_evaluate_flux_source_rhs;
  DECLARE_CCTK_PARAMETERS;

  GRHayLHDX_evaluate_flux_source_rhs_dir<0>(cctkGH);
  GRHayLHDX_evaluate_flux_source_rhs_dir<1>(cctkGH);
  GRHayLHDX_evaluate_flux_source_rhs_dir<2>(cctkGH);
}
