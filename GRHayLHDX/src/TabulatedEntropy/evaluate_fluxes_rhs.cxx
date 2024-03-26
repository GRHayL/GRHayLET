#include "GRHayLHDX.h"

template <int flux_dir>
void GRHayLHDX_tabulated_entropy_evaluate_fluxes_rhs_dir(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_tabulated_entropy_evaluate_fluxes_rhs;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D2layout ccc_layout(cctkGH, {1, 1, 1});

  constexpr std::array<int, Loop::dim> facetype = {flux_dir!=0, flux_dir!=1, flux_dir!=2};
  const Loop::GF3D2layout flux_layout(cctkGH, facetype);

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

  Loop::GF3D2<const CCTK_REAL> ent_flux = flux_dir==0 ? ent_flux_x :
                                          flux_dir==1 ? ent_flux_y : ent_flux_z;

  Loop::GF3D2<const CCTK_REAL> Ye_flux = flux_dir==0 ? Ye_flux_x :
                                         flux_dir==1 ? Ye_flux_y : Ye_flux_z;

  const CCTK_REAL dxi = 1.0/CCTK_DELTA_SPACE(flux_dir);

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index index(ccc_layout, p.I);
    const Loop::GF3D2index ind_flux(flux_layout, p.I);
    const Loop::GF3D2index ind_flp1(flux_layout, p.I + p.DI[flux_dir]);

    rho_star_rhs(index) += dxi*(rho_star_flux(ind_flux) - rho_star_flux(ind_flp1));
    tau_rhs(index)      += dxi*(tau_flux(ind_flux)      - tau_flux(ind_flp1));
    Stildex_rhs(index)  += dxi*(Sx_flux(ind_flux)       - Sx_flux(ind_flp1));
    Stildey_rhs(index)  += dxi*(Sy_flux(ind_flux)       - Sy_flux(ind_flp1));
    Stildez_rhs(index)  += dxi*(Sz_flux(ind_flux)       - Sz_flux(ind_flp1));
    ent_star_rhs(index) += dxi*(ent_flux(ind_flux)      - ent_flux(ind_flp1));
    Ye_star_rhs (index) += dxi*(Ye_flux (ind_flux)      - Ye_flux(ind_flp1));
  }); // ccc loop interior
}

extern "C" void GRHayLHDX_tabulated_entropy_evaluate_fluxes_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_tabulated_entropy_evaluate_fluxes_rhs;
  DECLARE_CCTK_PARAMETERS;

  GRHayLHDX_tabulated_entropy_evaluate_fluxes_rhs_dir<0>(cctkGH);
  GRHayLHDX_tabulated_entropy_evaluate_fluxes_rhs_dir<1>(cctkGH);
  GRHayLHDX_tabulated_entropy_evaluate_fluxes_rhs_dir<2>(cctkGH);
}
