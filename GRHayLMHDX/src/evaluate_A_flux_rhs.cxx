#include "GRHayLMHDX.hxx"

template <int flux_dir>
void GRHayLMHDX_A_flux_rhs_dir(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLMHDX_evaluate_A_flux_rhs;

  /*
     All of these conditionals are simply selecting the correct
     components from the cross product. To summarize,

           |    RHS       | A_x | A_y | A_z |
           | Component 1  |  y  |  z  |  x  |
           | Component 2  |  z  |  x  |  y  |
  */

  Loop::GF3D2<const CCTK_REAL> v_flux_dir = flux_dir==0 ? vx :
                                            flux_dir==1 ? vy : vz;

  Loop::GF3D2<CCTK_REAL> A_rhs = flux_dir==0 ? Ax_rhs :
                                 flux_dir==1 ? Ay_rhs : Az_rhs;

  // Selecting index "1" variables
  Loop::GF3D2<const CCTK_REAL> v1r = flux_dir==0 ? vy_yr :
                                     flux_dir==1 ? vz_xr : vx_xr;

  Loop::GF3D2<const CCTK_REAL> v1l = flux_dir==0 ? vy_yl :
                                     flux_dir==1 ? vz_xl : vx_xl;

  Loop::GF3D2<const CCTK_REAL> B1 = flux_dir==0 ? By_stagger :
                                    flux_dir==1 ? Bz_stagger : Bx_stagger;

  // Selecting index "2" variables
  Loop::GF3D2<const CCTK_REAL> v2r = flux_dir==0 ? vz_yr :
                                     flux_dir==1 ? vx_xr : vy_xr;

  Loop::GF3D2<const CCTK_REAL> v2l = flux_dir==0 ? vz_yl :
                                     flux_dir==1 ? vx_xl : vy_xl;

  Loop::GF3D2<const CCTK_REAL> B2 = flux_dir==0 ? Bz_stagger :
                                    flux_dir==1 ? Bx_stagger : By_stagger;

  // Selecting characteristic speed variables
  Loop::GF3D2<const CCTK_REAL> cmin1 = flux_dir==0 ? cmin_y :
                                       flux_dir==1 ? cmin_z : cmin_x;

  Loop::GF3D2<const CCTK_REAL> cmax1 = flux_dir==0 ? cmax_y :
                                       flux_dir==1 ? cmax_z : cmax_x;

  Loop::GF3D2<const CCTK_REAL> cmin2 = flux_dir==0 ? cmin_z :
                                       flux_dir==1 ? cmin_x : cmin_y;

  Loop::GF3D2<const CCTK_REAL> cmax2 = flux_dir==0 ? cmax_z :
                                       flux_dir==1 ? cmax_x : cmax_y;

  constexpr int recon1 = flux_dir==0 ? 2 :
                         flux_dir==1 ? 2 : 1;

  constexpr int recon2 = flux_dir==0 ? 1 :
                         flux_dir==1 ? 0 : 0;

  const Loop::GF3D2layout ccc_layout(cctkGH, {1, 1, 1});

  constexpr std::array<int, Loop::dim> edgetype = {flux_dir==0, flux_dir==1, flux_dir==2};
  const Loop::GF3D2layout edge_layout(cctkGH, edgetype);

  grid.loop_int_device<edgetype[0], edgetype[1], edgetype[2]>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index c_index(ccc_layout, p.I);
    const Loop::GF3D2index e_index(edge_layout, p.I);

    HLL_2D_vars vars;

    CCTK_REAL press_stencil[6], v_flux[6];
    CCTK_REAL vars_stencil[5][6], vars_r[5], vars_l[5];

    for(int ind=0; ind<6; ind++) {
      // Stencil from -3 to +2 reconstructs to e.g. i-1/2
      const auto stencil = p.I + (ind-3)*p.DI[recon1];
      const Loop::GF3D2index c_stencil(ccc_layout, stencil);
      v_flux[ind] = v_flux_dir(c_stencil); // Could be smaller; doesn't use full stencil
      press_stencil[ind] = pressure(c_stencil);
      vars_stencil[0][ind] = v1r(stencil);
      vars_stencil[1][ind] = v1l(stencil);
      vars_stencil[2][ind] = v2r(stencil);
      vars_stencil[3][ind] = v2l(stencil);
      //B swapped due to how y direction v_{r,l} are saved
      vars_stencil[4][ind] = flux_dir==1 ? B2(stencil) : B1(stencil);
    }

    // Compute Gamma_eff
    CCTK_REAL K, Gamma;
    ghl_hybrid_get_K_and_Gamma(ghl_eos, rho_b(c_index), &K, &Gamma);
    const CCTK_REAL P_cold = K*pow(rho_b(c_index), Gamma);
    const CCTK_REAL Gamma_eff = ghl_eos->Gamma_th + (Gamma - ghl_eos->Gamma_th)*P_cold/pressure(c_index);

    ghl_ppm_no_rho_P(
          press_stencil, vars_stencil,
          5, v_flux, Gamma_eff, vars_r, vars_l);

    vars.v1rr = vars_r[0];
    vars.v1ll = vars_l[1];
    vars.v2rr = vars_r[2];
    vars.v2ll = vars_l[3];
    /*
       v_rl = left recon. of v_r
       v_lr = right recon. of v_l
       Note, however, that this assumes a specific order of reconstruction
       that is not followed by the y coordinate, which should reconstruct
       z followed by x. As such, we have to swap the first and second
       letter for the mixed faces for this direction. Additionally,
       swapping the order of reconstruction means that B2 is being reconstructed
       in the same direction instead of B1.
    */
    if(flux_dir==1) {
      vars.v1rl = vars_l[1];
      vars.v1lr = vars_r[0];
      vars.v2rl = vars_l[3];
      vars.v2lr = vars_r[2];
      vars.B2r  = vars_r[4];
      vars.B2l  = vars_l[4];
    } else {
      vars.v1rl = vars_l[0];
      vars.v1lr = vars_r[1];
      vars.v2rl = vars_l[2];
      vars.v2lr = vars_r[3];
      vars.B1r  = vars_r[4];
      vars.B1l  = vars_l[4];
    }

    for(int ind=0; ind<6; ind++) {
      // Stencil from -3 to +2 reconstructs to e.g. i-1/2
      const auto stencil = p.I + (ind-3)*p.DI[recon2];
      const Loop::GF3D2index c_stencil(ccc_layout, stencil);
      v_flux[ind] = v_flux_dir(c_stencil); // Could be smaller; doesn't use full stencil
      press_stencil[ind] = pressure(c_stencil);
      //B swapped due to how y direction v_{r,l} are saved
      vars_stencil[0][ind] = flux_dir==1 ? B1(stencil) : B2(stencil);
    }

    ghl_ppm_no_rho_P(
          press_stencil, vars_stencil,
          1, v_flux, Gamma_eff, vars_r, vars_l);

    if(flux_dir==1) {
      vars.B1r = vars_r[0];
      vars.B1l = vars_l[0];
    } else {
      vars.B2r = vars_r[0];
      vars.B2l = vars_l[0];
    }

    /*
      Note that cmax/cmin (\alpha^{\pm}  as defined in Del Zanna et al) is at a slightly DIFFERENT
      point (e.g., (i+1/2,j,k) instead of (i+1/2,j+1/2,k) for F3). Yuk Tung Liu discussed this point
      with M. Shibata, who found that the effect is negligible.
    */
    vars.c1_min = cmin1(p.I);
    vars.c1_max = cmax1(p.I);
    vars.c2_min = cmin2(p.I);
    vars.c2_max = cmax2(p.I);

    // Reconstructed variables are declared vvv, but the RHS are properly set to edge centering
    A_rhs(e_index) = ghl_HLL_2D_flux_with_Btilde(&vars);
  });
}

extern "C" void GRHayLMHDX_evaluate_A_flux_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLMHDX_evaluate_A_flux_rhs;
  DECLARE_CCTK_PARAMETERS;

  GRHayLMHDX_A_flux_rhs_dir<0>(cctkGH);
  GRHayLMHDX_A_flux_rhs_dir<1>(cctkGH);
  GRHayLMHDX_A_flux_rhs_dir<2>(cctkGH);
}
