extern "C" void IllinoisGRMHDX_reconstruct_velocities(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_IllinoisGRMHDX_reconstruct_velocities;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D2layout ccc_layout(cctkGH, {1,1,1});
  const Loop::GF3D2layout vcc_layout(cctkGH, {0,1,1});
  const Loop::GF3D2layout cvc_layout(cctkGH, {1,0,1});

  /*
     Note that we need more than just the interior for these variables.
     This is because the second reconstruction will require some
     ghost zones to be filled. Hence, we loop 'everywhere', but
     short-circuit the loop with a continue if we are in the ghost zones
     for the direction of the reconstruction.
  */
  const int imin = cctkGH->cctk_nghostzones[0];
  const int imax = cctkGH->cctk_lsh[0] - cctkGH->cctk_nghostzones[0];
  grid.loop_all_device<0, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    if(p.i < imin || p.i > imax) continue;
    const Loop::GF3D2index vel_index(vcc_layout, p.I);

    double press_stencil[6], v_flux_stencil[6];
    double var_data[3][6], vars_r[3], vars_l[3];

    for(int ind=0; ind<6; ind++) {
      // Stencil from -3 to +2 reconstructs to e.g. i-1/2
      const Loop::GF3D2index stencil(ccc_layout, p.I + (ind-3)*p.DI[flux_dir]);
      v_flux_stencil[ind] = vy(stencil); // Could be smaller; doesn't use full stencil
      press_stencil[ind] = pressure(stencil);
      var_data[0][ind] = vx(stencil);
      var_data[1][ind] = vy(stencil);
      var_data[2][ind] = vz(stencil);
    }

    // Compute Gamma_eff
    const double rho_val = rho_b(p.I);
    double K, Gamma;
    ghl_hybrid_get_K_and_Gamma(ghl_eos, rho_val, &K, &Gamma);
    const double P_cold = K*pow(rho_val, Gamma);
    const double Gamma_eff = ghl_eos->Gamma_th + (Gamma - ghl_eos->Gamma_th)*P_cold/pressure(p.I);

    ghl_ppm_no_rho_P(
          press_stencil, var_data,
          3, v_flux_stencil, Gamma_eff,
          vars_r, vars_l);

    vx_xr(vel_index) = vars_r[0];
    vy_xr(vel_index) = vars_r[1];
    vz_xr(vel_index) = vars_r[2];

    vx_xl(vel_index) = vars_l[0];
    vy_xl(vel_index) = vars_l[1];
    vz_xl(vel_index) = vars_l[2];
  });

  const int jmin = cctkGH->cctk_nghostzones[1];
  const int jmax = cctkGH->cctk_lsh[1] - cctkGH->cctk_nghostzones[1];
  grid.loop_int_device<1, 0, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    if(p.j < jmin || p.j > jmax) continue;
    const Loop::GF3D2index vel_index(cvc_layout, p.I);

    double press_stencil[6], v_flux_stencil[6];
    double var_data[2][6], vars_r[2], vars_l[2];

    for(int ind=0; ind<6; ind++) {
      // Stencil from -3 to +2 reconstructs to e.g. i-1/2
      const Loop::GF3D2index stencil(ccc_layout, p.I + (ind-3)*p.DI[flux_dir]);
      v_flux_stencil[ind] = vy(stencil); // Could be smaller; doesn't use full stencil
      press_stencil[ind] = pressure(stencil);
      var_data[0][ind] = vy(stencil);
      var_data[1][ind] = vz(stencil);
    }

    // Compute Gamma_eff
    const double rho_val = rho_b(p.I);
    double K, Gamma;
    ghl_hybrid_get_K_and_Gamma(ghl_eos, rho_val, &K, &Gamma);
    const double P_cold = K*pow(rho_val, Gamma);
    const double Gamma_eff = ghl_eos->Gamma_th + (Gamma - ghl_eos->Gamma_th)*P_cold/pressure(p.I);

    ghl_ppm_no_rho_P(
          press_stencil, var_data,
          2, v_flux_stencil, Gamma_eff,
          vars_r, vars_l);

    vy_yr(vel_index) = vars_r[0];
    vz_yr(vel_index) = vars_r[1];

    vy_yl(vel_index) = vars_l[0];
    vz_yl(vel_index) = vars_l[1];
  });
}
