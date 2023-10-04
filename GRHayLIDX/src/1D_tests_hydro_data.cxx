#include "GRHayLIDX.h"

extern "C" void GRHayLIDX_1D_tests_hydro_data(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLIDX_1D_tests_hydro_data;
  DECLARE_CCTK_PARAMETERS;

  if(CCTK_EQUALS(EOS_type, "Tabulated"))
    CCTK_ERROR("1D test initial data is only defined for hybrid or ideal fluid EOS, and the standard comparison uses the ideal fluid EOS. Please change GRHayLib::EOS_type in the parfile.");

  const Loop::GF3D2layout layout(cctkGH, {1, 1, 1});

  double rho_l, rho_r;
  double press_l, press_r;
  double vx_l, vy_l, vz_l;
  double vx_r, vy_r, vz_r;
  if(CCTK_EQUALS(initial_data_1D,"Balsara1")) {
    rho_l = 1.0;
    rho_r = 0.125;
    press_l = 1.0;
    press_r = 0.1;
    vx_l = vy_l = vz_l = 0.0;
    vx_r = vy_r = vz_r = 0.0;
  } else if(CCTK_EQUALS(initial_data_1D,"Balsara2")) {
    rho_l = rho_r = 1.0;
    press_l = 30.0;
    press_r = 1.0;
    vx_l = vy_l = vz_l = 0.0;
    vx_r = vy_r = vz_r = 0.0;
  } else if(CCTK_EQUALS(initial_data_1D,"Balsara3")) {
    rho_l = rho_r = 1.0;
    press_l = 1000.0;
    press_r = 0.1;
    vx_l = vy_l = vz_l = 0.0;
    vx_r = vy_r = vz_r = 0.0;
  } else if(CCTK_EQUALS(initial_data_1D,"Balsara4")) {
    rho_l = rho_r = 1.0;
    press_l = press_r = 0.1;
    vx_l = 0.999;
    vx_r = -0.999;
    vy_l = vz_l = 0.0;
    vy_r = vz_r = 0.0;
  } else if(CCTK_EQUALS(initial_data_1D,"Balsara5")) {
    rho_l = 1.08;
    rho_r = 1.0;
    press_l = 0.95;
    press_r = 1.0;
    vx_l = 0.4;
    vy_l = 0.3;
    vx_r = -0.45;
    vy_r = -0.2;
    vz_l = vz_r = 0.2;
  } else if(CCTK_EQUALS(initial_data_1D,"equilibrium")) {
    rho_l = rho_r = 1.0;
    press_l = press_r = 1.0;
    vx_l = vy_l = vz_l = 0.0;
    vx_r = vy_r = vz_r = 0.0;
  /*
  } else if(CCTK_EQUALS(initial_data_1D,"sound wave")) {
    this case is handled in the loop because it isn't a
    step function but a sin() wave
  */
  } else if(CCTK_EQUALS(initial_data_1D,"shock tube")) {
    rho_l = 2.0;
    rho_r = 1.0;
    press_l = 2.0;
    press_r = 1.0;
    vx_l = vy_l = vz_l = 0.0;
    vx_r = vy_r = vz_r = 0.0;
  } else {
    CCTK_VERROR("Parameter initial_data_1D is not set "
                "to a valid test name. Something has gone wrong.");
  }

  //Above data assumes that the shock is in the x direction; we just have to rotate
  //the data for it to work in other directions.
  if(CCTK_EQUALS(shock_direction, "y")) {
    //x-->y, y-->z, z-->x
    const double vl[3] = {vz_l, vx_l, vy_l};
    vx_l = vl[0]; vy_l = vl[1]; vz_l = vl[2];

    const double vr[3] = {vz_r, vx_r, vy_r};
    vx_r = vr[0]; vy_r = vr[1]; vz_r = vr[2];
  } else if(CCTK_EQUALS(shock_direction, "z")) {
    //x-->z, y-->x, z-->y
    const double vl[3] = {vy_l, vz_l, vx_l};
    vx_l = vl[0]; vy_l = vl[1]; vz_l = vl[2];

    const double vr[3] = {vy_r, vz_r, vx_r};
    vx_r = vr[0]; vy_r = vr[1]; vz_r = vr[2];
  }

  grid.loop_all<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index index(layout, p.I);

    double step = p.x;
    if(CCTK_EQUALS(shock_direction, "y")) {
      step = p.y;
    } else if(CCTK_EQUALS(shock_direction, "z")) {
      step = p.z;
    }

    if(CCTK_EQUALS(initial_data_1D,"sound wave")) {
      rho(index)   = 1.0;
      press(index) = 1.0; // should add kinetic energy here
      velx(index)  = wave_amplitude * sin(M_PI * step);
      vely(index)  = velz(index)  = 0.0;
    } else if(step <= discontinuity_position) {
      rho(index)   = rho_l;
      press(index) = press_l;
      velx(index)  = vx_l;
      vely(index)  = vy_l;
      velz(index)  = vz_l;
    } else {
      rho(index)   = rho_r;
      press(index) = press_r;
      velx(index)  = vx_r;
      vely(index)  = vy_r;
      velz(index)  = vz_r;
    }
    const double Gamma = ghl_eos->Gamma_ppoly[
                                  ghl_hybrid_find_polytropic_index(
                                              ghl_eos, rho(index))];
    eps(index) = press(index)/( rho(index)*(Gamma-1) );
  });
}
