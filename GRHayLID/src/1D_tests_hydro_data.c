#include "GRHayLID.h"

void GRHayLID_1D_tests_hydro_data(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLID_1D_tests_hydro_data;
  DECLARE_CCTK_PARAMETERS;

  if(!(CCTK_EQUALS(EOS_type, "Simple") || CCTK_EQUALS(EOS_type, "Hybrid")))
    CCTK_ERROR("HydroTest1D requires EOS_type=Simple or Hybrid");
  GRHayLID_require_storage(cctkGH, "ADMBase::metric");
  if(CCTK_EQUALS(initial_data_1D, "sound wave") &&
     (!isfinite(wave_amplitude) || wave_amplitude < 0.0 || wave_amplitude >= 1.0))
    CCTK_ERROR("Sound-wave velocity amplitude must be finite and in [0,1)");

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
  } else if(CCTK_EQUALS(initial_data_1D,"sound wave")) {
    // The sinusoidal longitudinal velocity is filled pointwise below.
    rho_l = rho_r = press_l = press_r = 1.0;
    vx_l = vy_l = vz_l = vx_r = vy_r = vz_r = 0.0;
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

#pragma omp parallel for
  for(int k=0; k<cctk_lsh[2]; k++) {
    for(int j=0; j<cctk_lsh[1]; j++) {
      for(int i=0; i<cctk_lsh[0]; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        GRHayLID_check_flat_metric(gxx[index], gxy[index], gxz[index], gyy[index], gyz[index], gzz[index]);
        const int ind4x = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0);
        const int ind4y = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1);
        const int ind4z = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2);

        double step = x[index];
        if(CCTK_EQUALS(shock_direction, "y")) {
          step = y[index];
        } else if(CCTK_EQUALS(shock_direction, "z")) {
          step = z[index];
        }

        if(CCTK_EQUALS(initial_data_1D,"sound wave")) {
          rho[index]   = 1.0;
          press[index] = 1.0; // Thermal pressure excludes kinetic energy.
          vel[ind4x]   = 0.0;
          vel[ind4y]   = 0.0;
          vel[ind4z]   = 0.0;
          const int axis = CCTK_EQUALS(shock_direction, "y") ? ind4y :
                           CCTK_EQUALS(shock_direction, "z") ? ind4z : ind4x;
          vel[axis] = wave_amplitude * sin(M_PI * step);
        } else if(step <= discontinuity_position) {
          rho[index]   = rho_l;
          press[index] = press_l;
          vel[ind4x]   = vx_l;
          vel[ind4y]   = vy_l;
          vel[ind4z]   = vz_l;
        } else {
          rho[index]   = rho_r;
          press[index] = press_r;
          vel[ind4x]   = vx_r;
          vel[ind4y]   = vy_r;
          vel[ind4z]   = vz_r;
        }
        if(CCTK_EQUALS(EOS_type, "Hybrid")) {
          double p_cold, eps_cold;
          ghl_hybrid_compute_P_cold_and_eps_cold(ghl_eos, rho[index], &p_cold, &eps_cold);
          if(press[index] < p_cold)
            CCTK_ERROR("HydroTest1D pressure is below the Hybrid cold curve; choose admissible EOS parameters");
          eps[index] = ghl_hybrid_compute_epsilon(ghl_eos, rho[index], press[index]);
        } else {
          const double Gamma = ghl_eos->Gamma_ppoly[0];
          eps[index] = press[index]/(rho[index]*(Gamma-1.0));
        }
        if(!isfinite(eps[index]))
          CCTK_ERROR("HydroTest1D EOS produced nonfinite internal energy");
      }
    }
  }
  CCTK_VINFO("Finished writing hydrodynamic ID for %s test", initial_data_1D);
}
