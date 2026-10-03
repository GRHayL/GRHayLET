#include "GRHayLID.h"

void GRHayLID_compute_entropy_hybrid(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLID_compute_entropy_hybrid;
  DECLARE_CCTK_PARAMETERS;
  GRHayLID_require_storage(cctkGH, "HydroBase::entropy");
  if(!allow_native_entropy_proxy)
    CCTK_ERROR("Native Simple/Hybrid entropy proxy requires explicit compatible-consumer opt-in");

#pragma omp parallel for
  for(int k=0; k<cctk_lsh[2]; k++) {
    for(int j=0; j<cctk_lsh[1]; j++) {
      for(int i=0; i<cctk_lsh[0]; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        if(!isfinite(rho[index]) || rho[index] <= 0.0 ||
           !isfinite(press[index]) || press[index] < 0.0)
          CCTK_ERROR("Native entropy requires positive finite density and nonnegative finite pressure");
        const double value = ghl_hybrid_compute_entropy_function(ghl_eos, rho[index], press[index]);
        if(!isfinite(value))
          CCTK_ERROR("Native entropy calculation produced a nonfinite proxy");
        entropy[index] = value;
      }
    }
  }
}

void GRHayLID_compute_entropy_tabulated(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLID_compute_entropy_tabulated;
  DECLARE_CCTK_PARAMETERS;
  GRHayLID_require_storage(cctkGH, "HydroBase::Y_e");
  GRHayLID_require_storage(cctkGH, "HydroBase::temperature");
  GRHayLID_require_storage(cctkGH, "HydroBase::entropy");

#pragma omp parallel for
  for(int k=0; k<cctk_lsh[2]; k++) {
    for(int j=0; j<cctk_lsh[1]; j++) {
      for(int i=0; i<cctk_lsh[0]; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        double density = rho[index], ye = Y_e[index], temp = temperature[index];
        if(!isfinite(density) || !isfinite(ye) || !isfinite(temp))
          CCTK_ERROR("Tabulated entropy requires finite initialized rho, Ye, and temperature");
        if(impose_beta_equilibrium)
          GRHayLID_check_table_state(density, ye, temp, "post-beta entropy");
        else
          ghl_tabulated_enforce_bounds_rho_Ye_T(ghl_eos, &density, &ye, &temp);
        GRHayLID_check_table_state(density, ye, temp, "entropy");
        double pressure, energy, value;
        const ghl_error_codes_t err = ghl_tabulated_compute_P_eps_S_from_T(
              ghl_eos, density, ye, temp, &pressure, &energy, &value);
        if(err != ghl_success || !isfinite(pressure) || !isfinite(energy) || !isfinite(value))
          CCTK_VERROR("Tabulated entropy EOS failed for rho=%g, Ye=%g, T=%g (status %d)",
                      density, ye, temp, (int)err);
        rho[index] = density;
        Y_e[index] = ye;
        temperature[index] = temp;
        press[index] = pressure;
        eps[index] = energy;
        entropy[index] = value;
      }
    }
  }
}
