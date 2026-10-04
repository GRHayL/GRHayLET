#include "GRHayLID.h"

/*
 *
 * (c) 2021 Leo Werneck
 *
 * This is the thorn's driver function, responsible
 * for setting the initial data to that of an constant
 * density sphere in Minkowski space.
 */
void GRHayLID_ConstantDensitySphere(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLID_ConstantDensitySphere;
  DECLARE_CCTK_PARAMETERS;

  if(!CCTK_EQUALS(EOS_type, "Tabulated"))
    CCTK_ERROR("ConstantDensitySphere initial data is only defined for tabulated EOS. Please change GRHayLib::EOS_type to \"Tabulated\" in the parfile.");
  if(!CCTK_EQUALS(initial_Y_e, "GRHayLID"))
    CCTK_ERROR("To use ConstantDensitySphere initial data, please add initial_Y_e=\"GRHayLID\" to the parfile.");
  if(!CCTK_EQUALS(initial_temperature, "GRHayLID"))
    CCTK_ERROR("To use ConstantDensitySphere initial data, please add initial_temperature=\"GRHayLID\" to the parfile.");

  CHECK_PARAMETER(ConstantDensitySphere_sphere_radius);
  CHECK_PARAMETER(ConstantDensitySphere_rho_interior);
  CHECK_PARAMETER(ConstantDensitySphere_Y_e_interior);
  CHECK_PARAMETER(ConstantDensitySphere_T_interior);
  CHECK_PARAMETER(ConstantDensitySphere_rho_exterior);
  CHECK_PARAMETER(ConstantDensitySphere_Y_e_exterior);
  CHECK_PARAMETER(ConstantDensitySphere_T_exterior);

  GRHayLID_require_storage(cctkGH, "HydroBase::Y_e");
  GRHayLID_require_storage(cctkGH, "HydroBase::temperature");
  GRHayLID_require_storage(cctkGH, "ADMBase::metric");

  CCTK_INFO("Beginning ConstantDensitySphere initial data");

  // Compute hydro quantities inside and outside the sphere
  double P_interior, eps_interior;
  GRHayLID_check_table_state(ConstantDensitySphere_rho_interior, ConstantDensitySphere_Y_e_interior, ConstantDensitySphere_T_interior, "sphere interior");
  ghl_error_codes_t eos_error = ghl_tabulated_compute_P_eps_from_T(
        ghl_eos,
        ConstantDensitySphere_rho_interior,
        ConstantDensitySphere_Y_e_interior,
        ConstantDensitySphere_T_interior,
        &P_interior,
        &eps_interior);
  if(eos_error != ghl_success || !isfinite(P_interior) || !isfinite(eps_interior))
    CCTK_VERROR("sphere interior EOS failed for rho=%g, Ye=%g, T=%g (status %d)",
                ConstantDensitySphere_rho_interior, ConstantDensitySphere_Y_e_interior, ConstantDensitySphere_T_interior, (int)eos_error);

  double P_exterior, eps_exterior;
  GRHayLID_check_table_state(ConstantDensitySphere_rho_exterior, ConstantDensitySphere_Y_e_exterior, ConstantDensitySphere_T_exterior, "sphere exterior");
  eos_error = ghl_tabulated_compute_P_eps_from_T(
        ghl_eos,
        ConstantDensitySphere_rho_exterior,
        ConstantDensitySphere_Y_e_exterior,
        ConstantDensitySphere_T_exterior,
        &P_exterior,
        &eps_exterior);
  if(eos_error != ghl_success || !isfinite(P_exterior) || !isfinite(eps_exterior))
    CCTK_VERROR("sphere exterior EOS failed for rho=%g, Ye=%g, T=%g (status %d)",
                ConstantDensitySphere_rho_exterior, ConstantDensitySphere_Y_e_exterior,
                ConstantDensitySphere_T_exterior, (int)eos_error);

  const double vx = ConstantDensitySphere_vx_interior;
  const double vy = ConstantDensitySphere_vy_interior;
  const double vz = ConstantDensitySphere_vz_interior;
  if(!isfinite(vx) || !isfinite(vy) || !isfinite(vz) || vx*vx+vy*vy+vz*vz >= 1.0)
    CCTK_ERROR("Sphere interior Cartesian velocity must be finite and strictly subluminal");
  if(!isfinite(ConstantDensitySphere_sphere_radius) || ConstantDensitySphere_sphere_radius < 0.0)
    CCTK_ERROR("Sphere radius must be finite and nonnegative");

#pragma omp parallel for
  for(int k=0; k<cctk_lsh[2]; k++) {
    for(int j=0; j<cctk_lsh[1]; j++) {
      for(int i=0; i<cctk_lsh[0]; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        GRHayLID_check_flat_metric(gxx[index], gxy[index], gxz[index], gyy[index], gyz[index], gzz[index]);
        const int ind4x = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0);
        const int ind4y = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1);
        const int ind4z = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2);

        vel[ind4x] = 0;
        vel[ind4y] = 0;
        vel[ind4z] = 0;
        if( r[index] > ConstantDensitySphere_sphere_radius ) {
          // Outside the sphere
          rho        [index] = ConstantDensitySphere_rho_exterior;
          Y_e        [index] = ConstantDensitySphere_Y_e_exterior;
          temperature[index] = ConstantDensitySphere_T_exterior;
          vel        [ind4x] = 0;
          vel        [ind4y] = 0;
          vel        [ind4z] = 0;
          press      [index] = P_exterior;
          eps        [index] = eps_exterior;
        }
        else {
          // Inside the sphere
          rho        [index] = ConstantDensitySphere_rho_interior;
          Y_e        [index] = ConstantDensitySphere_Y_e_interior;
          temperature[index] = ConstantDensitySphere_T_interior;
          vel        [ind4x] = ConstantDensitySphere_vx_interior;
          vel        [ind4y] = ConstantDensitySphere_vy_interior;
          vel        [ind4z] = ConstantDensitySphere_vz_interior;
          press      [index] = P_interior;
          eps        [index] = eps_interior;
        }
      }
    }
  }

  CCTK_INFO("All done!");
}
