#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "GRHayLib.h"

#define velx (&vel[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely (&vel[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz (&vel[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])

#define CHECK_PARAMETER(par) if(par==-1) CCTK_VERROR("Please set %s::%s in your parfile",CCTK_THORNSTRING,#par);

/*
 *
 * (c) 2021 Leo Werneck
 *
 * This is the thorn's driver function, responsible
 * for setting the initial data to that of an constant
 * density sphere in Minkowski space.
 */
void ConstantDensitySphereID(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ConstantDensitySphereID;
  DECLARE_CCTK_PARAMETERS;

  if(!CCTK_EQUALS(initial_hydro, "ConstantDensitySphereID"))
    CCTK_VERROR("To use ConstantDensitySphereID, please add initial_hydro=\"ConstantDensitySphereID\" to the parfile.");
  if(!CCTK_EQUALS(initial_Y_e, "ConstantDensitySphereID"))
    CCTK_VERROR("To use ConstantDensitySphereID, please add initial_Y_e=\"ConstantDensitySphereID\" to the parfile.");
  if(!CCTK_EQUALS(initial_temperature, "ConstantDensitySphereID"))
    CCTK_VERROR("To use ConstantDensitySphereID, please add initial_temperature=\"ConstantDensitySphereID\" to the parfile.");

  CHECK_PARAMETER(ConstantDensitySphereID_sphere_radius);
  CHECK_PARAMETER(ConstantDensitySphereID_rho_interior);
  CHECK_PARAMETER(ConstantDensitySphereID_Y_e_interior);
  CHECK_PARAMETER(ConstantDensitySphereID_T_interior);
  CHECK_PARAMETER(ConstantDensitySphereID_rho_exterior);
  CHECK_PARAMETER(ConstantDensitySphereID_Y_e_exterior);
  CHECK_PARAMETER(ConstantDensitySphereID_T_exterior);

  CCTK_INFO("Beginning ConstantDensitySphereID initial data");

  // Compute hydro quantities inside and outside the sphere
  CCTK_REAL P_interior, eps_interior, S_interior;
  ghl_tabulated_compute_P_eps_S_from_T(
        ghl_eos,
        ConstantDensitySphereID_rho_interior,
        ConstantDensitySphereID_Y_e_interior,
        ConstantDensitySphereID_T_interior,
        &P_interior,
        &eps_interior,
        &S_interior );

  CCTK_REAL P_exterior, eps_exterior, S_exterior;
  ghl_tabulated_compute_P_eps_S_from_T(
        ghl_eos,
        ConstantDensitySphereID_rho_exterior,
        ConstantDensitySphereID_Y_e_exterior,
        ConstantDensitySphereID_T_exterior,
        &P_exterior,
        &eps_exterior,
        &S_exterior );

#pragma omp parallel for
  for(int k=0; k<cctk_lsh[2]; k++) {
    for(int j=0; j<cctk_lsh[1]; j++) {
      for(int i=0; i<cctk_lsh[0]; i++) {
        const int idx = CCTK_GFINDEX3D(cctkGH,i,j,k);

        velx[idx] = 0;
        vely[idx] = 0;
        velz[idx] = 0;
        if( r[idx] > ConstantDensitySphereID_sphere_radius ) {
          // Outside the sphere
          rho        [idx] = ConstantDensitySphereID_rho_exterior;
          Y_e        [idx] = ConstantDensitySphereID_Y_e_exterior;
          temperature[idx] = ConstantDensitySphereID_T_exterior;
          press      [idx] = P_exterior;
          eps        [idx] = eps_exterior;
          if(CCTK_EQUALS(initial_entropy, "ConstantDensitySphereID")) entropy[idx] = S_exterior;
        }
        else {
          // Inside the sphere
          rho        [idx] = ConstantDensitySphereID_rho_interior;
          Y_e        [idx] = ConstantDensitySphereID_Y_e_interior;
          temperature[idx] = ConstantDensitySphereID_T_interior;
          press      [idx] = P_interior;
          eps        [idx] = eps_interior;
          if(CCTK_EQUALS(initial_entropy, "ConstantDensitySphereID")) entropy[idx] = S_interior;
        }
      }
    }
  }

  CCTK_INFO("All done!");
}
