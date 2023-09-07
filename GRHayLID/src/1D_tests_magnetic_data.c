#include "GRHayLID.h"

// Note that we assume staggered vector potential
void GRHayLID_1D_tests_magnetic_data(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLID_1D_tests_magnetic_data;
  DECLARE_CCTK_PARAMETERS;

  if((!CCTK_EQUALS(initial_Avec, "GRHayLID")) && (!CCTK_EQUALS(initial_Bvec, "GRHayLID")))
    CCTK_VERROR("To use GRHayLID 1D initial magnetic data, please add either HydroBase::initial_Avec=\"GRHayLID\" or HydroBase::initial_Bvec=\"GRHayLID\" to the parfile.");

  double Bx_l, By_l, Bz_l;
  double Bx_r, By_r, Bz_r;
  if(CCTK_EQUALS(initial_data_1D,"Balsara1")) {
    Bx_l = Bx_r = 0.5;
    By_l = 1.0;
    By_r = -1.0;
    Bz_l = Bz_r = 0.0;
  } else if(CCTK_EQUALS(initial_data_1D,"Balsara2")) {
    Bx_l = Bx_r = 5.0;
    By_l = Bz_l = 6.0;
    By_r = Bz_r = 0.7;
  } else if(CCTK_EQUALS(initial_data_1D,"Balsara3")) {
    Bx_l = Bx_r = 10.0;
    By_l = Bz_l = 7.0;
    By_r = Bz_r = 0.7;
  } else if(CCTK_EQUALS(initial_data_1D,"Balsara4")) {
    Bx_l = Bx_r = 10.0;
    By_l = Bz_l = 7.0;
    By_r = Bz_r = -7.0;
  } else if(CCTK_EQUALS(initial_data_1D,"Balsara5")) {
    Bx_l = Bx_r = 2.0;
    By_l = Bz_l = 0.3;
    By_r = -0.7;
    Bz_r = 0.5;
  } else if(CCTK_EQUALS(initial_data_1D,"equilibrium")
         || CCTK_EQUALS(initial_data_1D,"sound wave")
         || CCTK_EQUALS(initial_data_1D,"shock tube")) {
    Bx_l = By_l = Bz_l = 0.0;
    Bx_r = By_r = Bz_r = 0.0;
  } else {
    CCTK_VERROR("Parameter initial_data_1D is not set "
                "to a valid test name. Something has gone wrong.");
  }

  //Above data assumes that the shock is in the x direction; we just have to rotate
  //the data for it to work in other directions.
  if(CCTK_EQUALS(shock_direction, "y")) {
    //x-->y, y-->z, z-->x
    double Bxtmp, Bytmp, Bztmp;
    Bxtmp = Bx_l; Bytmp = By_l; Bztmp = Bz_l;
    By_l = Bxtmp; Bz_l = Bytmp; Bx_l = Bztmp;

    Bxtmp = Bx_r; Bytmp = By_r; Bztmp = Bz_r;
    By_r = Bxtmp; Bz_r = Bytmp; Bx_r = Bztmp;
  } else if(CCTK_EQUALS(shock_direction, "z")) {
    //x-->z, y-->x, z-->y
    double Bxtmp, Bytmp, Bztmp;
    Bxtmp = Bx_l; Bytmp = By_l; Bztmp = Bz_l;
    Bz_l = Bxtmp; Bx_l = Bytmp; By_l = Bztmp;

    Bxtmp = Bx_r; Bytmp = By_r; Bztmp = Bz_r;
    Bz_r = Bxtmp; Bx_r = Bytmp; By_r = Bztmp;
  }

  CCTK_REAL dx[3] = { CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1), CCTK_DELTA_SPACE(2) };

#pragma omp parallel for
  for(int k=0; k<cctk_lsh[2]; k++) {
    for(int j=0; j<cctk_lsh[1]; j++) {
      for(int i=0; i<cctk_lsh[0]; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        const int ind4x = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0);
        const int ind4y = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1);
        const int ind4z = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2);

        double step = x[index];
        if(CCTK_EQUALS(shock_direction, "y")) {
          step = y[index];
        } else if(CCTK_EQUALS(shock_direction, "z")) {
          step = z[index];
        }

        if(step <= discontinuity_position) {
          Bvec[ind4x] = Bx_l;
          Bvec[ind4y] = By_l;
          Bvec[ind4z] = Bz_l;
        } else {
          Bvec[ind4x] = Bx_r;
          Bvec[ind4y] = By_r;
          Bvec[ind4z] = Bz_r;
        }
        const double x_stag = x[index] + stagger_A_fields*dx[0];
        const double y_stag = y[index] + stagger_A_fields*dx[1];
        const double z_stag = z[index] + stagger_A_fields*dx[2];

        step = x_stag;
        if(CCTK_EQUALS(shock_direction, "y")) {
          step = y_stag;
        } else if(CCTK_EQUALS(shock_direction, "z")) {
          step = z_stag;
        }

        if(step <= discontinuity_position) {
          if(CCTK_EQUALS(shock_direction, "x")) {
            Avec[ind4x] = By_l * z_stag - Bz_l * y_stag;
            Avec[ind4y] = 0.0;
            Avec[ind4z] = Bx_l * y_stag;
          } else if(CCTK_EQUALS(shock_direction, "y")) {
            Avec[ind4x] = By_l * z_stag;
            Avec[ind4y] = Bz_l * x_stag - Bx_l * z_stag;
            Avec[ind4z] = 0.0;
          } else {
            Avec[ind4x] = 0.0;
            Avec[ind4y] = Bz_l * x_stag;
            Avec[ind4z] = Bx_l * y_stag - By_l * x_stag;
          }
        } else {
          if(CCTK_EQUALS(shock_direction, "x")) {
            Avec[ind4x] = By_r * z_stag - Bz_r * y_stag;
            Avec[ind4y] = 0.0;
            Avec[ind4z] = Bx_r * y_stag;
          } else if(CCTK_EQUALS(shock_direction, "y")) {
            Avec[ind4x] = By_r * z_stag;
            Avec[ind4y] = Bz_r * x_stag - Bx_r * z_stag;
            Avec[ind4z] = 0.0;
          } else {
            Avec[ind4x] = 0.0;
            Avec[ind4y] = Bz_r * x_stag;
            Avec[ind4z] = Bx_r * y_stag - By_r * x_stag;
          }
        }
      }
    }
  }
  CCTK_VINFO("Finished writing magnetic ID for %s test", initial_data_1D);
}
