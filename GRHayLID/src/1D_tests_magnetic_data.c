#include "GRHayLID.h"

// Cartesian vector potential, with optional component-specific half-cell shifts.
void GRHayLID_1D_tests_magnetic_data(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLID_1D_tests_magnetic_data;
  DECLARE_CCTK_PARAMETERS;

  if(!CCTK_EQUALS(initial_Avec, "GRHayLID") || !CCTK_EQUALS(initial_Bvec, "GRHayLID"))
    CCTK_ERROR("GRHayLID magnetic initialization requires both initial_Avec and initial_Bvec=GRHayLID");
  GRHayLID_require_storage(cctkGH, "HydroBase::Avec");
  GRHayLID_require_storage(cctkGH, "HydroBase::Bvec");
  GRHayLID_require_storage(cctkGH, "ADMBase::metric");

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

        if(step <= discontinuity_position) {
          Bvec[ind4x] = Bx_l;
          Bvec[ind4y] = By_l;
          Bvec[ind4z] = Bz_l;
        } else {
          Bvec[ind4x] = Bx_r;
          Bvec[ind4y] = By_r;
          Bvec[ind4z] = Bz_r;
        }
        // A_x: (x,y+dy/2,z+dz/2), A_y: (x+dx/2,y,z+dz/2),
        // A_z: (x+dx/2,y+dy/2,z). Bvec remains at the base coordinates.
        const int axis = CCTK_EQUALS(shock_direction, "y") ? 1 :
                         CCTK_EQUALS(shock_direction, "z") ? 2 : 0;
        const double left[3] = {Bx_l, By_l, Bz_l};
        const double right[3] = {Bx_r, By_r, Bz_r};
        for(int component=0; component<3; component++) {
          const double x_stag = x[index] + (stagger_A_fields && component != 0 ? 0.5*CCTK_DELTA_SPACE(0) : 0.0);
          const double y_stag = y[index] + (stagger_A_fields && component != 1 ? 0.5*CCTK_DELTA_SPACE(1) : 0.0);
          const double z_stag = z[index] + (stagger_A_fields && component != 2 ? 0.5*CCTK_DELTA_SPACE(2) : 0.0);
          const double position[3] = {x_stag, y_stag, z_stag};
          const double *field = position[axis] <= discontinuity_position ? left : right;
          double value = 0.0;
          if(axis == 0) {
            if(component == 0) value = field[1]*z_stag - field[2]*y_stag;
            if(component == 2) value = field[0]*y_stag;
          } else if(axis == 1) {
            if(component == 0) value = field[1]*z_stag;
            if(component == 1) value = field[2]*x_stag - field[0]*z_stag;
          } else {
            if(component == 1) value = field[2]*x_stag;
            if(component == 2) value = field[0]*y_stag - field[1]*x_stag;
          }
          const int destination = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,component);
          Avec[destination] = value;
        }
      }
    }
  }
  CCTK_VINFO("Finished writing magnetic ID for %s test", initial_data_1D);
}
