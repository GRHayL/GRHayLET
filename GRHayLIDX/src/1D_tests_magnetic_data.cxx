#include "GRHayLIDX.h"

extern "C" void GRHayLIDX_1D_tests_magnetic_data(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLIDX_1D_tests_magnetic_data;
  DECLARE_CCTK_PARAMETERS;

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

  const Loop::GF3D2layout layout(cctkGH, {1, 1, 1});

  //This would probably be cleaner if CCTK_EQUALS was available on device

  if(CCTK_EQUALS(shock_direction, "x")) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      const Loop::GF3D2index index(layout, p.I);
      if(p.x <= discontinuity_position) {
        Bvecx(index) = Bx_l;
        Bvecy(index) = By_l;
        Bvecz(index) = Bz_l;
      } else {
        Bvecx(index) = Bx_r;
        Bvecy(index) = By_r;
        Bvecz(index) = Bz_r;
      }
    });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      if(p.x <= discontinuity_position) {
        Avecx(p.I) = By_l * (p.z) - Bz_l * (p.y);
      } else {
        Avecx(p.I) = By_r * (p.z) - Bz_r * (p.y);
      }
    });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      Avecy(p.I) = 0.0;
    });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      Avecz(p.I) = Bx_r * (p.y);
    });

  } else if(CCTK_EQUALS(shock_direction, "y")) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      const Loop::GF3D2index index(layout, p.I);

      if(p.y <= discontinuity_position) {
        Bvecx(index) = Bx_l;
        Bvecy(index) = By_l;
        Bvecz(index) = Bz_l;
      } else {
        Bvecx(index) = Bx_r;
        Bvecy(index) = By_r;
        Bvecz(index) = Bz_r;
      }
    });
  
    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      Avecx(p.I) = By_r * (p.z);
    });
  
    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
  
      if(p.y <= discontinuity_position) {
        Avecy(p.I) = Bz_l * (p.x) - Bx_l * (p.z);
      } else {
        Avecy(p.I) = Bz_r * (p.x) - Bx_r * (p.z);
      }
    });
  
    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      Avecz(p.I) = 0.0;
    });

  } else {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      const Loop::GF3D2index index(layout, p.I);

      if(p.z <= discontinuity_position) {
        Bvecx(index) = Bx_l;
        Bvecy(index) = By_l;
        Bvecz(index) = Bz_l;
      } else {
        Bvecx(index) = Bx_r;
        Bvecy(index) = By_r;
        Bvecz(index) = Bz_r;
      }
    });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      Avecx(p.I) = 0.0;
    });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      Avecy(p.I) = Bz_r * (p.x);
    });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

      if(p.z <= discontinuity_position) {
        Avecz(p.I) = Bx_l * (p.y) - By_l * (p.x);
      } else {
        Avecz(p.I) = Bx_r * (p.y) - By_r * (p.x);
      }
    });
  }
}
