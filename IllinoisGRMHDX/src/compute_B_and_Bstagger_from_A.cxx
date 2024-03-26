#include "GRHayLMHDX.hxx"

/*
   Set e.g. Bx_stagger = \partial_y A_z - partial_z A_y
   Ax has cvv staggering
   Ay has vcv staggering
   Az has vvc staggering
   while Bx_stagger has vcc staggering
   Therefore, the 2nd order derivative \partial_z A_y with
   the v->c conversion is:
          [Ay(i,j,k+1) - Ay(i,j,k)]/dZ
*/

extern "C" void GRHayLMHDX_compute_B_and_Bstagger_from_A(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLMHDX_compute_B_and_Bstagger_from_A;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D2layout layout(cctkGH, {1, 1, 1});

  const CCTK_REAL dxi = 1.0/CCTK_DELTA_SPACE(0);
  const CCTK_REAL dyi = 1.0/CCTK_DELTA_SPACE(1);
  const CCTK_REAL dzi = 1.0/CCTK_DELTA_SPACE(2);

  /**************/
  /* Bx_stagger */
  /**************/
  grid.loop_all_device<0, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const auto jplus = p.I + p.DI[1];
    const auto kplus = p.I + p.DI[2];
    Bx_stagger(p.I) = dyi*(Az(jplus) - Az(p.I)) - dzi*(Ay(kplus) - Ay(p.I));
if(Bx_stagger(p.I) > 50) CCTK_VINFO("B %e Az %e %e Ay %e %e", Bx_stagger(p.I), Az(jplus), Az(p.I), Ay(kplus), Ay(p.I));
  }); // vcc loop interior

  /**************/
  /* By_stagger */
  /**************/
  grid.loop_all_device<1, 0, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const auto iplus = p.I + p.DI[0];
    const auto kplus = p.I + p.DI[2];
    By_stagger(p.I) = dzi*(Ax(kplus) - Ax(p.I)) - dxi*(Az(iplus) - Az(p.I));
  }); // cvc loop interior

  /**************/
  /* Bz_stagger */
  /**************/
  grid.loop_all_device<1, 1, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const auto iplus = p.I + p.DI[0];
    const auto jplus = p.I + p.DI[1];
    Bz_stagger(p.I) = dxi*(Ay(iplus) - Ay(p.I)) - dyi*(Ax(jplus) - Ax(p.I));
  }); // ccv loop interior

  /*
     To compute cell-centered B, we average the staggered \tilde{B}:
        Bx[x,y,z] = 1/2 * ( Bx_stag[x+dx/2,y,z] + Bx_stag[x-dx/2,y,z] )
     However, since the staggered variable is densitized, we also divide
     by \sqrt{\gamma}
  */
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index index(layout, p.I);

    Bx_center(index) = 0.01; //0.5 * ( Bx_stagger(p.I + p.DI[0]) + Bx_stagger(p.I) )/sqrt_detgamma(index);
    By_center(index) = 0.01; //0.5 * ( By_stagger(p.I + p.DI[1]) + By_stagger(p.I) )/sqrt_detgamma(index);
    Bz_center(index) = 0.01; //0.5 * ( Bz_stagger(p.I + p.DI[2]) + Bz_stagger(p.I) )/sqrt_detgamma(index);
  }); // ccc loop interior
}
