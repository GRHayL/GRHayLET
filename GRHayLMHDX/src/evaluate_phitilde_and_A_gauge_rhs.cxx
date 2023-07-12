#include "GRHayLMHDX.hxx"

extern "C" void GRHayLMHDX_evaluate_phitilde_and_A_gauge_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLMHDX_evaluate_phitilde_and_A_gauge_rhs;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dxi[3] = { 1.0/CCTK_DELTA_SPACE(0), 1.0/CCTK_DELTA_SPACE(1), 1.0/CCTK_DELTA_SPACE(2) };

  // \partial_t A_i = [reconstructed stuff] + [gauge stuff],
  //    where [gauge stuff] = -\partial_i (\alpha \Phi - \beta^j A_j)
  grid.loop_int_device<1, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    Ax_rhs(p.I) += dxi[0]*(A_rhs_stencil(p.I) - A_rhs_stencil(p.I + p.DI[0]));
if(Ax_rhs(p.I) > 100) CCTK_VINFO("Ax %e stencil %e %e", Ax_rhs(p.I), A_rhs_stencil(p.I), A_rhs_stencil(p.I + p.DI[0]));
  });

  grid.loop_int_device<0, 1, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    Ay_rhs(p.I) += dxi[1]*(A_rhs_stencil(p.I) - A_rhs_stencil(p.I + p.DI[1]));
if(Ay_rhs(p.I) > 100) CCTK_VINFO("Ay %e stencil %e %e", Ay_rhs(p.I), A_rhs_stencil(p.I), A_rhs_stencil(p.I + p.DI[1]));
  });

  grid.loop_int_device<0, 0, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    Az_rhs(p.I) += dxi[2]*(A_rhs_stencil(p.I) - A_rhs_stencil(p.I + p.DI[2]));
if(Az_rhs(p.I) > 100) CCTK_VINFO("Az %e stencil %e %e", Az_rhs(p.I), A_rhs_stencil(p.I), A_rhs_stencil(p.I + p.DI[2]));
  });

  const Loop::GF3D2layout layout(cctkGH, {0, 0, 0});

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index index(layout, p.I);

    double betax_stencil[5], betay_stencil[5], betaz_stencil[5];
    double phitilde_stencil[3][5], sqrtg_Ai_stencil[3][2];

    sqrtg_Ai_stencil[0][0] = sqrtg_Ax(index);
    sqrtg_Ai_stencil[1][0] = sqrtg_Ay(index);
    sqrtg_Ai_stencil[2][0] = sqrtg_Az(index);

    sqrtg_Ai_stencil[0][1] = sqrtg_Ax(p.I + p.DI[0]);
    sqrtg_Ai_stencil[1][1] = sqrtg_Ay(p.I + p.DI[1]);
    sqrtg_Ai_stencil[2][1] = sqrtg_Az(p.I + p.DI[2]);

    for(int iter=-2; iter<3; iter++) {
      const Loop::GF3D2index indexx(layout, p.I + iter*p.DI[0]);
      const Loop::GF3D2index indexy(layout, p.I + iter*p.DI[1]);
      const Loop::GF3D2index indexz(layout, p.I + iter*p.DI[2]);
      betax_stencil[iter+2] = betax(indexx);
      betay_stencil[iter+2] = betay(indexy);
      betaz_stencil[iter+2] = betaz(indexz);
      phitilde_stencil[0][iter+2] = phitilde(indexx);
      phitilde_stencil[1][iter+2] = phitilde(indexy);
      phitilde_stencil[2][iter+2] = phitilde(indexz);
    }
    phitilde_rhs(index) += ghl_calculate_phitilde_rhs(dxi, ghl_params->Lorenz_damping_factor, alp(index), betax_stencil, betay_stencil, betaz_stencil, sqrtg_Ai_stencil, phitilde_stencil);
  });
}
