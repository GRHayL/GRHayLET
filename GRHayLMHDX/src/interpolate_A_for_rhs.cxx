#include "GRHayLMHDX.hxx"

extern "C" void GRHayLMHDX_interpolate_A_for_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLMHDX_interpolate_A_for_rhs;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D2layout layout(cctkGH, {0, 0, 0});

  /* Compute \partial_t psi6phi = -\partial_i (  \alpha psi^6 A^i - psi6phi \beta^i)
   *    (Eq 13 of http://arxiv.org/pdf/1110.4633.pdf), using Lorenz gauge.
   * Note that the RHS consists of a shift advection term on psi6phi and
   *    a term depending on the vector potential.
   * psi6phi is defined at (i+1/2,j+1/2,k+1/2), but instead of reconstructing
   *    to compute the RHS of \partial_t psi6phi, we instead use standard
   *    interpolations.
   */
  const int loop_min[3] = {cctkGH->cctk_nghostzones[0],
                           cctkGH->cctk_nghostzones[1],
                           cctkGH->cctk_nghostzones[2]};
  const int loop_max[3] = {cctkGH->cctk_lsh[0] - cctkGH->cctk_nghostzones[0],
                           cctkGH->cctk_lsh[1] - cctkGH->cctk_nghostzones[1],
                           cctkGH->cctk_lsh[2] - cctkGH->cctk_nghostzones[2]};
  grid.loop_all_device<1, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      sqrtg_Ax(p.I) = 1e300;
  });
  grid.loop_all_device<0, 1, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      sqrtg_Ay(p.I) = 1e300;
  });
  grid.loop_all_device<0, 0, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      sqrtg_Az(p.I) = 1e300;
  });
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      A_rhs_stencil(p.I) = 1e300;
  });
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

    if((p.i >= loop_min[0] && p.i < loop_max[0]+2) &&
       (p.j >= loop_min[1] && p.j < loop_max[1]+2) &&
       (p.k >= loop_min[2] && p.k < loop_max[2]+2)) {
      // First compute \partial_j \alpha \sqrt{\gamma} A^j (RHS of \partial_i psi6phi)
      // FIXME: Would be much cheaper & easier to unstagger A_i, raise, then interpolate A^i.
      //        However, we keep it this way to be completely compatible with the original
      //        Illinois GRMHD thorn, called mhd_evolve.
      //
      //Step 1) j=x: Need to raise A_i, but to do that, we must have all variables at the same gridpoints:
      // The goal is to compute \partial_j (\alpha \sqrt{\gamma} A^j) at (i+1/2,j+1/2,k+1/2)
      //    We do this by first interpolating (RHS1x) = (\alpha \sqrt{\gamma} A^x) at
      //    (i,j+1/2,k+1/2)and (i+1,j+1/2,k+1/2), then taking \partial_x (RHS1x) =
      //    [ RHS1x(i+1,j+1/2,k+1/2) - RHS1x(i,j+1/2,k+1/2) ]/dX.
      // First bring gtup's, psi, and alpha to (i,j+1/2,k+1/2):
      ghl_metric_quantities metric_stencil[2][2][2];
      double Ax_stencil[3][3][3];
      double Ay_stencil[3][3][3];
      double Az_stencil[3][3][3];
      induction_interp_vars interp_vars;

      // Read in variable at interpolation stencil points from main memory.
      for(int iterz=0; iterz<2; iterz++) {
        for(int itery=0; itery<2; itery++) {
          for(int iterx=0; iterx<2; iterx++) {
            const Loop::GF3D2index ind(layout, p.I + iterx*p.DI[0] + itery*p.DI[1] + iterz*p.DI[2]);
            ghl_initialize_metric(
                  alp(ind), betax(ind), betay(ind), betaz(ind),
                  gxx(ind), gxy(ind), gxz(ind),
                  gyy(ind), gyz(ind), gzz(ind),
                  &metric_stencil[iterz][itery][iterx]);
          }
        }
      }
      // A_x needs a stencil s.t. interp_limits={ 0,1,-1,1,-1,1}.
      // A_y needs a stencil s.t. interp_limits={-1,1, 0,1,-1,1}.
      // A_z needs a stencil s.t. interp_limits={-1,1,-1,1, 0,1}.
      // We could fill only the needed elements, but it is cleaner
      // to fill in the whole 3x3x3 array.
      // TODO: the old code explicitly only filled in the necessary
      // elements. If we want to remove ~15 memcopies, do that here.
      for(int iterz=-1; iterz<2; iterz++) {
        for(int itery=-1; itery<2; itery++) {
          for(int iterx=-1; iterx<2; iterx++) {
            const auto ind = p.I + iterx*p.DI[0] + itery*p.DI[1] + iterz*p.DI[2];
            Ax_stencil[iterz+1][itery+1][iterx+1] = Ax(ind);
            Ay_stencil[iterz+1][itery+1][iterx+1] = Ay(ind);
            Az_stencil[iterz+1][itery+1][iterx+1] = Az(ind);
          }
        }
      }

      ghl_interpolate_with_vertex_centered_ADM(metric_stencil, Ax_stencil, Ay_stencil, Az_stencil, phitilde(p.I), &interp_vars);

      A_rhs_stencil(p.I) = isnan(interp_vars.alpha_Phi_minus_betaj_A_j) ? 1e300 : interp_vars.alpha_Phi_minus_betaj_A_j;
      sqrtg_Ax(p.I) = isnan(interp_vars.sqrtg_Ai[0]) ? 1e300 : interp_vars.sqrtg_Ai[0];
      sqrtg_Ay(p.I) = isnan(interp_vars.sqrtg_Ai[1]) ? 1e300 : interp_vars.sqrtg_Ai[1];
      sqrtg_Az(p.I) = isnan(interp_vars.sqrtg_Ai[2]) ? 1e300 : interp_vars.sqrtg_Ai[2];
    }
  });
}
