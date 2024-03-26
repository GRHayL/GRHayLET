#include "IllinoisGRMHDX.hxx"

/*
 * Compute \partial_t psi6phi = -\partial_i (  \alpha psi^6 A^i - psi6phi \beta^i)
 *    (Eq 13 of http://arxiv.org/pdf/1110.4633.pdf), using Lorenz gauge.
 * Note that the RHS consists of a shift advection term on psi6phi and
 *    a term depending on the vector potential.
 * psi6phi is defined at (i+1/2,j+1/2,k+1/2), but instead of reconstructing
 *    to compute the RHS of \partial_t psi6phi, we instead use standard
 *    interpolations.
*/
//CCTK_REAL A_rhs_stencil type = GF CENTERING={vvv} TAGS='prolongation="none" Checkpoint="no"' "interpolated quantity alpha Phi - beta^j A_j"
//CCTK_REAL sqrtg_Ax type = GF CENTERING={cvv} TAGS='prolongation="none" Checkpoint="no"' "interpolated quantity sqrt(g) A^x"
//CCTK_REAL sqrtg_Ay type = GF CENTERING={vcv} TAGS='prolongation="none" Checkpoint="no"' "interpolated quantity sqrt(g) A^y"
//CCTK_REAL sqrtg_Az type = GF CENTERING={vvc} TAGS='prolongation="none" Checkpoint="no"' "interpolated quantity sqrt(g) A^z"

extern "C" void IllinoisGRMHDX_evaluate_phitilde_and_A_gauge_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_IllinoisGRMHDX_evaluate_phitilde_and_A_gauge_rhs;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D2layout vvv_layout(cctkGH, {0, 0, 0});

  Arith::vect<int, Loop::dim> imin, imax;
  // Set up temporary variable with vvv layout
  grid.box_int<0, 0, 0>(grid.nghostzones, imin, imax);
  imin[0] -= 1; imin[1] -= 1; imin[2] -= 1;
  imax[0] += 1; imax[1] += 1; imax[2] += 1;
  const Loop::GF3D5layout local_vvv_layout(imin, imax);
  Loop::GF3D5vector<CCTK_REAL> tmps(local_vvv_layout, 1);
  int itmp = 0; const auto make_gf = [&]() { return Loop::GF3D5<CCTK_REAL>(tmps(itmp++)); };
  Loop::GF3D5<CCTK_REAL> A_rhs_stencil(make_gf());

  // Set up temporary variable with cvv layout including one +x ghost zone
  grid.box_int<1, 0, 0>(grid.nghostzones, imin, imax);
  imin[0] -= 1; imin[1] -= 1; imin[2] -= 1;
  imax[0] += 2; imax[1] += 1; imax[2] += 1;
  const Loop::GF3D5layout cvv_layout(imin, imax);
  Loop::GF3D5vector<CCTK_REAL> tmpsx(cvv_layout, 1);
  itmp = 0; const auto make_gfx = [&]() { return Loop::GF3D5<CCTK_REAL>(tmpsx(itmp++)); };
  Loop::GF3D5<CCTK_REAL> sqrtg_Ax(make_gfx());

  // Set up temporary variable with vcv layout including one +y ghost zone
  grid.box_int<0, 1, 0>(grid.nghostzones, imin, imax);
  imin[0] -= 1; imin[1] -= 1; imin[2] -= 1;
  imax[0] += 1; imax[1] += 2; imax[2] += 1;
  const Loop::GF3D5layout vcv_layout(imin, imax);
  Loop::GF3D5vector<CCTK_REAL> tmpsy(vcv_layout, 1);
  itmp = 0; const auto make_gfy = [&]() { return Loop::GF3D5<CCTK_REAL>(tmpsy(itmp++)); };
  Loop::GF3D5<CCTK_REAL> sqrtg_Ay(make_gfy());

  // Set up temporary variable with vvc layout including one +z ghost zone
  grid.box_int<0, 0, 1>(grid.nghostzones, imin, imax);
  imin[0] -= 1; imin[1] -= 1; imin[2] -= 1;
  imax[0] += 1; imax[1] += 1; imax[2] += 2;
  const Loop::GF3D5layout vvc_layout(imin, imax);
  Loop::GF3D5vector<CCTK_REAL> tmpsz(vvc_layout, 1);
  itmp = 0; const auto make_gfz = [&]() { return Loop::GF3D5<CCTK_REAL>(tmpsz(itmp++)); };
  Loop::GF3D5<CCTK_REAL> sqrtg_Az(make_gfz());

  Arith::vect<int, 3> bnd_min, bnd_max;
  grid.boundary_box<0, 0, 0>(grid.nghostzones, bnd_min, bnd_max);
  grid.box_int<0, 0, 0>(grid.nghostzones, imin, imax);

  // Extend the loop box. You must make sure that you extend by at most the
  // number of ghosts.
  for (int d = 0; d < 3; d++) {
    imin[d] -= 1;
    imax[d] += 1;
  }

  grid.loop_box_device<0, 0, 0>(
      bnd_min, bnd_max, imin, imax,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

    // First compute \partial_j \alpha \sqrt{\gamma} A^j (RHS of \partial_i psi6phi)
    // FIXME: Would be much cheaper & easier to unstagger A_i, raise, then interpolate A^i.
    //        However, we keep it this way to be completely compatible with the original
    //        Illinois GRMHD thorn, called mhd_evolve.

    // First bring gtup's, psi, and alpha to (i,j+1/2,k+1/2):
    ghl_metric_quantities metric_stencil[2][2][2];
    double Ax_stencil[3][3][3];
    double Ay_stencil[3][3][3];
    double Az_stencil[3][3][3];
    induction_interp_vars interp_vars;

    // Read in variable at interpolation stencil points from main memory.
    for(int iterz=-1; iterz<1; iterz++) {
      for(int itery=-1; itery<1; itery++) {
        for(int iterx=-1; iterx<1; iterx++) {
          const Loop::GF3D2index ind(vvv_layout, p.I + iterx*p.DI[0] + itery*p.DI[1] + iterz*p.DI[2]);
          ghl_initialize_metric(
                alp(ind), betax(ind), betay(ind), betaz(ind),
                gxx(ind), gxy(ind), gxz(ind),
                gyy(ind), gyz(ind), gzz(ind),
                &metric_stencil[iterz+1][itery+1][iterx+1]);
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
    for(int iterz=-2; iterz<1; iterz++) {
      for(int itery=-2; itery<1; itery++) {
        for(int iterx=-2; iterx<1; iterx++) {
          const auto ind = p.I + iterx*p.DI[0] + itery*p.DI[1] + iterz*p.DI[2];
          Ax_stencil[iterz+2][itery+2][iterx+2] = Ax(ind);
          Ay_stencil[iterz+2][itery+2][iterx+2] = Ay(ind);
          Az_stencil[iterz+2][itery+2][iterx+2] = Az(ind);
        }
      }
    }

    ghl_interpolate_with_vertex_centered_ADM(metric_stencil, Ax_stencil, Ay_stencil, Az_stencil, phitilde(p.I), &interp_vars);

    const Loop::GF3D5index vvv_ind(local_vvv_layout, p.I);
    const Loop::GF3D5index cvv_ind(cvv_layout, p.I);
    const Loop::GF3D5index vcv_ind(vcv_layout, p.I);
    const Loop::GF3D5index vvc_ind(vvc_layout, p.I);
    A_rhs_stencil(vvv_ind) = isnan(interp_vars.alpha_Phi_minus_betaj_A_j) ? 1e300 : interp_vars.alpha_Phi_minus_betaj_A_j;
    sqrtg_Ax(cvv_ind) = isnan(interp_vars.sqrtg_Ai[0]) ? 1e300 : interp_vars.sqrtg_Ai[0];
    sqrtg_Ay(vcv_ind) = isnan(interp_vars.sqrtg_Ai[1]) ? 1e300 : interp_vars.sqrtg_Ai[1];
    sqrtg_Az(vvc_ind) = isnan(interp_vars.sqrtg_Ai[2]) ? 1e300 : interp_vars.sqrtg_Ai[2];
  });

  const CCTK_REAL dxi[3] = { 1.0/CCTK_DELTA_SPACE(0), 1.0/CCTK_DELTA_SPACE(1), 1.0/CCTK_DELTA_SPACE(2) };

  // \partial_t A_i = [reconstructed stuff] + [gauge stuff],
  //    where [gauge stuff] = -\partial_i (\alpha \Phi - \beta^j A_j)
//  grid.loop_int_device<1, 0, 0>(
//      grid.nghostzones,
//      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
//    const Loop::GF3D5index vvv_index(local_vvv_layout, p.I);
//    const Loop::GF3D5index vvv_indp1(local_vvv_layout, p.I + p.DI[0]);
//    Ax_rhs(p.I) += dxi[0]*(A_rhs_stencil(vvv_index) - A_rhs_stencil(vvv_indp1));
//  });
//
//  grid.loop_int_device<0, 1, 0>(
//      grid.nghostzones,
//      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
//    const Loop::GF3D5index vvv_index(local_vvv_layout, p.I);
//    const Loop::GF3D5index vvv_indp1(local_vvv_layout, p.I + p.DI[1]);
//    Ay_rhs(p.I) += dxi[1]*(A_rhs_stencil(vvv_index) - A_rhs_stencil(vvv_indp1));
//  });
//
//  grid.loop_int_device<0, 0, 1>(
//      grid.nghostzones,
//      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
//    const Loop::GF3D5index vvv_index(local_vvv_layout, p.I);
//    const Loop::GF3D5index vvv_indp1(local_vvv_layout, p.I + p.DI[2]);
//    Az_rhs(p.I) += dxi[2]*(A_rhs_stencil(vvv_index) - A_rhs_stencil(vvv_indp1));
//  });
//
  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index index(vvv_layout, p.I);
//
//    double betax_stencil[5], betay_stencil[5], betaz_stencil[5];
//    double phitilde_stencil[3][5], sqrtg_Ai_stencil[3][2];
//
//    const Loop::GF3D5index cvv_indm1(cvv_layout, p.I - p.DI[0]);
//    const Loop::GF3D5index vcv_indm1(vcv_layout, p.I - p.DI[1]);
//    const Loop::GF3D5index vvc_indm1(vvc_layout, p.I - p.DI[2]);
//    sqrtg_Ai_stencil[0][0] = sqrtg_Ax(cvv_indm1);
//    sqrtg_Ai_stencil[1][0] = sqrtg_Ay(vcv_indm1);
//    sqrtg_Ai_stencil[2][0] = sqrtg_Az(vvc_indm1);
//
//    const Loop::GF3D5index cvv_ind(cvv_layout, p.I);
//    const Loop::GF3D5index vcv_ind(vcv_layout, p.I);
//    const Loop::GF3D5index vvc_ind(vvc_layout, p.I);
//    sqrtg_Ai_stencil[0][1] = sqrtg_Ax(cvv_ind);
//    sqrtg_Ai_stencil[1][1] = sqrtg_Ay(vcv_ind);
//    sqrtg_Ai_stencil[2][1] = sqrtg_Az(vvc_ind);
//
//    for(int iter=-2; iter<3; iter++) {
//      const Loop::GF3D2index indexx(vvv_layout, p.I + iter*p.DI[0]);
//      const Loop::GF3D2index indexy(vvv_layout, p.I + iter*p.DI[1]);
//      const Loop::GF3D2index indexz(vvv_layout, p.I + iter*p.DI[2]);
//      betax_stencil[iter+2] = betax(indexx);
//      betay_stencil[iter+2] = betay(indexy);
//      betaz_stencil[iter+2] = betaz(indexz);
//      phitilde_stencil[0][iter+2] = phitilde(indexx);
//      phitilde_stencil[1][iter+2] = phitilde(indexy);
//      phitilde_stencil[2][iter+2] = phitilde(indexz);
//    }
//    phitilde_rhs(index) = ghl_calculate_phitilde_rhs(dxi, ghl_params->Lorenz_damping_factor, alp(index), betax_stencil, betay_stencil, betaz_stencil, sqrtg_Ai_stencil, phitilde_stencil);
    phitilde_rhs(index) = 0.0;
  });
}
