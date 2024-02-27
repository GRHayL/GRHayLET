static inline int GRHayLMHD_local_avg(
    const cGH *restrict cctkGH, const int i, const int j, const int k,
    const int weight, const CCTK_INT *restrict needs_average,
    const CCTK_REAL *restrict rho_star, const CCTK_REAL *restrict tau,
    const CCTK_REAL *restrict Stildex, const CCTK_REAL *restrict Stildey,
    const CCTK_REAL *restrict Stildez, const CCTK_REAL *restrict ent_star,
    const CCTK_REAL *restrict Ye_star,
    ghl_conservative_quantities *restrict cons) {

  // This sets how many neighbors must exist to consider
  // the averaging scheme a success. 1 means we only
  // require at least 1 valid neighboring point to continue
  // with the averaging method.
  const int min_neighbors = 1;

  const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);
  const double wfac = weight / 4.0;

  const int iavg_min = MAX(0, i - 1);
  const int javg_min = MAX(0, j - 1);
  const int kavg_min = MAX(0, k - 1);
  const int iavg_max = MIN(cctkGH->cctk_lsh[0], i + 2);
  const int javg_max = MIN(cctkGH->cctk_lsh[1], j + 2);
  const int kavg_max = MIN(cctkGH->cctk_lsh[2], k + 2);
  int num_avg = 0;

  const double cfac = 1.0 - wfac;
  cons->rho = cfac * rho_star[index];
  cons->tau = cfac * tau[index];
  cons->SD[0] = cfac * Stildex[index];
  cons->SD[1] = cfac * Stildey[index];
  cons->SD[2] = cfac * Stildez[index];
  cons->entropy = cfac * ent_star[index];
  cons->Y_e = cfac * Ye_star[index];

  double rhotmp, tautmp, Sxtmp, Sytmp, Sztmp, enttmp, yetmp;
  rhotmp = tautmp = Sxtmp = Sytmp = Sztmp = enttmp = yetmp = 0;

  for (int kavg = kavg_min; kavg < kavg_max; kavg++) {
    for (int javg = javg_min; javg < javg_max; javg++) {
      for (int iavg = iavg_min; iavg < iavg_max; iavg++) {
        // We're only averaging over neighbors, so skip this point.
        // Also, we skip over any points that also failed.
        const int inavg = CCTK_GFINDEX3D(cctkGH, iavg, javg, kavg);
        if ((i == iavg && j == javg && k == kavg) || needs_average[inavg])
          continue;
        rhotmp += rho_star[inavg];
        tautmp += tau[inavg];
        Sxtmp += Stildex[inavg];
        Sytmp += Stildey[inavg];
        Sztmp += Stildez[inavg];
        enttmp += ent_star[inavg];
        yetmp += Ye_star[inavg];
        num_avg++;
      }
    }
  }

  if (num_avg < min_neighbors)
    return 1;

  cons->rho += wfac * rhotmp;
  cons->tau += wfac * tautmp;
  cons->SD[0] += wfac * Sxtmp;
  cons->SD[1] += wfac * Sytmp;
  cons->SD[2] += wfac * Sztmp;
  cons->entropy += wfac * enttmp;
  cons->Y_e += wfac * yetmp;

  // If weight=4, then the central point isn't part of the average
  num_avg += (weight != 4);

  cons->rho /= num_avg;
  cons->tau /= num_avg;
  cons->SD[0] /= num_avg;
  cons->SD[1] /= num_avg;
  cons->SD[2] /= num_avg;
  cons->entropy /= num_avg;
  cons->Y_e /= num_avg;

  return 0;
}

