#include "GRHayLMHD.h"

static void compute_metric_derivs(const CCTK_REAL inv_ds,
                                  const ghl_metric_quantities *restrict metric,
                                  const ghl_metric_quantities *restrict metric_p1,
                                  ghl_metric_quantities *restrict derivs)
{
    derivs->lapse         = (metric_p1->lapse - metric->lapse) * inv_ds;
    derivs->betaU[0]      = (metric_p1->betaU[0] - metric->betaU[0]) * inv_ds;
    derivs->betaU[1]      = (metric_p1->betaU[1] - metric->betaU[1]) * inv_ds;
    derivs->betaU[2]      = (metric_p1->betaU[2] - metric->betaU[2]) * inv_ds;
    derivs->gammaDD[0][0] = (metric_p1->gammaDD[0][0] - metric->gammaDD[0][0]) * inv_ds;
    derivs->gammaDD[0][1] = (metric_p1->gammaDD[0][1] - metric->gammaDD[0][1]) * inv_ds;
    derivs->gammaDD[0][2] = (metric_p1->gammaDD[0][2] - metric->gammaDD[0][2]) * inv_ds;
    derivs->gammaDD[1][1] = (metric_p1->gammaDD[1][1] - metric->gammaDD[1][1]) * inv_ds;
    derivs->gammaDD[1][2] = (metric_p1->gammaDD[1][2] - metric->gammaDD[1][2]) * inv_ds;
    derivs->gammaDD[2][2] = (metric_p1->gammaDD[2][2] - metric->gammaDD[2][2]) * inv_ds;
}

void GRHayLMHD_Compute_RHSs_Sources(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    CCTK_REAL inv_dx = 1.0 / CCTK_DELTA_SPACE(0);
    CCTK_REAL inv_dy = 1.0 / CCTK_DELTA_SPACE(1);
    CCTK_REAL inv_dz = 1.0 / CCTK_DELTA_SPACE(2);

    const int imin = cctk_nghostzones[0], imax = cctk_lsh[0] - cctk_nghostzones[0];
    const int jmin = cctk_nghostzones[1], jmax = cctk_lsh[1] - cctk_nghostzones[1];
    const int kmin = cctk_nghostzones[2], kmax = cctk_lsh[2] - cctk_nghostzones[2];

    OMPLOOP3D(imin, imax, jmin, jmax, kmin, kmax)
    {
        const int ip1jk = CCTK_GFINDEX3D(cctkGH, i + 1, j, k);
        const int ijp1k = CCTK_GFINDEX3D(cctkGH, i, j + 1, k);
        const int ijkp1 = CCTK_GFINDEX3D(cctkGH, i, j, k + 1);

        ghl_metric_quantities metric_ijk = { 0 };
        ghl_initialize_metric(alp[ijk],
                              betax[ijk],
                              betay[ijk],
                              betaz[ijk],
                              gxx[ijk],
                              gxy[ijk],
                              gxz[ijk],
                              gyy[ijk],
                              gyz[ijk],
                              gzz[ijk],
                              &metric_ijk);

        ghl_metric_quantities metric_ip1jk = { 0 };
        ghl_initialize_metric(alp[ip1jk],
                              betax[ip1jk],
                              betay[ip1jk],
                              betaz[ip1jk],
                              gxx[ip1jk],
                              gxy[ip1jk],
                              gxz[ip1jk],
                              gyy[ip1jk],
                              gyz[ip1jk],
                              gzz[ip1jk],
                              &metric_ip1jk);

        ghl_metric_quantities metric_ijp1k = { 0 };
        ghl_initialize_metric(alp[ijp1k],
                              betax[ijp1k],
                              betay[ijp1k],
                              betaz[ijp1k],
                              gxx[ijp1k],
                              gxy[ijp1k],
                              gxz[ijp1k],
                              gyy[ijp1k],
                              gyz[ijp1k],
                              gzz[ijp1k],
                              &metric_ijp1k);

        ghl_metric_quantities metric_ijkp1 = { 0 };
        ghl_initialize_metric(alp[ijkp1],
                              betax[ijkp1],
                              betay[ijkp1],
                              betaz[ijkp1],
                              gxx[ijkp1],
                              gxy[ijkp1],
                              gxz[ijkp1],
                              gyy[ijkp1],
                              gyz[ijkp1],
                              gzz[ijkp1],
                              &metric_ijkp1);

        ghl_metric_quantities metric_derivs_x = { 0 };
        compute_metric_derivs(inv_dx, &metric_ijk, &metric_ip1jk, &metric_derivs_x);

        ghl_metric_quantities metric_derivs_y = { 0 };
        compute_metric_derivs(inv_dy, &metric_ijk, &metric_ijp1k, &metric_derivs_y);

        ghl_metric_quantities metric_derivs_z = { 0 };
        compute_metric_derivs(inv_dz, &metric_ijk, &metric_ijkp1, &metric_derivs_z);

        ghl_primitive_quantities prims = { 0 };
        ghl_initialize_primitives(rho[ijk],
                                  press[ijk],
                                  eps[ijk],
                                  vx[ijk],
                                  vy[ijk],
                                  vz[ijk],
                                  Bvec[ijkx],
                                  Bvec[ijky],
                                  Bvec[ijkz],
                                  entropy[ijk],
                                  Y_e[ijk],
                                  temperature[ijk],
                                  &prims);

        ghl_extrinsic_curvature curv = { 0 };
        ghl_initialize_extrinsic_curvature(kxx[ijk],
                                           kxy[ijk],
                                           kxz[ijk],
                                           kyy[ijk],
                                           kyz[ijk],
                                           kzz[ijk],
                                           &curv);

        ghl_conservative_quantities sources = { 0 };
        ghl_calculate_source_terms(ghl_eos,
                                   &prims,
                                   &metric_ijk,
                                   &metric_derivs_x,
                                   &metric_derivs_y,
                                   &metric_derivs_z,
                                   &curv,
                                   &sources);

        rho_rhs[ijk] = sources.rho;
        tau_rhs[ijk] = sources.tau;
        S_x_rhs[ijk] = sources.SD[0];
        S_y_rhs[ijk] = sources.SD[1];
        S_z_rhs[ijk] = sources.SD[2];
        Y_e_rhs[ijk] = sources.Y_e;
        ent_rhs[ijk] = sources.entropy;
    }
    ENDLOOP3D
}
