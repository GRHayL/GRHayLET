#include "GRHayLMHD.h"

#define GRHAYLMHD_COMPUTE_DERIV(gf, m2, m1, p1, p2, inv_h) \
    (((1.0 / 12.0) * (gf[m2] - gf[p2]) + (2.0 / 3.0) * (gf[p1] - gf[m1])) * inv_h)

#define GRHAYLMHD_COMPUTE_METRIC_DERIVS(derivs_struct, inv_h, di, dj, dk)                    \
    {                                                                                        \
        const int m2 = CCTK_GFINDEX3D(cctkGH, i - 2 * di, j - 2 * dj, k - 2 * dk);           \
        const int m1 = CCTK_GFINDEX3D(cctkGH, i - 1 * di, j - 1 * dj, k - 1 * dk);           \
        const int p1 = CCTK_GFINDEX3D(cctkGH, i + 1 * di, j + 1 * dj, k + 1 * dk);           \
        const int p2 = CCTK_GFINDEX3D(cctkGH, i + 2 * di, j + 2 * dj, k + 2 * dk);           \
                                                                                             \
        const CCTK_REAL deriv_alp   = GRHAYLMHD_COMPUTE_DERIV(alp, m2, m1, p1, p2, inv_h);   \
        const CCTK_REAL deriv_betax = GRHAYLMHD_COMPUTE_DERIV(betax, m2, m1, p1, p2, inv_h); \
        const CCTK_REAL deriv_betay = GRHAYLMHD_COMPUTE_DERIV(betay, m2, m1, p1, p2, inv_h); \
        const CCTK_REAL deriv_betaz = GRHAYLMHD_COMPUTE_DERIV(betaz, m2, m1, p1, p2, inv_h); \
        const CCTK_REAL deriv_gxx   = GRHAYLMHD_COMPUTE_DERIV(gxx, m2, m1, p1, p2, inv_h);   \
        const CCTK_REAL deriv_gxy   = GRHAYLMHD_COMPUTE_DERIV(gxy, m2, m1, p1, p2, inv_h);   \
        const CCTK_REAL deriv_gxz   = GRHAYLMHD_COMPUTE_DERIV(gxz, m2, m1, p1, p2, inv_h);   \
        const CCTK_REAL deriv_gyy   = GRHAYLMHD_COMPUTE_DERIV(gyy, m2, m1, p1, p2, inv_h);   \
        const CCTK_REAL deriv_gyz   = GRHAYLMHD_COMPUTE_DERIV(gyz, m2, m1, p1, p2, inv_h);   \
        const CCTK_REAL deriv_gzz   = GRHAYLMHD_COMPUTE_DERIV(gzz, m2, m1, p1, p2, inv_h);   \
        ghl_initialize_metric(deriv_alp,                                                     \
                              deriv_betax,                                                   \
                              deriv_betay,                                                   \
                              deriv_betaz,                                                   \
                              deriv_gxx,                                                     \
                              deriv_gxy,                                                     \
                              deriv_gxz,                                                     \
                              deriv_gyy,                                                     \
                              deriv_gyz,                                                     \
                              deriv_gzz,                                                     \
                              &derivs_struct);                                               \
    }

void GRHayLMHD_Compute_RHSs_Sources(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    const CCTK_REAL inv_dx = 1.0 / CCTK_DELTA_SPACE(0);
    const CCTK_REAL inv_dy = 1.0 / CCTK_DELTA_SPACE(1);
    const CCTK_REAL inv_dz = 1.0 / CCTK_DELTA_SPACE(2);

    const int imin = cctk_nghostzones[0], imax = cctk_lsh[0] - cctk_nghostzones[0];
    const int jmin = cctk_nghostzones[1], jmax = cctk_lsh[1] - cctk_nghostzones[1];
    const int kmin = cctk_nghostzones[2], kmax = cctk_lsh[2] - cctk_nghostzones[2];

    OMPLOOP3D(imin, imax, jmin, jmax, kmin, kmax)
    {
        ghl_metric_quantities metric = { 0 };
        GRHAYLMHD_PACK_METRIC(metric);

        ghl_metric_quantities metric_derivs[3] = { { 0 }, { 0 }, { 0 } };
        GRHAYLMHD_COMPUTE_METRIC_DERIVS(metric_derivs[0], inv_dx, 1, 0, 0);
        GRHAYLMHD_COMPUTE_METRIC_DERIVS(metric_derivs[1], inv_dy, 0, 1, 0);
        GRHAYLMHD_COMPUTE_METRIC_DERIVS(metric_derivs[2], inv_dz, 0, 0, 1);

        ghl_primitive_quantities prims = { 0 };
        GRHAYLMHD_PACK_PRIMS(prims);

        ghl_extrinsic_curvature curv = { 0 };
        GRHAYLMHD_PACK_EXTRINSIC_CURVATURE(curv);

        ghl_conservative_quantities sources = { 0 };
        ghl_calculate_source_terms(ghl_eos,
                                   &prims,
                                   &metric,
                                   &metric_derivs[0],
                                   &metric_derivs[1],
                                   &metric_derivs[2],
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
