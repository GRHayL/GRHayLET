#include "GRHayLMHD.h"

typedef void (*ghl_cmin_cmax_func_t)(ghl_primitive_quantities *restrict prims_r,
                                     ghl_primitive_quantities *restrict prims_l,
                                     const ghl_eos_parameters *restrict eos,
                                     const ghl_metric_quantities *restrict ADM_metric_face,
                                     double *cmin,
                                     double *cmax);

typedef void (*ghl_flux_func_t)(ghl_primitive_quantities *restrict prims_r,
                                ghl_primitive_quantities *restrict prims_l,
                                const ghl_eos_parameters *restrict eos,
                                const ghl_metric_quantities *restrict ADM_metric_face,
                                const double cmin,
                                const double cmax,
                                ghl_conservative_quantities *restrict cons_fluxes);

#define PPM_STENCIL_SIZE (5)
#define STENCIL_SIZE     (PPM_STENCIL_SIZE + 1)

#define GRHAYLMHD_COMPUTE_REMAINING_PRIMS_AT_FACE(m, p)                           \
    {                                                                             \
        p.temperature      = temperature[ijk];                                    \
        bool speed_limited = false;                                               \
        (void)ghl_limit_v_and_compute_u0(ghl_params, &m, &p, &speed_limited);     \
        ghl_tabulated_enforce_bounds_rho_Ye_P(ghl_eos, &p.rho, &p.Y_e, &p.press); \
        ghl_tabulated_compute_eps_T_from_P(ghl_eos,                               \
                                           p.rho,                                 \
                                           p.Y_e,                                 \
                                           p.press,                               \
                                           &p.eps,                                \
                                           &p.temperature);                       \
    }

#define GRHAYLMHD_RECONSTRUCT_HD(m, pr, pl, ps, dir)                          \
    {                                                                         \
        CCTK_REAL ftilde[2];                                                  \
        ghl_compute_ftilde(ghl_params, ps.press, ps.vU[dir], ftilde);         \
        ghl_ppm_reconstruction_with_steepening(ghl_params,                    \
                                               ps.press,                      \
                                               1.0,                           \
                                               ftilde,                        \
                                               ps.rho,                        \
                                               &pr.rho,                       \
                                               &pl.rho);                      \
        ghl_ppm_reconstruction(ftilde, ps.press, &pr.press, &pl.press);       \
        ghl_ppm_reconstruction(ftilde, ps.vU[0], &pr.vU[0], &pl.vU[0]);       \
        ghl_ppm_reconstruction(ftilde, ps.vU[1], &pr.vU[1], &pl.vU[1]);       \
        ghl_ppm_reconstruction(ftilde, ps.vU[2], &pr.vU[2], &pl.vU[2]);       \
        ghl_ppm_reconstruction(ftilde, ps.Y_e, &pr.Y_e, &pl.Y_e);             \
        ghl_ppm_reconstruction(ftilde, ps.entropy, &pr.entropy, &pl.entropy); \
        GRHAYLMHD_COMPUTE_REMAINING_PRIMS_AT_FACE(m, pr);                     \
        GRHAYLMHD_COMPUTE_REMAINING_PRIMS_AT_FACE(m, pl);                     \
    }

typedef struct {
    double rho[STENCIL_SIZE];
    double press[STENCIL_SIZE];
    double Y_e[STENCIL_SIZE];
    double vU[3][STENCIL_SIZE];
    double entropy[STENCIL_SIZE];
} ppm_reconstruction_stencil;

#define GRHAYLMHD_INTERP_TO_FACE(gf, m2, m1, p1) \
    ((-1.0 / 16.0) * (gf[p1] + gf[m2]) + (9.0 / 16.0) * (gf[ijk] + gf[m1]))
#define GRHAYLMHD_INTERP_METRIC_TO_FACE(face_struct, flux_dir)                     \
    {                                                                              \
        const int di = kronecker_delta[flux_dir][0];                               \
        const int dj = kronecker_delta[flux_dir][1];                               \
        const int dk = kronecker_delta[flux_dir][2];                               \
        const int m2 = CCTK_GFINDEX3D(cctkGH, i - 2 * di, j - 2 * dj, k - 2 * dk); \
        const int m1 = CCTK_GFINDEX3D(cctkGH, i - 1 * di, j - 1 * dj, k - 1 * dk); \
        const int p1 = CCTK_GFINDEX3D(cctkGH, i + 1 * di, j + 1 * dj, k + 1 * dk); \
                                                                                   \
        const CCTK_REAL face_alp   = GRHAYLMHD_INTERP_TO_FACE(alp, m2, m1, p1);    \
        const CCTK_REAL face_betax = GRHAYLMHD_INTERP_TO_FACE(betax, m2, m1, p1);  \
        const CCTK_REAL face_betay = GRHAYLMHD_INTERP_TO_FACE(betay, m2, m1, p1);  \
        const CCTK_REAL face_betaz = GRHAYLMHD_INTERP_TO_FACE(betaz, m2, m1, p1);  \
        const CCTK_REAL face_gxx   = GRHAYLMHD_INTERP_TO_FACE(gxx, m2, m1, p1);    \
        const CCTK_REAL face_gxy   = GRHAYLMHD_INTERP_TO_FACE(gxy, m2, m1, p1);    \
        const CCTK_REAL face_gxz   = GRHAYLMHD_INTERP_TO_FACE(gxz, m2, m1, p1);    \
        const CCTK_REAL face_gyy   = GRHAYLMHD_INTERP_TO_FACE(gyy, m2, m1, p1);    \
        const CCTK_REAL face_gyz   = GRHAYLMHD_INTERP_TO_FACE(gyz, m2, m1, p1);    \
        const CCTK_REAL face_gzz   = GRHAYLMHD_INTERP_TO_FACE(gzz, m2, m1, p1);    \
        ghl_initialize_metric(face_alp,                                            \
                              face_betax,                                          \
                              face_betay,                                          \
                              face_betaz,                                          \
                              face_gxx,                                            \
                              face_gxy,                                            \
                              face_gxz,                                            \
                              face_gyy,                                            \
                              face_gyz,                                            \
                              face_gzz,                                            \
                              &face_struct);                                       \
    }

void GRHayLMHD_Compute_RHSs_Fluxes(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    const int imin = cctk_nghostzones[0], imax = cctk_lsh[0] - cctk_nghostzones[0];
    const int jmin = cctk_nghostzones[1], jmax = cctk_lsh[1] - cctk_nghostzones[1];
    const int kmin = cctk_nghostzones[2], kmax = cctk_lsh[2] - cctk_nghostzones[2];

    const CCTK_REAL inv_h[3] = {
        1.0 / CCTK_DELTA_SPACE(0),
        1.0 / CCTK_DELTA_SPACE(1),
        1.0 / CCTK_DELTA_SPACE(2),
    };

    static const CCTK_REAL kronecker_delta[3][3] = {
        { 1, 0, 0 },
        { 0, 1, 0 },
        { 0, 0, 1 },
    };

    ghl_cmin_cmax_func_t ghl_cmin_cmax_funcs[3] = {
        ghl_calculate_characteristic_speed_dirn0,
        ghl_calculate_characteristic_speed_dirn1,
        ghl_calculate_characteristic_speed_dirn2,
    };

    ghl_flux_func_t ghl_flux_funcs[3] = {
        ghl_calculate_HLLE_fluxes_dirn0_tabulated_entropy,
        ghl_calculate_HLLE_fluxes_dirn1_tabulated_entropy,
        ghl_calculate_HLLE_fluxes_dirn2_tabulated_entropy,
    };

    for(int flux_dir = 0; flux_dir < 3; flux_dir++) {
        OMPLOOP3D(imin, imax + 1, jmin, jmax + 1, kmin, kmax + 1)
        {
            ghl_metric_quantities metric_face = { 0 };
            GRHAYLMHD_INTERP_METRIC_TO_FACE(metric_face, flux_dir)

            // Stencil is points -3, -2, -1, 0, +1, +2
            ppm_reconstruction_stencil prims_stencil = { 0 };
            for(int m = 0; m < STENCIL_SIZE; m++) {
                const int ii  = i + kronecker_delta[flux_dir][0] * (m - 3);
                const int jj  = j + kronecker_delta[flux_dir][1] * (m - 3);
                const int kk  = k + kronecker_delta[flux_dir][2] * (m - 3);
                const int idx = CCTK_GFINDEX3D(cctkGH, ii, jj, kk);

                prims_stencil.rho[m]     = rho[idx];
                prims_stencil.press[m]   = press[idx];
                prims_stencil.vU[0][m]   = vx[idx];
                prims_stencil.vU[1][m]   = vy[idx];
                prims_stencil.vU[2][m]   = vz[idx];
                prims_stencil.Y_e[m]     = Y_e[idx];
                prims_stencil.entropy[m] = entropy[idx];
            }

            ghl_primitive_quantities prims_r = { 0 }, prims_l = { 0 };
            GRHAYLMHD_RECONSTRUCT_HD(metric_face, prims_r, prims_l, prims_stencil, flux_dir);

            double cmin = 0, cmax = 0;
            ghl_cmin_cmax_funcs[flux_dir](&prims_r, &prims_l, ghl_eos, &metric_face, &cmin, &cmax);

            ghl_conservative_quantities fluxes = { 0 };
            ghl_flux_funcs[flux_dir](&prims_r,
                                     &prims_l,
                                     ghl_eos,
                                     &metric_face,
                                     cmin,
                                     cmax,
                                     &fluxes);

            rho_flux[ijk] = fluxes.rho;
            tau_flux[ijk] = fluxes.tau;
            S_x_flux[ijk] = fluxes.SD[0];
            S_y_flux[ijk] = fluxes.SD[1];
            S_z_flux[ijk] = fluxes.SD[2];
            Y_e_flux[ijk] = fluxes.Y_e;
            ent_flux[ijk] = fluxes.entropy;
        }
        ENDLOOP3D

        OMPLOOP3D(imin, imax, jmin, jmax, kmin, kmax)
        {
            const int ip     = i + kronecker_delta[flux_dir][0];
            const int jp     = j + kronecker_delta[flux_dir][1];
            const int kp     = k + kronecker_delta[flux_dir][2];
            const int ijk_p1 = CCTK_GFINDEX3D(cctkGH, ip, jp, kp);

            rho_rhs[ijk] -= (rho_flux[ijk_p1] - rho_flux[ijk]) * inv_h[flux_dir];
            tau_rhs[ijk] -= (tau_flux[ijk_p1] - tau_flux[ijk]) * inv_h[flux_dir];
            S_x_rhs[ijk] -= (S_x_flux[ijk_p1] - S_x_flux[ijk]) * inv_h[flux_dir];
            S_y_rhs[ijk] -= (S_y_flux[ijk_p1] - S_y_flux[ijk]) * inv_h[flux_dir];
            S_z_rhs[ijk] -= (S_z_flux[ijk_p1] - S_z_flux[ijk]) * inv_h[flux_dir];
            Y_e_rhs[ijk] -= (Y_e_flux[ijk_p1] - Y_e_flux[ijk]) * inv_h[flux_dir];
            ent_rhs[ijk] -= (ent_flux[ijk_p1] - ent_flux[ijk]) * inv_h[flux_dir];
        }
        ENDLOOP3D
    }
}
