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

#define GRHAYLMHD_COMPUTE_REMAINING_PRIMS_AT_FACE(struct_name) \
    struct_name.temperature = temperature[ijk];                \
    ghl_tabulated_compute_eps_T_from_P(ghl_eos,                \
                                       struct_name.rho,        \
                                       struct_name.Y_e,        \
                                       struct_name.press,      \
                                       &struct_name.eps,       \
                                       &struct_name.temperature);

#define GRHAYLMHD_RECONSTRUCT_HD(pr, pl, ps, dir)                       \
    {                                                                   \
        CCTK_REAL ftilde[2];                                            \
        ghl_compute_ftilde(ghl_params, ps.press, ps.vU[dir], ftilde);   \
        ghl_ppm_reconstruction_with_steepening(ghl_params,              \
                                               ps.press,                \
                                               1.0,                     \
                                               ftilde,                  \
                                               ps.rho,                  \
                                               &pr.rho,                 \
                                               &pl.rho);                \
        ghl_ppm_reconstruction(ftilde, ps.press, &pr.press, &pl.press); \
        ghl_ppm_reconstruction(ftilde, ps.Y_e, &pr.Y_e, &pl.press);     \
        ghl_ppm_reconstruction(ftilde, ps.vU[0], &pr.vU[0], &pl.vU[0]); \
        ghl_ppm_reconstruction(ftilde, ps.vU[1], &pr.vU[1], &pl.vU[1]); \
        ghl_ppm_reconstruction(ftilde, ps.vU[2], &pr.vU[2], &pl.vU[2]); \
        GRHAYLMHD_COMPUTE_REMAINING_PRIMS_AT_FACE(pr);                  \
        GRHAYLMHD_COMPUTE_REMAINING_PRIMS_AT_FACE(pl);                  \
    }

typedef struct {
    double rho[STENCIL_SIZE];
    double press[STENCIL_SIZE];
    double Y_e[STENCIL_SIZE];
    double vU[3][STENCIL_SIZE];
    double entropy[STENCIL_SIZE];
} ppm_reconstruction_stencil;

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
            ghl_metric_quantities adm_metric = { 0 };
            GRHAYLMHD_PACK_METRIC_ENFORCE_DETGAMMAEQ1(adm_metric);

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
            GRHAYLMHD_RECONSTRUCT_HD(prims_r, prims_l, prims_stencil, flux_dir);

            double cmin = 0, cmax = 0;
            ghl_cmin_cmax_funcs[flux_dir](&prims_r, &prims_l, ghl_eos, &adm_metric, &cmin, &cmax);

            ghl_conservative_quantities fluxes = { 0 };
            ghl_flux_funcs[flux_dir](&prims_r,
                                     &prims_l,
                                     ghl_eos,
                                     &adm_metric,
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
