#include "GRHayLMHD.h"

#include <assert.h>

#define GRHAYLMHD_COMPUTE_ERROR_CON(con)                 \
    errors_numer_##con = fabs(cons.con - cons_orig.con); \
    errors_denom_##con = fabs(cons_orig.con);

#define GRHAYLMHD_COMPUTE_ERRORS       \
    GRHAYLMHD_COMPUTE_ERROR_CON(rho)   \
    GRHAYLMHD_COMPUTE_ERROR_CON(tau)   \
    GRHAYLMHD_COMPUTE_ERROR_CON(SD[0]) \
    GRHAYLMHD_COMPUTE_ERROR_CON(SD[1]) \
    GRHAYLMHD_COMPUTE_ERROR_CON(SD[2]) \
    GRHAYLMHD_COMPUTE_ERROR_CON(Y_e)   \
    GRHAYLMHD_COMPUTE_ERROR_CON(entropy)

#define GRHAYLMHD_COMPUTE_FINAL_ERROR_OF(s, v) \
    s.v = errors_denom_##v == 0 ? 0 : errors_numer_##v / errors_denom_##v;

#define GRHAYLMHD_PRINT_CON2PRIM_DIAGNOSTICS                                                 \
    {                                                                                        \
        ghl_conservative_quantities cons_errors = { 0 };                                     \
        GRHAYLMHD_COMPUTE_FINAL_ERROR_OF(cons_errors, rho);                                  \
        GRHAYLMHD_COMPUTE_FINAL_ERROR_OF(cons_errors, tau);                                  \
        GRHAYLMHD_COMPUTE_FINAL_ERROR_OF(cons_errors, SD[0]);                                \
        GRHAYLMHD_COMPUTE_FINAL_ERROR_OF(cons_errors, SD[1]);                                \
        GRHAYLMHD_COMPUTE_FINAL_ERROR_OF(cons_errors, SD[2]);                                \
        GRHAYLMHD_COMPUTE_FINAL_ERROR_OF(cons_errors, Y_e);                                  \
        GRHAYLMHD_COMPUTE_FINAL_ERROR_OF(cons_errors, entropy);                              \
        const int ntotal = (imax - imin) * (jmax - jmin) * (kmax - kmin);                    \
        CCTK_VINFO("Con2Prim -- It. %d -- Ref. Lev. %d -- ATM Resets: %d / %d",              \
                   cctk_iteration,                                                           \
                   GetRefinementLevel(cctkGH),                                               \
                   atm_resets,                                                               \
                   ntotal);                                                                  \
        CCTK_VINFO("  rho %.2e -- tau %.2e -- S_x %.2e -- S_y %.2e -- S_z %.2e -- Y_e %.2e " \
                   "-- ent %.2e",                                                            \
                   cons_errors.rho,                                                          \
                   cons_errors.tau,                                                          \
                   cons_errors.SD[0],                                                        \
                   cons_errors.SD[1],                                                        \
                   cons_errors.SD[2],                                                        \
                   cons_errors.Y_e,                                                          \
                   cons_errors.entropy);                                                     \
    }

void GRHayLMHD_Con2Prim(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // TODO: move this to initialization
    assert(ghl_params->calc_prim_guess == true);

    const int imin = cctk_nghostzones[0];
    const int jmin = cctk_nghostzones[1];
    const int kmin = cctk_nghostzones[2];

    const int imax = cctk_lsh[0] - cctk_nghostzones[0];
    const int jmax = cctk_lsh[1] - cctk_nghostzones[1];
    const int kmax = cctk_lsh[2] - cctk_nghostzones[2];

    int atm_resets = 0;

    CCTK_REAL errors_numer_rho = 0.0, errors_numer_tau = 0.0, errors_numer_SD[3] = { 0 },
              errors_numer_Y_e = 0.0, errors_numer_entropy = 0.0;
    CCTK_REAL errors_denom_rho = 0.0, errors_denom_tau = 0.0, errors_denom_SD[3] = { 0 },
              errors_denom_Y_e = 0.0, errors_denom_entropy = 0.0;

    // clang-format off
#define ERROR_ARGS                                                              \
    errors_numer_rho, errors_numer_tau, errors_numer_Y_e, errors_numer_entropy, \
    errors_denom_rho, errors_denom_tau, errors_denom_Y_e, errors_denom_entropy, \
    errors_numer_SD[:3], errors_denom_SD[:3]
    // clang-format on

#pragma omp parallel for reduction(+ : atm_resets, ERROR_ARGS)
    LOOP3D(imin, imax, jmin, jmax, kmin, kmax)
    {
        ghl_con2prim_diagnostics diagnostics = { 0 };

        ghl_metric_quantities adm_metric = { 0 };
        GRHAYLMHD_LOAD_METRIC_ENFORCE_DETGAMMAEQ1(adm_metric);

        ghl_ADM_aux_quantities aux_metric = { 0 };
        ghl_compute_ADM_auxiliaries(&adm_metric, &aux_metric);

        ghl_conservative_quantities cons_orig = { 0 };
        GRHAYLMHD_LOAD_CONS(cons_orig);
        ghl_conservative_quantities cons = cons_orig;

        ghl_conservative_quantities cons_undens = { 0 };
        ghl_undensitize_conservatives(adm_metric.sqrt_detgamma, &cons, &cons_undens);

        ghl_primitive_quantities prims = { 0 };

        ghl_error_codes_t error = ghl_con2prim_tabulated_multi_method(ghl_params,
                                                                      ghl_eos,
                                                                      &adm_metric,
                                                                      &aux_metric,
                                                                      &cons_undens,
                                                                      &prims,
                                                                      &diagnostics);

        // TODO: averaging algorithm
        if(error) {
            atm_resets++;
            ghl_set_prims_to_constant_atm(ghl_eos, &prims);
        }

        bool speed_limited = false;
        ghl_enforce_primitive_limits_and_compute_u0(ghl_params,
                                                    ghl_eos,
                                                    &adm_metric,
                                                    &prims,
                                                    &speed_limited);

        ghl_compute_conservs(&adm_metric, &aux_metric, &prims, &cons);

        GRHAYLMHD_WRITE_PRIMS(prims);
        GRHAYLMHD_WRITE_CONS(cons);
        GRHAYLMHD_COMPUTE_ERRORS;
    }
    ENDLOOP3D

    GRHAYLMHD_PRINT_CON2PRIM_DIAGNOSTICS;
}
