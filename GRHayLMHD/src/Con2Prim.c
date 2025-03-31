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

#define GRHAYLMHD_PRINT_CON2PRIM_HEADER

static void print_header(const int it)
{
    DECLARE_CCTK_PARAMETERS;
    static int prev_iteration = 0;

    const bool new_iteration = (prev_iteration != it);
    const bool time_to_print = (it == 1 || (it % output_info_header_every) == 0);
    if(new_iteration && time_to_print) {
        // clang-format off
        prev_iteration = it;
        puts("--------.--------.----.-------..---------.---------.---------.---------.---------.---------.---------");
        puts("     It |  Time  | RL | ATM % ||   rho   |   tau   |   S_x   |   S_y   |   S_z   |   Y_e   |  Ent");
        puts("--------.--------.----.-------..---------.---------.---------.---------.---------.---------.---------");
        // clang-format on
    }
}

#define GRHAYLMHD_PRINT_CON2PRIM_DIAGNOSTICS                                            \
    {                                                                                   \
        const int ref_lev = GetRefinementLevel(cctkGH);                                 \
        print_header(cctk_iteration);                                                   \
        ghl_conservative_quantities cons_errors = { 0 };                                \
        GRHAYLMHD_COMPUTE_FINAL_ERROR_OF(cons_errors, rho);                             \
        GRHAYLMHD_COMPUTE_FINAL_ERROR_OF(cons_errors, tau);                             \
        GRHAYLMHD_COMPUTE_FINAL_ERROR_OF(cons_errors, SD[0]);                           \
        GRHAYLMHD_COMPUTE_FINAL_ERROR_OF(cons_errors, SD[1]);                           \
        GRHAYLMHD_COMPUTE_FINAL_ERROR_OF(cons_errors, SD[2]);                           \
        GRHAYLMHD_COMPUTE_FINAL_ERROR_OF(cons_errors, Y_e);                             \
        GRHAYLMHD_COMPUTE_FINAL_ERROR_OF(cons_errors, entropy);                         \
        const CCTK_INT ntotal = (imax - imin) * (jmax - jmin) * (kmax - kmin);          \
        printf("%8d %.2e %3d %7.1f  %9.1e   %.1e   %.1e   %.1e   %.1e   %.1e   %.1e\n", \
               cctk_iteration,                                                          \
               cctk_time,                                                               \
               ref_lev,                                                                 \
               100.0 * atm_resets / (CCTK_REAL)ntotal,                                  \
               cons_errors.rho,                                                         \
               cons_errors.tau,                                                         \
               cons_errors.SD[0],                                                       \
               cons_errors.SD[1],                                                       \
               cons_errors.SD[2],                                                       \
               cons_errors.Y_e,                                                         \
               cons_errors.entropy);                                                    \
    }

static ghl_error_codes_t
igm_entropy_con2prim_multi_method(const ghl_parameters *restrict params,
                                  const ghl_eos_parameters *restrict eos,
                                  const ghl_metric_quantities *restrict ADM_metric,
                                  const ghl_ADM_aux_quantities *restrict metric_aux,
                                  const ghl_conservative_quantities *restrict cons,
                                  ghl_primitive_quantities *restrict prims,
                                  ghl_con2prim_diagnostics *restrict diagnostics)
{

    if(params->calc_prim_guess) {
        ghl_guess_primitives(eos, ADM_metric, cons, prims);
    }

    // Store primitive guesses (used if con2prim fails)
    const ghl_primitive_quantities prims_guess = *prims;

    // Backup 1 triggered
    diagnostics->backup[1]  = true;
    ghl_error_codes_t error = ghl_tabulated_Palenzuela1D_entropy(params,
                                                                 eos,
                                                                 ADM_metric,
                                                                 metric_aux,
                                                                 cons,
                                                                 prims,
                                                                 diagnostics);

    if(error) {
        // Reset guesses
        *prims = prims_guess;

        // Backup routine #1
        error = ghl_tabulated_Palenzuela1D_energy(params,
                                                  eos,
                                                  ADM_metric,
                                                  metric_aux,
                                                  cons,
                                                  prims,
                                                  diagnostics);

        if(error) {
            diagnostics->backup[2] = false;

            // Reset guesses
            *prims = prims_guess;

            // Backup routine #2
            error = ghl_tabulated_Newman1D_entropy(params,
                                                   eos,
                                                   ADM_metric,
                                                   metric_aux,
                                                   cons,
                                                   prims,
                                                   diagnostics);

            if(error) {
                // Backup 3 triggered
                diagnostics->backup[0] = true;

                // Reset guesses
                *prims = prims_guess;

                // Backup routine #3
                error = ghl_tabulated_Newman1D_energy(params,
                                                      eos,
                                                      ADM_metric,
                                                      metric_aux,
                                                      cons,
                                                      prims,
                                                      diagnostics);
            }
        }
    }
    return error;
}

static inline bool
depsdT_check(const ghl_eos_parameters *eos, const double rho, const double Y_e, const double T)
{
    DECLARE_CCTK_PARAMETERS;

    double P      = 0.0;
    double eps    = 0.0;
    double depsdT = 0.0;
    ghl_tabulated_compute_P_eps_depsdT_from_T(eos, rho, Y_e, T, &P, &eps, &depsdT);

    return depsdT < depsdT_threshold;
}

static inline double compute_W(const double q,
                               const double r,
                               const double s,
                               const double t,
                               const double x,
                               const double inv_W_max_squared)
{
    double Wminus2 = 1.0 - (x * x * r + (2 * x + s) * t * t) / (x * x * (x + s) * (x + s));
    Wminus2        = fmin(fmax(Wminus2, inv_W_max_squared), 1.0);
    return pow(Wminus2, -0.5);
}

static inline bool use_entropy_check(const ghl_eos_parameters          *eos,
                                     const ghl_metric_quantities       *ADM_metric,
                                     const ghl_conservative_quantities *cons_undens,
                                     const double                       T_guess,
                                     const double                       inv_W_max_squared)
{
    double       SD[3]     = { cons_undens->SD[0], cons_undens->SD[1], cons_undens->SD[2] };
    double       S_squared = ghl_compute_vec2_from_vec3D(ADM_metric->gammaUU, SD);
    const double tau       = MAX(cons_undens->tau, 0.99 * eos->tau_atm);

    const double S_squared_max = SQR(tau + cons_undens->rho);
    if(S_squared > S_squared_max) {
        // Step 2.2: Rescale S_{i}
        const double rescale_factor = sqrt(0.9999 * S_squared_max / S_squared);
        for(int i = 0; i < 3; i++) {
            SD[i] *= rescale_factor;
        }

        // Step 2.3: Recompute S^{2}
        S_squared = ghl_compute_vec2_from_vec3D(ADM_metric->gammaUU, SD);
    }

    const double invD = 1.0 / cons_undens->rho;
    const double Y_e  = cons_undens->Y_e * invD;
    const double q    = tau * invD;
    const double r    = S_squared * invD * invD;
    // zero magnetic fields
    const double s    = 0.0;
    const double t    = 0.0;

    const double xlow     = 1.0 + q - s;
    const double xlow_W   = compute_W(q, r, s, t, xlow, inv_W_max_squared);
    const double xlow_rho = cons_undens->rho / xlow_W;
    if(depsdT_check(eos, xlow_rho, Y_e, T_guess)) {
        return true;
    }

    const double xup     = 2.0 + 2.0 * q - s;
    const double xup_W   = compute_W(q, r, s, t, xup, inv_W_max_squared);
    const double xup_rho = cons_undens->rho / xup_W;
    if(depsdT_check(eos, xup_rho, Y_e, T_guess)) {
        return true;
    }

    return false;
}

void GRHayLMHD_Con2Prim(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // TODO: move this to initialization
    assert(ghl_params->calc_prim_guess == true);

    // WARNING: these bounds *do not* match IGM's, but I think they are correct.
    // const int imin = cctk_nghostzones[0], imax = cctk_lsh[0] -
    // cctk_nghostzones[0]; const int jmin = cctk_nghostzones[1], jmax =
    // cctk_lsh[1] - cctk_nghostzones[1]; const int kmin = cctk_nghostzones[2],
    // kmax = cctk_lsh[2] - cctk_nghostzones[2];
    const int imin = 0, imax = cctk_lsh[0];
    const int jmin = 0, jmax = cctk_lsh[1];
    const int kmin = 0, kmax = cctk_lsh[2];

    CCTK_INT  atm_resets           = 0;
    CCTK_REAL errors_numer_rho     = 0.0;
    CCTK_REAL errors_numer_tau     = 0.0;
    CCTK_REAL errors_numer_Y_e     = 0.0;
    CCTK_REAL errors_numer_entropy = 0.0;
    CCTK_REAL errors_numer_SD[3]   = { 0.0, 0.0, 0.0 };
    CCTK_REAL errors_denom_rho     = 0.0;
    CCTK_REAL errors_denom_tau     = 0.0;
    CCTK_REAL errors_denom_Y_e     = 0.0;
    CCTK_REAL errors_denom_entropy = 0.0;
    CCTK_REAL errors_denom_SD[3]   = { 0.0, 0.0, 0.0 };

    // clang-format off
#define REDUCTION_VARIABLES                                                     \
    atm_resets, errors_numer_SD[:3], errors_denom_SD[:3],                       \
    errors_numer_rho, errors_numer_tau, errors_numer_Y_e, errors_numer_entropy, \
    errors_denom_rho, errors_denom_tau, errors_denom_Y_e, errors_denom_entropy
    // clang-format on

#pragma omp parallel for reduction(+ : REDUCTION_VARIABLES) schedule(static)
    LOOP3D(imin, imax, jmin, jmax, kmin, kmax)
    {
        ghl_con2prim_diagnostics diagnostics = { 0 };

        ghl_metric_quantities adm_metric = { 0 };
        GRHAYLMHD_PACK_METRIC(adm_metric);

        ghl_ADM_aux_quantities aux_metric = { 0 };
        ghl_compute_ADM_auxiliaries(&adm_metric, &aux_metric);

        ghl_conservative_quantities cons_orig = { 0 };
        GRHAYLMHD_PACK_CONS(cons_orig);
        ghl_conservative_quantities cons = cons_orig;

        ghl_conservative_quantities cons_undens = { 0 };
        ghl_undensitize_conservatives(adm_metric.sqrt_detgamma, &cons, &cons_undens);

        ghl_primitive_quantities prims = { 0 };

        // Check whether or not to use the entropy equation first, following
        // the logic inside IGM's con2prim. Note we hardcode T_guess for now.
        const double T_guess = ghl_eos->T_max;

        const bool use_entropy = use_entropy_check(ghl_eos,
                                                   &adm_metric,
                                                   &cons_undens,
                                                   T_guess,
                                                   ghl_params->inv_sq_max_Lorentz_factor);

        ghl_error_codes_t error = ghl_success;
        if(use_entropy) {
            error = igm_entropy_con2prim_multi_method(ghl_params,
                                                      ghl_eos,
                                                      &adm_metric,
                                                      &aux_metric,
                                                      &cons_undens,
                                                      &prims,
                                                      &diagnostics);
        }
        else {
            error = ghl_con2prim_tabulated_multi_method(ghl_params,
                                                        ghl_eos,
                                                        &adm_metric,
                                                        &aux_metric,
                                                        &cons_undens,
                                                        &prims,
                                                        &diagnostics);
        }

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

        GRHAYLMHD_UNPACK_PRIMS(prims);
        GRHAYLMHD_UNPACK_CONS(cons);
        GRHAYLMHD_COMPUTE_ERRORS;
    }
    ENDLOOP3D

    static int rk_substep = -1;
    if(cctk_iteration > 0) {
        rk_substep = (rk_substep + 1) % 4;
    }
    GRHAYLMHD_PRINT_CON2PRIM_DIAGNOSTICS;
}
