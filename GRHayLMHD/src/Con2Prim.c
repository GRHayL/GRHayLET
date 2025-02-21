#include "GRHayLMHD.h"

#include <assert.h>

// TODO: add more diagnostic information
static void print_diagnostics(const int imin,
                              const int imax,
                              const int jmin,
                              const int jmax,
                              const int kmin,
                              const int kmax,
                              const int iteration,
                              const int reflevel,
                              const int atm_resets)
{
    const int ntotal = (imax - imin) * (jmax - jmin) * (kmax - kmin);
    CCTK_VINFO("Con2Prim -- It. %d -- Ref. Lev. %d -- ATM Resets: %d / %d",
               iteration,
               reflevel,
               atm_resets,
               ntotal);
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

#pragma omp parallel for reduction(+ : atm_resets) schedule(static)
    LOOP3D(imin, imax, jmin, jmax, kmin, kmax)
    {
        ghl_con2prim_diagnostics diagnostics = { 0 };

        ghl_metric_quantities adm_metric = { 0 };
        GRHAYLMHD_LOAD_METRIC_ENFORCE_DETGAMMAEQ1(adm_metric);

        ghl_ADM_aux_quantities aux_metric = { 0 };
        ghl_compute_ADM_auxiliaries(&adm_metric, &aux_metric);

        ghl_conservative_quantities cons = { 0 };
        GRHAYLMHD_LOAD_CONS(cons);

        ghl_conservative_quantities cons_undens = { 0 };
        ghl_undensitize_conservatives(adm_metric.sqrt_detgamma, &cons, &cons_undens);

        ghl_primitive_quantities prims = { 0 };
        prims.BU[0]                    = Bvec[ijkx];
        prims.BU[1]                    = Bvec[ijky];
        prims.BU[2]                    = Bvec[ijkz];

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
    }
    ENDLOOP3D

    print_diagnostics(imin,
                      imax,
                      jmin,
                      jmax,
                      kmin,
                      kmax,
                      cctk_iteration,
                      GetRefinementLevel(cctkGH),
                      atm_resets);
}
