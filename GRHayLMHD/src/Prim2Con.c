#include "GRHayLMHD.h"

void GRHayLMHD_Prim2Con_SinglePoint(CCTK_ARGUMENTS,
                                    const int ijk,
                                    const int ijkx,
                                    const int ijky,
                                    const int ijkz)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    ghl_metric_quantities adm_metric = { 0 };
    GRHAYLMHD_LOAD_METRIC(adm_metric);

    ghl_ADM_aux_quantities aux_metric = { 0 };
    ghl_compute_ADM_auxiliaries(&adm_metric, &aux_metric);

    ghl_primitive_quantities prims = { 0 };
    GRHAYLMHD_LOAD_PRIMS(prims);

    bool speed_limited;
    (void)ghl_enforce_primitive_limits_and_compute_u0(ghl_params,
                                                      ghl_eos,
                                                      &adm_metric,
                                                      &prims,
                                                      &speed_limited);

    ghl_conservative_quantities cons = { 0 };
    ghl_compute_conservs(&adm_metric, &aux_metric, &prims, &cons);

    GRHAYLMHD_WRITE_PRIMS(prims);
    GRHAYLMHD_WRITE_CONS(cons);
}

void GRHayLMHD_Prim2Con(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    const int imin = 0, imax = cctk_lsh[0];
    const int jmin = 0, jmax = cctk_lsh[1];
    const int kmin = 0, kmax = cctk_lsh[2];

    OMPLOOP3D(imin, imax, jmin, jmax, kmin, kmax)
    {
        GRHayLMHD_Prim2Con_SinglePoint(CCTK_PASS_CTOC, ijk, ijkx, ijky, ijkz);
    }
    ENDLOOP3D
}
