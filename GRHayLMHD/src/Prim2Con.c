#include "GRHayLMHD.h"

void GRHayLMHD_Prim2Con_SinglePoint(CCTK_ARGUMENTS,
                                    const int ijk,
                                    const int ijkx,
                                    const int ijky,
                                    const int ijkz)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    ghl_metric_quantities ADM_metric = { 0 };
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
                          &ADM_metric);

    ghl_ADM_aux_quantities AUX_metric = { 0 };
    ghl_compute_ADM_auxiliaries(&ADM_metric, &AUX_metric);

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

    bool speed_limited;
    (void)ghl_enforce_primitive_limits_and_compute_u0(ghl_params,
                                                      ghl_eos,
                                                      &ADM_metric,
                                                      &prims,
                                                      &speed_limited);

    ghl_conservative_quantities cons = { 0 };
    ghl_compute_conservs(&ADM_metric, &AUX_metric, &prims, &cons);

    rho[ijk]         = prims.rho;
    press[ijk]       = prims.press;
    eps[ijk]         = prims.eps;
    entropy[ijk]     = prims.entropy;
    Y_e[ijk]         = prims.Y_e;
    temperature[ijk] = prims.temperature;
    u0[ijk]          = prims.u0;
    vx[ijk]          = prims.vU[0];
    vy[ijk]          = prims.vU[1];
    vz[ijk]          = prims.vU[2];

    rho_tilde[ijk] = cons.rho;
    tau_tilde[ijk] = cons.tau;
    S_x_tilde[ijk] = cons.SD[0];
    S_y_tilde[ijk] = cons.SD[1];
    S_z_tilde[ijk] = cons.SD[2];
    Y_e_tilde[ijk] = cons.Y_e;
    ent_tilde[ijk] = cons.entropy;
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
