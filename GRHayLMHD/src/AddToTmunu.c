#include "GRHayLMHD.h"

void GRHayLMHD_AddToTmunu(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    const int imin = 0, imax = cctk_lsh[0];
    const int jmin = 0, jmax = cctk_lsh[1];
    const int kmin = 0, kmax = cctk_lsh[2];

    OMPLOOP3D(imin, imax, jmin, jmax, kmin, kmax)
    {
        ghl_metric_quantities adm_metric = { 0 };
        GRHAYLMHD_PACK_METRIC_ENFORCE_DETGAMMAEQ1(adm_metric);

        ghl_ADM_aux_quantities aux_metric = { 0 };
        ghl_compute_ADM_auxiliaries(&adm_metric, &aux_metric);

        ghl_primitive_quantities prims = { 0 };
        GRHAYLMHD_PACK_PRIMS(prims);

        ghl_stress_energy Tmunu = { 0 };
        ghl_compute_TDNmunu(&adm_metric, &aux_metric, &prims, &Tmunu);

        eTtt[ijk] += Tmunu.T4[0][0];
        eTtx[ijk] += Tmunu.T4[0][1];
        eTty[ijk] += Tmunu.T4[0][2];
        eTtz[ijk] += Tmunu.T4[0][3];
        eTxx[ijk] += Tmunu.T4[1][1];
        eTxy[ijk] += Tmunu.T4[1][2];
        eTxz[ijk] += Tmunu.T4[1][3];
        eTyy[ijk] += Tmunu.T4[2][2];
        eTyz[ijk] += Tmunu.T4[2][3];
        eTzz[ijk] += Tmunu.T4[3][3];
    }
    ENDLOOP3D
}
