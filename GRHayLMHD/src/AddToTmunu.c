#include "GRHayLMHD.h"

#define DEBUG_PRINT_METRIC(struct_name)                             \
    CCTK_VINFO("lapse         = %.15e", struct_name.lapse);         \
    CCTK_VINFO("betaU[0]      = %.15e", struct_name.betaU[0]);      \
    CCTK_VINFO("betaU[1]      = %.15e", struct_name.betaU[1]);      \
    CCTK_VINFO("betaU[2]      = %.15e", struct_name.betaU[2]);      \
    CCTK_VINFO("gammaDD[0][0] = %.15e", struct_name.gammaDD[0][0]); \
    CCTK_VINFO("gammaDD[0][1] = %.15e", struct_name.gammaDD[0][1]); \
    CCTK_VINFO("gammaDD[0][2] = %.15e", struct_name.gammaDD[0][2]); \
    CCTK_VINFO("gammaDD[1][1] = %.15e", struct_name.gammaDD[1][1]); \
    CCTK_VINFO("gammaDD[1][2] = %.15e", struct_name.gammaDD[1][2]); \
    CCTK_VINFO("gammaDD[2][2] = %.15e", struct_name.gammaDD[2][2]); \
    CCTK_VINFO("gammaUU[0][0] = %.15e", struct_name.gammaUU[0][0]); \
    CCTK_VINFO("gammaUU[0][1] = %.15e", struct_name.gammaUU[0][1]); \
    CCTK_VINFO("gammaUU[0][2] = %.15e", struct_name.gammaUU[0][2]); \
    CCTK_VINFO("gammaUU[1][1] = %.15e", struct_name.gammaUU[1][1]); \
    CCTK_VINFO("gammaUU[1][2] = %.15e", struct_name.gammaUU[1][2]); \
    CCTK_VINFO("gammaUU[2][2] = %.15e", struct_name.gammaUU[2][2]);

#define DEBUG_PRINT_PRIMS(struct_name)                          \
    CCTK_VINFO("rho         = %.15e", struct_name.rho);         \
    CCTK_VINFO("press       = %.15e", struct_name.press);       \
    CCTK_VINFO("eps         = %.15e", struct_name.eps);         \
    CCTK_VINFO("u0          = %.15e", struct_name.u0);          \
    CCTK_VINFO("vU[0]       = %.15e", struct_name.vU[0]);       \
    CCTK_VINFO("vU[1]       = %.15e", struct_name.vU[1]);       \
    CCTK_VINFO("vU[2]       = %.15e", struct_name.vU[2]);       \
    CCTK_VINFO("Y_e         = %.15e", struct_name.Y_e);         \
    CCTK_VINFO("temperature = %.15e", struct_name.temperature); \
    CCTK_VINFO("entropy     = %.15e", struct_name.entropy);

void GRHayLMHD_AddToTmunu(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    const int imin = 0, imax = cctk_lsh[0];
    const int jmin = 0, jmax = cctk_lsh[1];
    const int kmin = 0, kmax = cctk_lsh[2];
    bool      stop = false;

    OMPLOOP3D(imin, imax, jmin, jmax, kmin, kmax)
    {
        ghl_metric_quantities adm_metric = { 0 };
        GRHAYLMHD_PACK_METRIC(adm_metric);

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
