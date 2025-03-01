#ifndef GRHAYLET_GRHAYLMHD_H
#define GRHAYLET_GRHAYLMHD_H

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <math.h>
#include <stdbool.h>

#include "GRHayLib.h"

enum GRHayLMHD_bc_dir {
    GRHayLMHD_xmin,
    GRHayLMHD_xmax,
    GRHayLMHD_ymin,
    GRHayLMHD_ymax,
    GRHayLMHD_zmin,
    GRHayLMHD_zmax
};

void GRHayLMHD_Prim2Con_SinglePoint(CCTK_ARGUMENTS, const int ijk);

#define LOOP3D(_imin, _imax, _jmin, _jmax, _kmin, _kmax)                                   \
    for(int k = _kmin; k < _kmax; k++) {                                                   \
        for(int j = _jmin; j < _jmax; j++) {                                               \
            for(int i = _imin; i < _imax; i++) {                                           \
                CCTK_ATTRIBUTE_UNUSED const int ijk  = CCTK_GFINDEX3D(cctkGH, i, j, k);    \
                CCTK_ATTRIBUTE_UNUSED const int ijkx = CCTK_GFINDEX4D(cctkGH, i, j, k, 0); \
                CCTK_ATTRIBUTE_UNUSED const int ijky = CCTK_GFINDEX4D(cctkGH, i, j, k, 1); \
                CCTK_ATTRIBUTE_UNUSED const int ijkz = CCTK_GFINDEX4D(cctkGH, i, j, k, 2);

#define OMPLOOP3D(_imin, _imax, _jmin, _jmax, _kmin, _kmax) \
    _Pragma("omp parallel for") LOOP3D(_imin, _imax, _jmin, _jmax, _kmin, _kmax)

#define ENDLOOP3D \
    }             \
    }             \
    }

// Useful macros for packing/unpacking GRHayL's structs from/to gridfunctions.
#define GRHAYLMHD_PACK_PRIMS(struct_name)       \
    struct_name.u0 = u0[ijk];                   \
    ghl_initialize_primitives(rho[ijk],         \
                              press[ijk],       \
                              eps[ijk],         \
                              vx[ijk],          \
                              vy[ijk],          \
                              vz[ijk],          \
                              0.0,              \
                              0.0,              \
                              0.0,              \
                              entropy[ijk],     \
                              Y_e[ijk],         \
                              temperature[ijk], \
                              &struct_name);

#define GRHAYLMHD_UNPACK_PRIMS(struct_name)       \
    u0[ijk] = struct_name.u0;                     \
    {                                             \
        double BU[3] = { 0.0, 0.0, 0.0 };         \
        ghl_return_primitives(&struct_name,       \
                              &rho[ijk],          \
                              &press[ijk],        \
                              &eps[ijk],          \
                              &vx[ijk],           \
                              &vy[ijk],           \
                              &vz[ijk],           \
                              &BU[0],             \
                              &BU[1],             \
                              &BU[2],             \
                              &entropy[ijk],      \
                              &Y_e[ijk],          \
                              &temperature[ijk]); \
    }

#define GRHAYLMHD_PACK_CONS(struct_name)         \
    ghl_initialize_conservatives(rho_tilde[ijk], \
                                 tau_tilde[ijk], \
                                 S_x_tilde[ijk], \
                                 S_y_tilde[ijk], \
                                 S_z_tilde[ijk], \
                                 ent_tilde[ijk], \
                                 Y_e_tilde[ijk], \
                                 &struct_name);

#define GRHAYLMHD_UNPACK_CONS(struct_name)    \
    ghl_return_conservatives(&struct_name,    \
                             &rho_tilde[ijk], \
                             &tau_tilde[ijk], \
                             &S_x_tilde[ijk], \
                             &S_y_tilde[ijk], \
                             &S_z_tilde[ijk], \
                             &ent_tilde[ijk], \
                             &Y_e_tilde[ijk]);

#define GRHAYLMHD_PACK_METRIC(struct_name) \
    ghl_initialize_metric(alp[ijk],        \
                          betax[ijk],      \
                          betay[ijk],      \
                          betaz[ijk],      \
                          gxx[ijk],        \
                          gxy[ijk],        \
                          gxz[ijk],        \
                          gyy[ijk],        \
                          gyz[ijk],        \
                          gzz[ijk],        \
                          &struct_name);

#define GRHAYLMHD_PACK_METRIC_ENFORCE_DETGAMMAEQ1(struct_name) \
    ghl_enforce_detgtij_and_initialize_ADM_metric(alp[ijk],    \
                                                  betax[ijk],  \
                                                  betay[ijk],  \
                                                  betaz[ijk],  \
                                                  gxx[ijk],    \
                                                  gxy[ijk],    \
                                                  gxz[ijk],    \
                                                  gyy[ijk],    \
                                                  gyz[ijk],    \
                                                  gzz[ijk],    \
                                                  &struct_name);

#define GRHAYLMHD_PACK_EXTRINSIC_CURVATURE(struct_name) \
    ghl_initialize_extrinsic_curvature(kxx[ijk],        \
                                       kxy[ijk],        \
                                       kxz[ijk],        \
                                       kyy[ijk],        \
                                       kyz[ijk],        \
                                       kzz[ijk],        \
                                       &struct_name);

#endif // GRHAYLET_GRHAYLMHD_H
