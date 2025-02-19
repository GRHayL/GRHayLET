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

// clang-format off
#define OBLOOPX(_dir)                                                    \
    if(cctk_bbox[_dir]) {                                                \
        for(int k = 0; k < cctk_lsh[2]; k++) {                           \
            for(int j = 0; j < cctk_lsh[1]; j++) {                       \
                const int i_dst = (_dir == GRHayLMHD_xmin) * imin        \
                                + (_dir == GRHayLMHD_xmax) * imax;       \
                const int i_src = (_dir == GRHayLMHD_xmin) * (imin + 1)  \
                                + (_dir == GRHayLMHD_xmax) * (imax - 1); \
                const int src_idx = CCTK_GFINDEX3D(cctkGH, i_src, j, k); \
                const int dst_idx = CCTK_GFINDEX3D(cctkGH, i_dst, j, k);

#define OBLOOPY(_dir)                                                    \
    if(cctk_bbox[_dir]) {                                                \
        for(int k = 0; k < cctk_lsh[2]; k++) {                           \
            for(int i = 0; i < cctk_lsh[0]; i++) {                       \
                const int j_dst = (_dir == GRHayLMHD_zmin) * jmin        \
                                + (_dir == GRHayLMHD_zmax) * jmax;       \
                const int j_src = (_dir == GRHayLMHD_ymin) * (jmin + 1)  \
                                + (_dir == GRHayLMHD_ymax) * (jmax - 1); \
                const int src_idx = CCTK_GFINDEX3D(cctkGH, i, j_src, k); \
                const int dst_idx = CCTK_GFINDEX3D(cctkGH, i, j_dst, k);

#define OBLOOPZ(_dir)                                                    \
    if(cctk_bbox[_dir]) {                                                \
        for(int j = 0; j < cctk_lsh[1]; j++) {                           \
            for(int i = 0; i < cctk_lsh[0]; i++) {                       \
                const int k_dst = (_dir == GRHayLMHD_zmin) * kmin        \
                                + (_dir == GRHayLMHD_zmax) * kmax;       \
                const int k_src = (_dir == GRHayLMHD_zmin) * (kmin + 1)  \
                                + (_dir == GRHayLMHD_zmax) * (kmax - 1); \
                const int src_idx = CCTK_GFINDEX3D(cctkGH, i, j, k_src); \
                const int dst_idx = CCTK_GFINDEX3D(cctkGH, i, j, k_dst);

#define ENDLOOP }}}
// clang-format on


#endif // GRHAYLET_GRHAYLMHD_H
