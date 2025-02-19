#ifndef GRHAYLET_GRHAYLMHD_H
#define GRHAYLET_GRHAYLMHD_H

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <math.h>
#include <stdbool.h>

#define LOOP3D(_imin, _imax, _jmin, _jmax, _kmin, _kmax)             \
    for(int k = _kmin; k < _kmax; k++) {                             \
        for(int j = _jmin; j < _jmax; j++) {                         \
            for(int i = _imin; i < _imax; i++) {                     \
                const int ijk  = CCTK_GFINDEX3D(cctkGH, i, j, k);    \
                const int ijkx = CCTK_GFINDEX4D(cctkGH, i, j, k, 0); \
                const int ijky = CCTK_GFINDEX4D(cctkGH, i, j, k, 1); \
                const int ijkz = CCTK_GFINDEX4D(cctkGH, i, j, k, 2);

#define OMPLOOP3D(_imin, _imax, _jmin, _jmax, _kmin, _kmax) \
    _Pragma("omp parallel for")                             \
        LOOP3D(_imin, _imax, _jmin, _jmax, _kmin, _kmax)

#define ENDLOOP3D \
    }             \
    }             \
    }

#endif // GRHAYLET_GRHAYLMHD_H
