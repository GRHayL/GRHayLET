#ifndef GRHAYLIDX_H_
#define GRHAYLIDX_H_

#include "loop_device.hxx"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "GRHayLib.h"

/* GRHayL's API is double precision. */
#ifndef CCTK_REAL_PRECISION_8
#error "GRHayLIDX requires Cactus REAL_PRECISION=8 (double)"
#endif

#define CHECK_PARAMETER(par) if(par==-1) CCTK_VERROR("Please set %s::%s in your parfile",CCTK_THORNSTRING,#par);

#endif // GRHAYLIDX_H_
