#ifndef GRHAYLIDX_H_
#define GRHAYLIDX_H_

#include "loop_device.hxx"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "GRHayLib.h"

#define CHECK_PARAMETER(par) if(par==-1) CCTK_VERROR("Please set %s::%s in your parfile",CCTK_THORNSTRING,#par);

#endif // GRHAYLIDX_H_
