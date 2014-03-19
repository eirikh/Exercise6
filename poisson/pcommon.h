#ifndef COMMON_H_
#define COMMON_H_

#ifdef HAVE_MPI
#include "mpif.h"
integer world_comm
common /pcommon  /world_comm
#endif

#endif
