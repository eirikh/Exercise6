#ifndef COMMON_H_
#define COMMON_H_

#define index(i,j) max(i,j)*(max(i,j)-3)/2 + i + j

#ifdef HAVE_MPI
integer world_comm
common / pcommon / world_comm
#include "mpif.h"
#endif

#endif
