#define index(i,j) max(i,j)*(max(i,j)-3)/2 + i + j

#ifdef HAVE_MPI
#include "mpif.h"
integer world_comm
common / pcommon / world_comm
#endif


#ifdef HAVE_OPENMP
#include "omp_lib.h"
#endif
