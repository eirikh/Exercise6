subroutine mpistuff(rank,mpi_size,ierror)
   implicit none

#include "pcommon.h"

   integer rank, mpi_size
   integer ierror

#ifdef HAVE_MPI
    call mpi_init(ierror)
    call mpi_comm_size(mpi_comm_world,mpi_size,ierror)
    call mpi_comm_rank(mpi_comm_world,rank,ierror)
    call mpi_comm_dup(mpi_comm_world,world_comm,ierror)
#else
    rank = 0
    mpi_size = 1
#endif
end

         
subroutine mpiendstuff(ierror)
   implicit none

   integer ierror

#ifdef HAVE_MPI
   call mpi_finalize(ierror)
#endif
end


subroutine transp (at, a, m, mp)
   implicit none
!====================================================
!     set at equal to the transpose of a 
!====================================================
   real(kind=8)   :: a(m,mp)
   integer(kind=8), intent(in) :: m, mp
   integer i,j
#ifdef HAVE_MPI
   real(kind=8)   :: at(m*mp)
#else
   real(kind=8)   :: at(m,mp)
#endif

#ifdef HAVE_MPI

   do i = 1,m
      
   enddo

#else

do j=1,m
   do i=1,m
      at(j,i) = a(i,j)
   enddo
enddo

#endif
   return
end
