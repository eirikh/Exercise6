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
#include "pcommon.h"

   integer ierror

#ifdef HAVE_MPI
   call mpi_finalize(ierror)
#endif
end


subroutine transp(at, a, m, mp, mpi_size, rank, ierror)
   implicit none
#include "pcommon.h"

!====================================================
!     set at equal to the transpose of a 
!====================================================
   integer i,j,k, mpi_size, rank, ierror
   integer(kind=8), intent(in) :: m
   integer(kind=8), intent(in) :: mp(mpi_size)
#ifdef HAVE_MPI
   integer :: koff, joff
   integer :: group(mpi_size), offset(mpi_size)
   real(kind=8)   :: a(m*mp(rank+1))
   real(kind=8)   :: at(m*mp(rank+1))
#else
   real(kind=8)   :: a(m,mp(rank+1))
   real(kind=8)   :: at(m,mp(rank+1))
#endif

#ifdef HAVE_MPI
! Packs a into sendbuf at
   koff = 0
   joff = 0
   offset(1)=0
   do i = 1,mpi_size
      group(i) = mp(i)*mp(rank+1)
      if (i .ne. 1) offset(i) = offset(i-1) + group(i-1)
      do j = 1, mp(rank+1)
         do k = 1, mp(i)
            at(koff+k) = a(joff+(j-1)*m+k)
         enddo
         koff = koff + mp(i)
      enddo
      joff = joff + mp(i)
   enddo
   call mpi_barrier(world_comm,ierror)
! sends at to a
   call mpi_alltoallv(at,group,offset,mpi_double_precision,a,group,offset,mpi_double_precision,world_comm,ierror) 
! unwraps the received a into at correctly
   koff = 0
!   if (rank .eq. 0) then
!      write(*,*)
!      write(*,*) "a in"
!      write(*,"(F8.5)") a
!      write(*,*)
 !  endif
!   endif
! loop over all the rows
      do j = 1,mp(rank+1)
! loop over the columns owned by this proc
         do k = 1,m
!            if (rank .eq. 0) write(*,*) "koff + k,k+(j-1)*mp(rank+1)",koff +k,j + (k-1)*mp(rank+1)
!            if (rank .eq. 0) write(*,*)
            at(koff+k) = a(j+(k-1)*mp(rank+1))
         enddo
         koff = koff + m
      enddo
!   if (rank .eq. 3) then
!      write(*,*) "at in"
!      write(*,"(1F8.5)") at
!      do i=1,7
!      write(*,"(2F8.5)") at(i), at(i+7)
!      enddo
!      write(*,*)
!   endif
#else
   do j=1,m
      do i=1,m
         at(j,i) = a(i,j)
      enddo
   enddo

#endif
   return
end



subroutine write_matrix(a,m, mp, mpi_size, rank, ierror)
   implicit none
#include "pcommon.h"

integer mpi_size, rank, ierror
integer(kind=8) :: m
integer(kind=8) :: mp(mpi_size)
real(kind=8) :: a(m*mp(rank+1))
#ifdef HAVE_MPI
integer i,filehand, procsize,stat(2)
integer(kind=mpi_offset_kind) offset, totsize

filehand=4
totsize = m*m*sizeof(a(1))

offset = 0
do i = 1,rank
   offset = offset + mp(i)*m*sizeof(a(1))
enddo
procsize = mp(rank+1)*m

call mpi_file_open(world_comm,"matrix_output.dat",mpi_mode_wronly + mpi_mode_create,mpi_info_null,filehand,ierror)

call mpi_file_set_size(filehand,totsize,ierror)
!call mpi_file_set_view(filehand,offset,mpi_double_precision,mpi_double_precision,"natve",mpi_info_null,ierror)
call mpi_file_write_at_all(filehand,offset,a,procsize,mpi_double_precision,stat,ierror)
call mpi_file_close(filehand,ierror)
#endif
   return
end
