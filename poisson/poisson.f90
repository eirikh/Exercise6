program poisson

      implicit none
#include "pcommon.h"
!==================================================================
!
!     solve the two-dimensional Poisson equation on a unit square 
!     using one-dimensional eigenvalue decompositions
!     and fast sine transforms
!
!     note: n needs to be a power of 2
!
!     einar m. r�nquist
!     ntnu, october 2000
!
!===================================================================
   integer(kind=8), parameter :: n  = 128
   integer(kind=8), parameter :: m  = n-1
   integer(kind=8), parameter :: nn = 4*n
! b is G=(TU + UT) 
! diag are the eigenvalues lambda
 
   real(kind=8), dimension(:,:),  allocatable ::  b,bt 
   integer(kind=8), dimension(:), allocatable ::  m_per_p 
   real(kind=8)    ::  z(0:nn-1)
   real(kind=8)    ::  diag(m)
   real(kind=8)    ::  pi, umax, h
   real(kind=4)    ::  tarray(2), t1, t2, dt
   integer   ::  i,j
   integer   ::  ierror,comm,mpi_size,rank,r,mpi_m,rankp1
   integer   ::  astatus
      
   call second(t1)

   call mpistuff(rank,mpi_size,ierror)
   rankp1 = rank + 1

   allocate (m_per_p(mpi_size),stat=astatus)
   if (astatus /= 0) stop "Not enough memory for m_per_p"

!  set the sizes, with the first CPUs having more columns
   mpi_m = m/mpi_size
   m_per_p = mpi_m
   r= mod(m,mpi_size)
   do i = 1,r 
      m_per_p(i) = m_per_p(i) + 1
   enddo 

   allocate (b(m,m_per_p(rankp1)),stat=astatus)
   if (astatus /= 0) stop "Not enough memory for b"

   allocate (bt(m,m_per_p(rankp1)),stat=astatus)
   if (astatus /= 0) stop "Not enough memory for bt"

   h    = 1./n
   pi   = 4.*atan(1.)

!     full diag needed on all CPUs. Easiest to just calculate it on each
!     CPU
   do i=1,m
      diag(i) = 2*(1-cos(i*pi/n))
   enddo

!  see equation 6, page 6            
   do j=1,m_per_p(rankp1)
      do i=1,m
         b(i,j) = h*h
      enddo
   enddo
      
!  do the first sine transforms, no communication needed
   do j=1,m_per_p(rankp1)
      call fst(b(1,j), n, z, nn)
   enddo 

!  transpose function must be rewritten
   call transp (b, bt, m, m_per_p(rankp1)i,mpi_size)

!  transform back
   do i=1,m_per_p(rankp1)
      call fstinv(bt(1,i), n, z, nn)
   enddo 

!  Divide by diag elements. All elements available on each node.
   do j=1,m_per_p(rankp1)
      do i=1,m
         bt(i,j) = bt(i,j)/(diag(i)+diag(j))
      enddo
   enddo

!  transform again
   do i=1,m_per_p(rankp1)
      call fst (bt(1,i), n, z, nn)
   enddo 

!  transpose again
   call transp (b, bt, m, m_per_p)

!  last back transform
   do j=1,m_per_p(rankp1)
      call fstinv (b(1,j), n, z, nn)
   enddo 

!  now, all cpu will print largest u
!  some I/O required here
   umax = 0.0
   do j=1,m_per_p(rankp1)
      do i=1,m
         if (b(i,j) .gt. umax) umax = b(i,j)
      enddo
   enddo

   write(6,*) ' ' 
   write(6,*) umax

   call mpiendstuff(ierror)
   stop
end
