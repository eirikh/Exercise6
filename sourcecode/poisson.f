      program poisson

      implicit none
      include 'mpif.h'
c==================================================================
c
c     solve the two-dimensional Poisson equation on a unit square 
C     using one-dimensional eigenvalue decompositions
c     and fast sine transforms
c
c     note: n needs to be a power of 2
c
c     einar m. rønquist
c     ntnu, october 2000
c
c===================================================================
      integer*8,parameter (n  = 128)
      integer*8,parameter (m  = n-1)
      integer*8,parameter (nn = 4*n)
c b is G=(TU + UT) 
c diag are the eigenvalues lambda
c 
      real*8,dimension(:),   allocatable ::  diag
      real*8,dimension(:,:), allocatable ::  b,bt 
      real*8,dimension(:),   allocatable ::  z
      real*8,dimension(:),   allocatable ::  m_per_p 
      real*8      pi,wstart,wstop
      real*4      tarray(2), t1, t2, dt
      integer*8   ierror,comm,mpi_size,rank,r,mpi_m
      
      call second(wstart)

      call mpi_init(ierror)
      call mpi_comm_size(mpi_comm_world,mpi_size,ierror)
      call mpi_comm_rank(mpi_comm_world,rank,ierror)

      mpi_m = m/mpi_size
      m_per_p = mpi_m
      r= mod(m,mpi_size)
      do i = 1,r 
         m_per_p(i) = m_per_p(i) + 1
      enddo 

      h    = 1./n
      pi   = 4.*atan(1.)

      do i=1,m
         diag(i) = 2*(1-cos((rank*mpi_m+i)*pi/n))
      enddo
c see equation 6, page 6            
      do j=1,m_per_p(rank)
         do i=1,m
            b(i,j) = h*h
         enddo
      enddo
      
      do j=1,m_per_p(rank)
         call fst(b(1,j), n, z, nn)
      enddo 

      call transp(bt, b, m)
      do i=1,m_per_p(rank)
         call fstinv(bt(1,i), n, z, nn)
      enddo 

      do j=1,m_per_p(rank)
         do i=1,m
            bt(i,j) = bt(i,j)/(diag(i)+diag(j))
         enddo
      enddo

      do i=1,m
         call fst (bt(1,i), n, z, nn)
      enddo 
      call transp (b, bt, m)
      do j=1,m
         call fstinv (b(1,j), n, z, nn)
      enddo 

      umax = 0.0
      do j=1,m
         do i=1,m
            if (b(i,j) .gt. umax) umax = b(i,j)
         enddo
      enddo

      write(6,*) ' ' 
      write(6,*) umax

      stop
      end

      subroutine transp (at, a, m)
c====================================================
c     set at equal to the transpose of a 
c====================================================
      real*8 a(m,m), at(m,m)

      do j=1,m
         do i=1,m
            at(j,i) = a(i,j)
         enddo
      enddo
      return
      end



