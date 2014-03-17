      program poisson

      implicit none
      include 'poisson.h'
#ifdef HAVE_MPI
      include 'mpif.h'
#endif
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
      integer(kind=8), parameter :: n  = 128
      integer(kind=8), parameter :: m  = n-1
      integer(kind=8), parameter :: nn = 4*n
c b is G=(TU + UT) 
c diag are the eigenvalues lambda
c 
      real(kind=8), dimension(:,:),  allocatable ::  b,bt 
      integer(kind=8), dimension(:), allocatable ::  m_per_p 
      real(kind=8)    ::  z(0:nn-1)
      real(kind=8)    ::  diag(m)
      real(kind=8)    ::  pi, umax, h
      real(kind=4)    ::  tarray(2), t1, t2, dt
      integer   ::  i,j
      integer   ::  ierror,comm,mpi_size,rank,r,mpi_m,rankp1
      integer   ::  astatus
      
      write(*,*) 'lalal'
      call second(t1)

      write(*,*) 'jijijij'
      call mpistuff(rank,mpi_size,ierror)
      rankp1 = rank + 1

      allocate (m_per_p(mpi_size),stat=astatus)
      if (astatus /= 0) stop "Not enough memory for m_per_p"

c     set the sizes, with the first CPUs having more columns
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

c     full diag needed on all CPUs. Easiest to just calculate it on each
c     CPU
      do i=1,m
         diag(i) = 2*(1-cos((rank*mpi_m+i)*pi/n))
      enddo

c     see equation 6, page 6            
      do j=1,m_per_p(rankp1)
         do i=1,m
            b(i,j) = h*h
         enddo
      enddo
      
c     do the first sine transforms, no communication needed
      do j=1,m_per_p(rankp1)
         call fst(b(1,j), n, z, nn)
      enddo 

c     transpose function must be rewritten
      call transp(bt, b, m)

c     transform back
      do i=1,m_per_p(rankp1)
         call fstinv(bt(1,i), n, z, nn)
      enddo 

c     Divide by diag elements. All elements available on each node.
      do j=1,m_per_p(rankp1)
         do i=1,m
            bt(i,j) = bt(i,j)/(diag(i)+diag(j))
         enddo
      enddo

c     transform again
      do i=1,m_per_p(rankp1)
         call fst (bt(1,i), n, z, nn)
      enddo 

c     transpose again
      call transp (b, bt, m)

c     last back transform
      do j=1,m_per_p(rankp1)
         call fstinv (b(1,j), n, z, nn)
      enddo 

c     now, all cpu will print largest u
c     some I/O required here
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
c
      subroutine mpistuff(rank,mpi_size,ierror)
      implicit none
c
#ifdef HAVE_MPI
      include 'mpif.h'
#endif
c
         integer rank, mpi_size
         integer ierror, world_comm
c
         write(*,*) 'blusdasasdf'
#ifdef HAVE_MPI
c
         write(*,*) 'bluggi'
         call mpi_init(ierror)
         call mpi_comm_size(mpi_comm_world,mpi_size,ierror)
         call mpi_comm_rank(mpi_comm_world,rank,ierror)
c         call mpi_comm_dup(mpi_comm_world,world_comm,ierror)
         write(*,*) 'dkashfkasdfk'
#else
c         write(*,*) 'kakskkkas'
c         rank = 0
c         mpi_size = 1
#endif
c
c
      end
         
      subroutine mpiendstuff(ierror)
      implicit none
c
         integer ierror
c
#ifdef HAVE_MPI
      call mpi_finalize(ierror)
#endif
      end

      subroutine transp (at, a, m)
      implicit none
c====================================================
c     set at equal to the transpose of a 
c====================================================
      real(kind=8)   :: a(m,m), at(m,m)
      integer(kind=8), intent(in) :: m
      integer i,j

      do j=1,m
         do i=1,m
            at(j,i) = a(i,j)
         enddo
      enddo
      return
      end



