program main
  use iso_fortran_env, only: dp=>real64
  use fsparse

  implicit none
  type (COO_dp) :: sparse
  external :: dnaupd
  external ::dneupd
  integer :: IDO = 0
  character(len = 1) :: BMAT = 'I'
  integer :: N = 2
  character(len = 2) ::  WHICH = 'LM'
  integer :: NEV = 1 !number of eigenvalues comptued
  double precision :: TOL = 1e-3
  double precision,allocatable :: RESID(:)
  integer :: NCV = 3
  double precision, allocatable :: V(:,:)
  integer :: IPARAM(11)
  integer :: ishfts = 1
  integer :: maxitr = 300
  integer :: mode = 1
  integer :: LDV
  double precision, allocatable :: workd(:), workl(:)
  integer :: LWORKL
  integer :: INFO = 1
  logical :: converged = .FALSE.
  integer :: IPNTR(14)
  double precision, dimension(:,:), allocatable :: A
  LOGICAL :: RVEC
  integer :: i
  character(len = 1) :: HOWMNY = 'A'
  logical, allocatable :: SELECT(:)
  allocate(SELECT(NCV))

  LWORKL = 3*NCV**2 + 6 *NCV
  allocate(workl(LWORKL))
  allocate(workd(3*N))
  allocate(RESID(N))
  allocate(A(N,N))
  allocate(V(N,NCV))
  IPARAM = 0
  IPARAM(1) = ishfts
  IPARAM(3) = maxitr !Number of arnoldi update iterations
  IPARAM (4) = 1 !needs to be 1
  IPARAM(7) = mode !specifies eigenvalue problem A*x = lambda*x

  
  LDV = N
  RESID = [4.0, 3.0]
  A = 0
  A(1,1) = 1
  A(2,2) = 2

  
  call dense2coo(A, sparse)

!  do while(converged)
!      call dnaupd(IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV,IPARAM,IPNTR, workd, workl, lworkl,info)!
!
!      if(IDO == -1 .or. IDO == 1) then
!          workd(ipntr(2) : ipntr(2) + N - 1) = matmul(A, workd(ipntr(1) : ipntr(1) + N - 1))
!      else 
!          converged = .TRUE.
!      end if
!  end do
!  if ( info .lt. 0 ) then
!      print *, ' '
!      print *, ' Error with _naupd, info = ', info
!      print *, ' Check the documentation of _naupd'
!      print *, ' '
!  else
!      RVEC = .TRUE. !Calculate eigenvector.!
!
!  end if



end program main