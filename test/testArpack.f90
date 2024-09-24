program testArpack
    use iso_fortran_env, only: dp=>real64
    use fsparse

    implicit none
    external :: dnaupd
    external ::dneupd
    integer :: IDO = 0
    character(len = 1) :: BMAT = 'I'
    integer :: N = 10
    character(len = 2) ::  WHICH = 'LM'
    integer :: NEV = 1 !number of eigenvalues comptued
    double precision :: TOL = -1
    double precision,allocatable :: RESID(:)
    integer :: NCV = 10
    double precision, allocatable :: V(:,:)
    integer :: IPARAM(11)
    integer :: ishfts = 1
    integer :: maxitr = 300
    integer :: mode = 1
    integer :: LDV
    double precision, allocatable :: workd(:), workl(:), workev(:)
    integer :: LWORKL
    integer :: INFO = 0 !0 = randomized initial vector, 1 = resid initial vector
    logical :: converged = .FALSE.
    integer :: IPNTR(14)
    double precision, dimension(:,:), allocatable :: A
    LOGICAL :: RVEC
    character(len = 1) :: HOWMNY = 'A'
    logical, allocatable :: SELECT(:)
    double precision, dimension(:), allocatable :: DR, DI, temp
    double precision, dimension(:,:), allocatable :: Z, SIGMAR, SIGMAI
    TYPE(COO_dp) :: sparseMat
    integer :: LDZ
    integer :: i
    LDZ = N
    allocate(Z(N, NEV + 1))
    allocate(DR(NEV + 1), DI(NEV + 1))
    allocate(SELECT(NCV))
    allocate(temp(N))

    LWORKL = 3*NCV**2 + 6 *NCV
    allocate(workl(LWORKL))
    allocate(workd(3*N), workev(3*NCV))
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
    do i = 1,SIZE(A,1)
        A(i,i) = i
    end do
    call dense2coo(A, sparseMat)
    do while(.not. converged)
        call dnaupd(IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV,IPARAM,IPNTR, workd, workl, lworkl,info)

        if(IDO .eq. -1 .or. IDO .eq. 1) then
            temp = 0
            call matvec(sparseMat, workd(ipntr(1):ipntr(1) + N - 1),temp)
            workd(ipntr(2) : ipntr(2) + N - 1) = temp
            !workd(ipntr(2) : ipntr(2) + N - 1) = matmul(A, workd(ipntr(1) : ipntr(1) + N - 1))
        else 
            converged = .TRUE.
        end if
    end do
    if ( info .lt. 0 ) then
        print *, ' '
        print *, ' Error with _naupd, info = ', info
        print *, ' Check the documentation of _naupd'
        print *, ' '
    else
        RVEC = .TRUE. !Calculate eigenvector.
        call dneupd(RVEC,HOWMNY, SELECT, DR, DI, Z, LDZ, SIGMAR, SIGMAI, WORKEV, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, &
                    LDV, IPARAM, IPNTR, WORKD, WORKL, LWORKL, INFO)


        if(INFO .ne. 0) then
            print *, ' '
            print *, ' Error with _neupd, info = ', INFO
            print *, ' Check the documentation of _neupd. '
            print *, ' '
        else
            print *, ' '
            print *, 'Eigenvalues:', DR
            print* , 'Eigenvector:', Z(:,1)
        endif

    endif



end program testArpack