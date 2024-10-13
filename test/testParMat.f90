PROGRAM testParMat

    USE potsurf
    USE sparseFromPotSurf
    use iso_fortran_env, only: dp=>real64
    use, intrinsic :: omp_lib
    
    IMPLICIT NONE
  
  
    ! --------------------------------------------------------------------
    ! Input parameters
    ! --------------------------------------------------------------------
  
    INTEGER(kind=i_kind) :: Z = 102, A = 256                                  ! Proton number and mass number of compound nucleus
    REAL(kind=r_kind) :: Eexc = 30.0_r_kind                                  ! Excitation energy relative ground state
    REAL(kind=r_kind), parameter :: Egs = -4.026                              ! Ground-state energy for nucleus (Z,A)=(102,256) 
    INTEGER(kind=i_kind) :: II_fusion = 3                                   ! Fusion occurs if index II becomes less than or equal than II=3
    REAL(kind=r_kind) :: Rneck_fission = 1.5                                  ! Fission occurs if neck radius becomes less than or equal than 1.5 fm
    
    ! --------------------------------------------------------------------
    ! Input files
    ! --------------------------------------------------------------------
  
    CHARACTER(100) :: filename_pot5D = '../input/pot5D102256.dat'                      ! Total potential energy (Emac + Emic)
    CHARACTER(100) :: filename_emac5D = '../input/emacr102256.dat'                     ! Macroscopic potential energy Emac
    CHARACTER(100) :: filename_rneck5D = '../input/neck5D2024-09-11.dat'               ! Neck radii
  
    ! --------------------------------------------------------------------
  
    REAL(kind=r_kind) :: Etotal                                               
    REAL(kind=r_kind) :: prob_test, fusion_prob, fission_prob
    INTEGER(kind=i_kind) :: II_1=19, JJ_1=7, KK_1=5,LL_1=8,MM_1=30 
    INTEGER(kind=i_kind) :: II_2=20, JJ_2=6, KK_2=5,LL_2=8,MM_2=30
  
  
    ! ---------------------------------------------------------------------
    type(COO_dp) :: COO
    type(CSR_dp) :: CSR
    INTEGER(kind=i_kind), dimension(:), allocatable :: fusionIndices, fissionIndices
    real(kind=r_kind), dimension(1) :: energies
    Integer(kind=i_kind), dimension(5) :: startCoord = [17,7,3,13,30]
    INTEGER(kind=i_kind), dimension(1) :: startIndices
    real(kind=r_kind), allocatable :: colVec(:), temp(:)
    integer(8) :: i, nmuls, count1, count2, rate
    
  
    WRITE(*,*) ''
    WRITE(*,*) ''
    WRITE(*,*) '************************************************************'
    WRITE(*,*) '                     Input data                             '
    WRITE(*,*) '************************************************************'
    WRITE(*,*) ''
    WRITE(*,'(A20,A9,I3,A2,I3,A1)') 'Compound nucleus: ', '(Z,A) = (', Z,', ',A,')'
    WRITE(*,'(A20,A8,F7.3,A4)') 'Excitation energy: ', 'Eexc = ', Eexc, ' MeV'
    WRITE(*,*) 'filename_pot5D = ', filename_pot5D
    WRITE(*,*) 'filename_emac5D = ', filename_emac5D
    WRITE(*,*) 'filename_rneck5D = ', filename_rneck5D
    WRITE(*,*) ''
  
    Etotal = Eexc + Egs
   
    WRITE(*,*) ''
    WRITE(*,*) '************************************************************'
    WRITE(*,*) '*************** Potential-energy surface *******************'
    WRITE(*,*) '************************************************************'
  
    
    ! Read potential-energy surface
    ! Stored in matrix Epot(I,J,K,L,M)
    CALL read_5d_epot(Z, A, filename_pot5D)
    
    ! Read macroscopic energy
    ! Stored in matrix Emac(I,J,K,L,M)
    CALL read_5d_emac(Z, A, filename_emac5D) 
  
      ! Read neck radii
    ! Stored in matrix Rneck(I,J,K,L,M)
    CALL read_5d_rneck(Z,A,filename_rneck5D)
  
    fusion_prob = 0.5_r_kind
    fission_prob = 0.5_r_kind
    !Set probability to fusion/fission if landed on fusion/fission grid point
  
    nmuls = 50
    !Generate sparse matrix.
    COO = COOfromPotSurf(Z,A,Etotal, II_fusion, fusion_prob, Rneck_fission, fission_prob)
    allocate(colVec(COO%nrows), temp(COO%nrows))

    colVec = 1.0_r_kind/COO%nrows
    CALL system_clock(count_rate = rate)
    CALL system_clock(count=count1)
    do i = 1,nmuls
        temp = 0
        call matvec(COO, colVec,temp)
    end do
    CALL system_clock(count=count2)
    print*, nmuls, " matrix-vector multiplication in series done in ", real((count2-count1))/real(rate), " s"

    call coo2ordered(COO, .TRUE.)
    call coo2csr(COO, CSR)

    !Initialize column vector:
    
   
    
    do i = 1,nmuls
        temp = 0
        

    
        call matrixVec(CSR, colVec, temp)
    end do
    CALL system_clock(count=count1)
    print*, nmuls, " matrix-vector multiplication in parallel (4 threads) done in ", real((count1-count2))/real(rate), " s"

    
  
    contains
        subroutine matrixVec(matrix,vec_x,vec_y)
            

            implicit none
            type(CSR_dp), intent(in) :: matrix
            real(dp), intent(in)    :: vec_x(:)
            real(dp), intent(inout) :: vec_y(:)
            real(dp) :: vec_x_private(size(vec_x))
            ! real(dp) :: dat(matrix%nnz)
            ! integer :: rowptr(matrix%nrows + 1), col(matrix%nnz)
            integer :: i, j
            real(dp) :: aux
            call omp_set_num_threads(4)
            ! dat = matrix%data
            ! col = matrix%col
            ! rowptr = matrix%rowptr
            associate(data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, nnz => matrix%nnz, nrows => matrix%nrows, &
                &ncols => matrix%ncols, sym => matrix%sym )
                if( sym == k_NOSYMMETRY) then
                    !$omp parallel shared(matrix, vec_y, vec_x) private(i,j)

                    !$omp do schedule(static)
                    
                    do i = 1, nrows
                        !print*, "i: ", i
                        !print '("Thread: ", i0, " handling i = ", i0)', omp_get_thread_num(), i
                        ! if(mod(i, 100000) == 0) then
! 
                            ! print '("Thread: ", i0, " handling i = ", i0)', omp_get_thread_num(), i
                        ! endif
                        do j = matrix%rowptr(i), matrix%rowptr(i+1)-1
                            
                            vec_y(i) = vec_y(i) + matrix%data(j) * vec_x(matrix%col(j))
                        end do
                    end do
                    !$omp end do
                    !$omp end parallel
                endif
            end associate
        end subroutine matrixVec
        
END PROGRAM testParMat
  
