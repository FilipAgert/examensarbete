PROGRAM main

  USE potsurf
  USE sparseFromPotSurf
  use sparse_solver_arnoldi_module
  use sparse_solver_module
  use denseMatrix, only: printMatrixToFileDp
  
  
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
  type(sparse_arnoldi_solver) :: solv
  INTEGER(kind=i_kind), dimension(:), allocatable :: fusionIndices, fissionIndices
  real(kind=r_kind), dimension(1) :: energies
  Integer(kind=i_kind), dimension(5) :: startCoord = [17,7,3,13,30]
  INTEGER(kind=i_kind), dimension(1) :: startIndices
  real(kind=r_kind), allocatable :: plotGrid(:,:), printPd(:,:)
  

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


  !Generate sparse matrix.
  COO = COOfromPotSurf(Z,A,Etotal, II_fusion, fusion_prob, Rneck_fission, fission_prob)

  !Gets all linear indexes for the fusion grid points
  fusionIndices = getFusionIndices(II_fusion)
  print*, "Number of fusion indices: ", SIZE(fusionIndices)

  !Gets all linear indexes for the fission grid points
  fissionIndices = getFissionIndices(Rneck_fission)
  print*, "Number of fission indicies: ", SIZE(fissionIndices)
  print*, "Rneck at start: ", Rneck(startCoord(1), startCoord(2), startCoord(3), startCoord(4), startCoord(5)), " fm"
  print*, "Rneck lim: ", Rneck_fission, " fm"

  energies = Eexc !!Excitation energy above ground state
  startIndices(1) = linearIdxFromCoordShifted(startCoord, dimSize)
  call solv%init(COO, real(grid(),4), startIndices, real(energies,4), fusionIndices, fissionIndices&
    , dimSize, fusion_prob, fission_prob)
  call solv%solve()
  call solv%printResult()
  allocate(printPd(COO%ncols,1))
  printPd(:,1) = solv%result%getProbabilityDensity(1)
  call printMatrixToFileDp('PD', printPd)
  plotGrid = gridFromColumnVectorSliced(solv%result%getProbabilityDensity(1), dimSize, startCoord(3), startCoord(4), startCoord(5))
  call printMatrixToFileDp('Fullsize', plotGrid)


  
  
END PROGRAM MAIN
