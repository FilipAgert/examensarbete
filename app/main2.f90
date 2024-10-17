PROGRAM main

    USE potsurf
    use markovSolver, only: results, setupSolver, solveAllEnergies
    
    
    IMPLICIT NONE
  
  
    ! --------------------------------------------------------------------
    ! Input parameters
    ! --------------------------------------------------------------------
  
    INTEGER(kind=i_kind) :: Z = 102, A = 256                                ! Excitation energy relative ground state
    REAL(kind=r_kind), parameter :: Egs = -4.026                              ! Ground-state energy for nucleus (Z,A)=(102,256) 
    INTEGER(kind=i_kind) :: II_fusion = 3                                   ! Fusion occurs if index II becomes less than or equal than II=3
    REAL(kind=r_kind) :: Rneck_fission = 1.5                                  ! Fission occurs if neck radius becomes less than or equal than 1.5 fm
    REAL(kind=r_kind) :: TOL
    INTEGER(kind=i_kind) :: NUM_THREADS = 4
    INTEGER(kind=i_kind) :: NCV
    real(kind=r_kind), dimension(15) :: energies = [13.0, 14.0, 15.0,17.0,20.0,22.0,25.0,30.0,35.0,40.0,45.0,50.0,55.0,60.0,65.0] !excitation energies relative gnd state
    real(kind=r_kind), dimension(11) :: e3s = [13.0,13.1,13.2,13.3,13.4,13.5,13.6,13.7,13.8,13.9,14.0]
    Integer(kind=i_kind), dimension(5,15) :: startCoords = RESHAPE( (/ &
    21, 6, 3, 14, 29,        20, 6, 3, 14, 29,& !13, 14
    20, 6, 3, 14, 29,        20, 6, 3, 14, 29,& !15, 17
    19, 6, 3, 13, 30,        19, 6, 3, 13, 30,& !20, 22
    18, 7, 3, 14, 29,        17, 7, 3, 13, 30,& !25, 30
    17, 7, 3, 13, 30,        16, 7, 3, 13, 30,& !35, 40
    16, 7, 3, 13, 30,        16, 7, 3, 13, 30,&
    15, 7, 3, 14, 30,        15, 7, 3, 14, 30,&
    15, 7, 3, 14, 30 /), (/5,15/) )

    real(kind=r_kind), dimension(1) :: E3
    integer(kind=i_kind), dimension(5,11) :: C3
    LOGICAL :: useFullMM = .FALSE.
    integer :: i
    
    ! --------------------------------------------------------------------
    ! Input files
    ! --------------------------------------------------------------------
  
    !MAX_MM = 30
    !CHARACTER(100) :: filename_pot5D = '../input/old_potential_surface/pot5D102256.dat'                      ! Total potential energy (Emac + Emic)
    !CHARACTER(100) :: filename_emac5D = '../input/old_potential_surface/emacr102256.dat'                     ! Macroscopic potential energy Emac
    !CHARACTER(100) :: filename_rneck5D = '../input/old_potential_surface/neck5D.dat'               ! Neck radii


    !MAX_MM = 40
    CHARACTER(100) :: filename_pot5D = '../input/pot5D102256.dat'                      ! Total potential energy (Emac + Emic)
    CHARACTER(100) :: filename_emac5D = '../input/emacr102256.dat'                     ! Macroscopic potential energy Emac
    CHARACTER(100) :: filename_rneck5D = '../input/neck5D.dat'               ! Neck radii
  
    ! --------------------------------------------------------------------
  
    REAL(kind=r_kind) :: Etotal                                          
  
  
    ! ---------------------------------------------------------------------
  
    WRITE(*,*) ''
    WRITE(*,*) ''
    WRITE(*,*) '************************************************************'
    WRITE(*,*) '                     Input data                             '
    WRITE(*,*) '************************************************************'
    WRITE(*,*) ''
    WRITE(*,'(A20,A9,I3,A2,I3,A1)') 'Compound nucleus: ', '(Z,A) = (', Z,', ',A,')'
    WRITE(*,*) 'filename_pot5D = ', filename_pot5D
    WRITE(*,*) 'filename_emac5D = ', filename_emac5D
    WRITE(*,*) 'filename_rneck5D = ', filename_rneck5D
    WRITE(*,*) ''
  
    energies = energies + Egs !Total energy.
    e3s = e3s + Egs
    TOL = 6.9e-9
    NCV = 15

    E3 = energies(1)
    do i = 1,5
      C3(:,i) = startCoords(:,1)
    end do
    do i = 6,11
      C3(:,i) = startCoords(:,2)
    end do
    call setupSolver(TOL,NCV,NUM_THREADS,Z,A, Rneck_fission, II_fusion, Egs, energies, startCoords ,useFullMM, filename_emac5D,&
                    filename_pot5D, filename_rneck5D)

    call solveAllEnergies()
    call results%printResult()
    call results%printResultToFile()


  
    
    
  END PROGRAM MAIN
  