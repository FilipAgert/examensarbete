module solver_module
    use iso_fortran_env, only: dp=>real64
    use result_module
    
    implicit none
    private
    public :: solver, create_solver, init_solver

    ! Define the base solver class as a derived type
    type, abstract :: solver
        real(dp), allocatable :: matrix(:,:)  ! The matrix is dynamically sized. Matrix should have no connections beyond neigbours (no fusion/fission)
        type(result) :: result
        integer, allocatable :: startIdxs(:)
        integer, allocatable :: fissionIdxs(:)
        integer, allocatable :: fusionIdxs(:)
    contains
        procedure(solve_interface), deferred, public :: solve  ! Define abstract method 'solve'
        procedure, public :: printResult => print_result_method
        procedure, public :: generateConnections => generate_connections_method
        procedure :: fusionFrac => fusion_frac_method
        procedure :: fissionFrac => fission_frac_method
    end type solver

    ! Define an interface for the solve methodse iso_fortran_env, only: dp=>real64
    abstract interface
        subroutine solve_interface(self)
            import :: solver
            class(solver), intent(inout) :: self
        end subroutine solve_interface
    end interface

    ! A generic constructor to initialize a solver object
    interface create_solver
        module procedure init_solver
    end interface

contains

    ! Constructor to initialize the matrix in the solver base class
    subroutine init_solver(self, input_matrix, start_idxs, fusion_idxs, fission_idxs)
        class(solver), intent(inout) :: self
        real(dp), intent(in) :: input_matrix(:,:)
        integer, intent(in) :: start_idxs(:), fusion_idxs(:), fission_idxs(:)
        
        self%matrix = input_matrix  ! Store the matrix
        self%startIdxs = start_idxs
        self%fusionIdxs = fusion_idxs
        self%fissionIdxs = fission_idxs
    end subroutine init_solver

    subroutine print_result_method(self) !Prints result
        class(solver), intent(in) :: self
        call self%result%printResult()
    end subroutine print_result_method

    function generate_connections_method(self, start_idx) result(mat) !generates 25% probability to fusion
        class(solver), intent(in) :: self
        real(dp), allocatable :: mat(:,:)
        integer, intent(in) :: start_idx
        integer, allocatable :: fusion_idxs(:), fission_idxs(:)
        integer :: i
        real :: fusionFrac = 0.25 !Prob to fusion
        real :: fissionFrac = 0.25 !prob to fission
        fusion_idxs = self%fusionIdxs
        fission_idxs = self%fissionIdxs

        mat = self%matrix
        do i = 1, size(fusion_idxs) !For each fusion_idxs, generate connection to start_idx
            mat(:, fusion_idxs(i)) = mat(:, fusion_idxs(i))*(1.0 - fusionFrac)
            mat(start_idx, fusion_idxs(i)) = fusionFrac
        end do

        do i = 1, size(fission_idxs) !For each fusion_idxs, generate connection to start_idx
            mat(:, fission_idxs(i)) = mat(:, fission_idxs(i))*(1.0 - fissionFrac)
            mat(start_idx, fission_idxs(i)) = fissionFrac
        end do        
    end function generate_connections_method

    function fusion_frac_method(self, pd) result (fusionFrac)
        class(solver),intent(in) :: self
        real(dp), intent(in) :: pd(:)
        real :: fusionFrac, fusionProb, fissionProb
        integer :: i
        fusionFrac = 0
        fusionProb = 0
        fissionProb = 0
        do i = 1, SIZE(self%fusionIdxs)
            fusionProb = fusionProb + pd(self%fusionIdxs(i))
        end do

        do i = 1, SIZE(self%fissionIdxs)
            fissionProb = fissionProb + pd(self%fissionIdxs(i))
        end do

        fusionFrac = fusionProb/(fusionProb + fissionProb)
    end function fusion_frac_method

    function fission_frac_method(self, pd) result(fissionFrac)
        class(solver),intent(in) :: self
        real(dp), intent(in) :: pd(:) !probability density
        real :: fissionFrac
        fissionFrac = (1.0 - fusion_frac_method(self, pd))
    end function fission_frac_method

    subroutine add_result_method(self, pd, time, multiplications, startCoord)
        class(solver), intent(inout) :: self
        real(dp), intent(in) :: pd(:)
        real, intent(in) :: time
        integer, intent(in) :: multiplications
        real, intent(in) :: startCoord(:)

        call self%result%addResult(pd, startCoord, time, multiplications, self%fusionFrac(pd), self%fissionFrac(pd))
    end subroutine add_result_method

end module solver_module



module sparse_solver_module
    use solver_module
    use fsparse
    use iso_fortran_env, only: dp=>real64
    implicit none
    public :: sparse_solver

    ! Define a subclass for sparse matrix solver
    type, extends(solver) :: sparse_solver
    contains
        procedure :: solve => sparse_solve_method
        procedure :: init => sparse_init_solver
        procedure :: sparseGenerateConnections => sparse_generate_connections_method
    end type sparse_solver

contains

    
    subroutine sparse_init_solver(self, input_matrix, start_idxs, fusion_idxs, fission_idxs)
        class (sparse_solver), intent(inout) :: self
        real(dp), intent(in) :: input_matrix(:,:)
        integer, intent(in) :: start_idxs(:), fusion_idxs(:), fission_idxs(:)

        
        call init_solver(self, input_matrix, start_idxs, fusion_idxs, fission_idxs) !Call base initializer
    end subroutine sparse_init_solver

    function sparse_generate_connections_method(self, start_idx) result (sparseMat)
        class (sparse_solver), intent(inout) :: self
        integer, intent(in) :: start_idx
        real(dp), allocatable :: denseMat(:,:)
        type(COO_dp) :: sparseMat
        
        denseMat = self%generateConnections(start_idx)
        call dense2coo(denseMat, sparseMat)
    end function sparse_generate_connections_method

    ! Sparse matrix solver method
    subroutine sparse_solve_method(self)
        class(sparse_solver), intent(inout) :: self

    end subroutine sparse_solve_method

    

end module sparse_solver_module