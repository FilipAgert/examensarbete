module solver_module
    use iso_fortran_env, only: dp=>real64
    use result_module
    use denseMatrix, only: coordFromLinearIdx
    
    implicit none
    private
    public :: solver, create_solver, init_solver

    ! Define the base solver class as a derived type
    type, abstract :: solver
        real(dp), allocatable :: matrix(:,:)  ! The sparseMat is dynamically sized. sparseMat should have no connections beyond neigbours (no fusion/fission)
        type(convergedResult) :: result
        integer, allocatable :: startIdxs(:)
        integer, allocatable :: fissionIdxs(:)
        integer, allocatable :: fusionIdxs(:)
        integer, allocatable :: dimSize(:)
    contains
        procedure(solve_interface), deferred, public :: solve  ! Define abstract method 'solve'
        procedure, public :: printResult => print_result_method
        procedure, public :: generateConnections => generate_connections_method
        procedure :: fusionFrac => fusion_frac_method
        procedure :: fissionFrac => fission_frac_method
        procedure, public :: addResult => add_result_method
        procedure, public :: numberOfGridPoints => get_number_of_coordinates
        procedure, public :: idxToCoords => idx_to_coords_method
        procedure, public :: startingGuess => starting_guess_method
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
    function idx_to_coords_method(self, idx) result(coordinate)
        real, allocatable, dimension(:) :: coordinate !Coordinate in variable space
        integer :: idx
        integer, allocatable, dimension(:) :: coordinateIdx  !coordinate in index space
        class(solver) :: self
        allocate(coordinate(SIZE(self%dimSize)))
        allocate(coordinateIdx(SIZE(self%dimSize)))
        coordinateIdx = coordFromLinearIdx(idx, self%dimSize)

        coordinate = coordinateIdx!Change this to conversion from idx -> coordinate
    end function idx_to_coords_method

    function get_number_of_coordinates(self) result(coordinateN)
        class(solver) :: self
        integer :: coordinateN
        coordinateN = SIZE(self%matrix,1)
    end function get_number_of_coordinates
    ! Constructor to initialize the sparseMat in the solver base class
    subroutine init_solver(self, input_matrix, start_idxs, fusion_idxs, fission_idxs, dim_size)
        class(solver), intent(inout) :: self
        real(dp), intent(in) :: input_matrix(:,:)
        integer, intent(in) :: start_idxs(:), fusion_idxs(:), fission_idxs(:), dim_size(:)
        
        self%matrix = input_matrix  ! Store the sparseMat
        self%startIdxs = start_idxs
        self%fusionIdxs = fusion_idxs
        self%fissionIdxs = fission_idxs
        self%dimSize = dim_size
    end subroutine init_solver

    subroutine print_result_method(self) !Prints result
        class(solver), intent(in) :: self
        call self%result%printResult()
    end subroutine print_result_method

    function starting_guess_method(self) result(guess) !Guesses the initial probability distribution
        class(solver), intent(inout) :: self
        real(dp), allocatable :: guess(:), prevGuess(:)
        integer :: storedResults
        allocate(guess(self%numberOfGridPoints()))
        allocate(prevGuess(self%numberOfGridPoints())) !Allocate array to correct size
        guess = 0
        storedResults = self%result%numResults

        if(storedResults > 0) then !Use previous prob distribution as initial guess
            prevGuess = self%result%getProbabilityDensity(storedResults) !Index latest guess
            guess = prevGuess
        else
            guess = 1.0_dp/self%numberOfGridPoints() !even distribution as first guess
        end if
    end function starting_guess_method

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

module dense_solver_module
    use solver_module
    use iso_fortran_env, only: dp => real64
    use markovSparse, only: norm, convergence
    implicit none

    type, extends(solver) :: dense_solver
    contains
        procedure :: solve => dense_solve_method
        procedure :: runUntilConverged => run_until_converged_method
        procedure :: timeStep => time_step_method
        procedure :: init => dense_init_solver
    end type dense_solver
    contains

    subroutine dense_init_solver(self, denseMat, start_idxs, fusion_idxs, fission_idxs, dimSize)
        class(dense_solver), intent(inout) :: self
        real(dp), dimension(:,:), intent(inout) :: denseMat
        integer, intent(inout) :: start_idxs(:), fusion_idxs(:), fission_idxs(:), dimSize(:)
        call create_solver(self, denseMat, start_idxs, fusion_idxs, fission_idxs, dimSize)
    end subroutine dense_init_solver
    
    subroutine dense_solve_method(self)
        class(dense_solver), intent(inout) :: self
        real(dp), dimension(:,:), allocatable :: denseMat
        integer :: i, startIdx, matMultiplications
        real :: T1, T2, elapsedTime, fusionFrac, fissionFrac
        real(dp), allocatable :: startPd(:), endPd(:)
        allocate(startPd(self%numberOfGridPoints()))
        allocate(endPd(self%numberOfGridPoints()))
        allocate(denseMat(self%numberOfGridPoints(), self%numberOfGridPoints()))
        startPd = 0.0_dp
        endPd = 0.0_dp

        do i = 1, SIZE(self%startIdxs)
            call cpu_time(T1)
            startIdx = self%startIdxs(i)
            denseMat = self%generateConnections(startIdx) !Connects fission and fusion points to starting index
            startPd = self%startingGuess() !starting guess
            call self%runUntilConverged(denseMat, startPd, endPd, matMultiplications)
            call cpu_time(T2)

            fusionFrac = self%fusionFrac(endPd)
            fissionFrac = self%fissionFrac(endPd)
            elapsedTime = T2-T1
            call self%addResult(endPd, elapsedTime, matMultiplications, self%idxToCoords(startIdx))
        end do
    end subroutine dense_solve_method

    function time_step_method(self, matrix, probabilityCoeffs, steps)result (timeStep) !This computes repeated matrix multiplication (A*(A*...*(A*V)
        class(dense_solver), intent(inout) :: self
        integer, intent(in) :: steps                     !In order to only have to store one matrix
        double precision, dimension(:), intent(in) :: probabilityCoeffs
        double precision, dimension(SIZE(probabilityCoeffs)) :: timeStep
        double precision, dimension(SIZE(probabilityCoeffs), SIZE(probabilityCoeffs)), intent(in):: matrix 
        integer :: i
        timeStep = probabilityCoeffs
        do i = 1, steps
            timeStep = MATMUL(matrix, timeStep)
        end do
        
    end function time_step_method
    subroutine run_until_converged_method(self, denseMat, startPd, endPd, numberOfMultiplications)
        class(dense_solver), intent(inout) :: self
        real(dp), dimension(:,:), intent(in) :: denseMat
        real(dp), dimension(:), intent(in) :: startPd
        real(dp), dimension(:), intent(inout) :: endPd
        real(dp), dimension(SIZE(startPd)) :: prevPd, diff
        integer, intent(inout) :: numberOfMultiplications
        real(dp) :: tol, normVal
        integer :: multiplicationSteps
        logical :: converged
        character (len = 2):: normType = "L2"

        numberOfMultiplications = 0
        tol = 1.0/(1e2*SIZE(startPd)) !This tolerance is one part in one hundered assuming even spread of probability.
        multiplicationSteps = 20
        tol = tol * multiplicationSteps / 10 !Dynamically change tolerance based on number of sparseMat multiplications in a row
        print*, tol
        
        converged = .FALSE.
        prevPd = startPd
        do while (.not. converged)
            endPd = self%timeStep(denseMat, prevPd, multiplicationSteps) !Dont check for convergence after every step, rather take a few time steps at a time.
            converged = convergence(endPd, prevPd, tol, normType)
            
            if(mod(numberOfMultiplications, 500) == 0) then

                diff = endPd - prevPd
                normVal = norm(diff, normType)
                print *, 'Matrix multiplications: ', numberOfMultiplications, " ", normType, " norm: ", normVal
            end if
            prevPd = endPd
            numberOfMultiplications = numberOfMultiplications + multiplicationSteps
        end do
        print*, 'sparseMat multiplications: ', numberOfMultiplications
    end subroutine run_until_converged_method

end module dense_solver_module

module sparse_solver_module
    use solver_module
    use fsparse
    use markovSparse, only: norm, convergence
    use sparseMatrix
    use iso_fortran_env, only: dp=>real64
    implicit none
    public :: sparse_solver

    ! Define a subclass for sparse sparseMat solver
    type, extends(solver) :: sparse_solver
        type(COO_dp) :: sparseMatrix
    contains
        procedure :: solve => sparse_solve_method
        procedure :: init => sparse_init_solver
        procedure :: sparseGenerateConnections => sparse_generate_connections_method
        procedure :: runUntilConverged => run_until_converged_method
        procedure, public :: numberOfGridPoints => sparse_number_of_grid_points_method
        procedure :: timeStep => time_step_sparse_method
    end type sparse_solver

contains
    subroutine sparse_init_solver(self, input_sparseMat, start_idxs, fusion_idxs, fission_idxs, dim_size)
        class (sparse_solver), intent(inout) :: self
        type(COO_dp), intent(in) :: input_sparseMat
        integer, intent(in) :: start_idxs(:), fusion_idxs(:), fission_idxs(:), dim_size(:)

        self%sparseMatrix = input_sparseMat
        self%startIdxs = start_idxs
        self%fusionIdxs = fusion_idxs
        self%fissionIdxs = fission_idxs
        self%dimSize = dim_size
    end subroutine sparse_init_solver

    function time_step_sparse_method(self, matrixSparse, probabilityCoeffs, steps) result (timeStepSparse)
        integer, intent(in) :: steps
        class(sparse_solver) :: self
        real(dp), dimension(:), intent(in) :: probabilityCoeffs
        real(dp), dimension(SIZE(probabilityCoeffs)) :: timeStepSparse, temp
        type(COO_dp) :: matrixSparse
        integer :: i
        timeStepSparse = probabilityCoeffs
        do i = 1, steps
            temp = 0
            call matvec(matrixSparse, timeStepSparse, temp)
            timeStepSparse = temp ! new timeStepSparse = matrixSparse * timeStepSparse
        end do
    end function time_step_sparse_method

    function sparse_number_of_grid_points_method(self) result (nGridPoints)
        class (sparse_solver) :: self
        integer :: nGridPoints
        nGridPoints = self%sparseMatrix%nrows
    end function sparse_number_of_grid_points_method

    function sparse_generate_connections_method(self, start_idx) result(sparseMat) !creates copy of base sparsematrix and links fusion/fission to start_idx chosen
        class (sparse_solver), intent(inout) :: self
        integer, intent(in) :: start_idx
        integer :: i, j, fusionIdx, neighbourIdx, fissionIdx
        type(COO_dp) :: sparseMat
        real(dp) :: fusionFrac = 0.25_dp !Prob to fusion
        real(dp) :: fissionFrac = 0.25_dp !prob to fission

        call sparseMat%malloc(self%sparseMatrix%nrows, self%sparseMatrix%ncols, self%sparseMatrix%nnz)
        sparseMat%index = self%sparseMatrix%index
        sparseMat%data = self%sparseMatrix%data
        do i = 1,size(self%fusionIdxs)
            fusionIdx = self%fusionIdxs(i)
            do j = 1, sparseMat%nnz
                if(sparseMat%index(2, j) == fusionIdx) then
                    neighbourIdx = sparseMat%index(1,j)
                    call sparseMat%set(sparseMat%data(j)*(1-fusionFrac), neighbourIdx , fusionIdx) !Reduce probability for all existing neighbours
                end if
            end do
            call linkStatesIdx(sparseMat, fusionIdx, start_idx, fusionFrac) !Add link between fusion coord to starting coordinate.
        end do

        do i = 1,size(self%fissionIdxs)
            fissionIdx = self%fissionIdxs(i)
            do j = 1, sparseMat%nnz
                if(sparseMat%index(2, j) == fissionIdx) then
                    neighbourIdx = sparseMat%index(1,j)
                    call sparseMat%set(sparseMat%data(j)*(1-fissionFrac), neighbourIdx , fissionIdx) !Reduce probability for all existing neighbours
                end if
            end do
            call linkStatesIdx(sparseMat, fissionIdx, start_idx, fissionFrac)
        end do
    end function sparse_generate_connections_method

    ! Sparse sparseMat solver method
    subroutine sparse_solve_method(self)
        class(sparse_solver), intent(inout) :: self
        type(COO_dp) :: sparseMat
        integer :: i, startIdx, sparseMatMultiplications
        real :: T1, T2, elapsedTime, fusionFrac, fissionFrac
        real(dp), allocatable :: startPd(:), endPd(:)
        allocate(startPd(self%numberOfGridPoints()))
        allocate(endPd(self%numberOfGridPoints()))
        startPd = 0.0_dp
        endPd = 0.0_dp

        do i = 1, SIZE(self%startIdxs)
            call cpu_time(T1)
            startIdx = self%startIdxs(i)
            sparseMat = self%sparseGenerateConnections(startIdx) !Connects fission and fusion points to starting index
            startPd = self%startingGuess() !starting guess
            call self%runUntilConverged(sparseMat, startPd, endPd, sparseMatMultiplications)
            call cpu_time(T2)

            fusionFrac = self%fusionFrac(endPd)
            fissionFrac = self%fissionFrac(endPd)
            elapsedTime = T2-T1
            call self%addResult(endPd, elapsedTime, sparseMatMultiplications, self%idxToCoords(startIdx))
        end do
    end subroutine sparse_solve_method

    

    subroutine run_until_converged_method(self, sparseMat, startPd, endPd, numberOfMultiplications)
        class(sparse_solver), intent(inout) :: self
        type(COO_dp) :: sparseMat
        real(dp), dimension(:), intent(in) :: startPd
        real(dp), dimension(:), intent(inout) :: endPd
        real(dp), dimension(SIZE(startPd)) :: prevPd, diff
        integer, intent(inout) :: numberOfMultiplications
        real(dp) :: tol, normVal
        integer :: multiplicationSteps
        logical :: converged
        character (len = 2):: normType = "L2"

        numberOfMultiplications = 0
        tol = 1.0/(1e2*SIZE(startPd)) !This tolerance is one part in one hundered assuming even spread of probability.
        multiplicationSteps = 20
        tol = tol * multiplicationSteps / 10 !Dynamically change tolerance based on number of sparseMat multiplications in a row
        print*, tol
        
        converged = .FALSE.
        prevPd = startPd
        do while (.not. converged)
            endPd = self%timeStep(sparseMat, prevPd, multiplicationSteps) !Dont check for convergence after every step, rather take a few time steps at a time.
            converged = convergence(endPd, prevPd, tol, normType)
            
            if(mod(numberOfMultiplications, 500) == 0) then

                diff = endPd - prevPd
                normVal = norm(diff, normType)
                print *, 'Matrix multiplications: ', numberOfMultiplications, " ", normType, " norm: ", normVal
            end if
            prevPd = endPd
            numberOfMultiplications = numberOfMultiplications + multiplicationSteps
        end do
        print*, 'sparseMat multiplications: ', numberOfMultiplications
    end subroutine run_until_converged_method
end module sparse_solver_module

module sparse_solver_arnoldi_module
    use sparse_solver_module
    use fsparse
    use markovSparse
    use iso_fortran_env, only: dp=> real64

    type, extends(sparse_solver) :: sparse_arnoldi_solver
    contains
        procedure :: runUntilConverged => run_until_converged_method_arnoldi
        
    end type sparse_arnoldi_solver

    contains
        subroutine run_until_converged_method_arnoldi(self, sparseMat, startPd, endPd, numberOfMultiplications)
            class(sparse_arnoldi_solver), intent(inout) :: self
            type(COO_dp) :: sparseMat
            real(dp), dimension(:), intent(in) :: startPd
            real(dp), dimension(:), intent(inout) :: endPd
            real(dp), dimension(SIZE(startPd)) :: prevPd, diff
            integer, intent(inout) :: numberOfMultiplications
            real(dp) :: tol, normVal
            integer :: multiplicationSteps
            logical :: converged
            character (len = 2):: normType = "L2"


        end subroutine run_until_converged_method_arnoldi
end module sparse_solver_arnoldi_module

module sparse_solver_linear_interpolator_module
    use sparse_solver_module
    use fsparse
    use markovSparse
    use iso_fortran_env, only: dp=>real64

    type, extends(sparse_solver) :: sparse_linearint_solver
    contains
        procedure :: startingGuess => starting_guess_linear_interpolation_method
    end type sparse_linearint_solver

contains
    function starting_guess_linear_interpolation_method(self) result(guess)
        class(sparse_linearint_solver), intent(inout) :: self
        real(dp), allocatable :: guess(:), newerSol(:), olderSol(:)
        integer :: storedResults
        real, allocatable :: newerSolCoord(:), olderSolCoord(:), startCoord(:)
        allocate(guess(self%numberOfGridPoints()))
        allocate(olderSol(self%numberOfGridPoints()))
        allocate(newerSol(self%numberOfGridPoints()))
        allocate(newerSolCoord(SIZE(self%dimSize)))
        allocate(olderSolCoord(SIZE(self%dimSize)))
        allocate(startCoord(SIZE(self%dimSize)))
        guess = 0
        storedResults = self%result%numResults
        
        if(storedResults == 0) then !Guess even probability distribution
            guess = 1.0_dp/self%numberOfGridPoints()
        else
            guess = self%result%getProbabilityDensity(storedResults) !guess last solution as solution
             !even distribution as first guess
        end if
        !TODO:do linear interpolation


    end function starting_guess_linear_interpolation_method

end module sparse_solver_linear_interpolator_module