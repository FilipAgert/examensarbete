module solver_module
    use iso_fortran_env, only: dp=>real64
    use result_module
    
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
    contains
        procedure(solve_interface), deferred, public :: solve  ! Define abstract method 'solve'
        procedure, public :: printResult => print_result_method
        procedure, public :: generateConnections => generate_connections_method
        procedure :: fusionFrac => fusion_frac_method
        procedure :: fissionFrac => fission_frac_method
        procedure, public :: addResult => add_result_method
        procedure, public :: coordinateNumber => get_coordinate_number_method
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
    function get_coordinate_number_method(self)result(coordinateN)
        class(solver), intent(inout) :: self
        integer :: coordinateN
        coordinateN = SIZE(self%matrix, 1)
    end function get_coordinate_number_method
    ! Constructor to initialize the sparseMat in the solver base class
    subroutine init_solver(self, input_matrix, start_idxs, fusion_idxs, fission_idxs)
        class(solver), intent(inout) :: self
        real(dp), intent(in) :: input_matrix(:,:)
        integer, intent(in) :: start_idxs(:), fusion_idxs(:), fission_idxs(:)
        
        self%matrix = input_matrix  ! Store the sparseMat
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
    use markovSparse
    use iso_fortran_env, only: dp=>real64
    implicit none
    public :: sparse_solver

    ! Define a subclass for sparse sparseMat solver
    type, extends(solver) :: sparse_solver
    contains
        procedure :: solve => sparse_solve_method
        procedure :: init => sparse_init_solver
        procedure :: sparseGenerateConnections => sparse_generate_connections_method
        procedure :: runUntilConverged => run_until_converged_method
        procedure :: startingGuess => starting_guess_method
    end type sparse_solver

contains

    
    subroutine sparse_init_solver(self, input_sparseMat, start_idxs, fusion_idxs, fission_idxs)
        class (sparse_solver), intent(inout) :: self
        real(dp), intent(in) :: input_sparseMat(:,:)
        integer, intent(in) :: start_idxs(:), fusion_idxs(:), fission_idxs(:)

        
        call init_solver(self, input_sparseMat, start_idxs, fusion_idxs, fission_idxs) !Call base initializer
    end subroutine sparse_init_solver

    function sparse_generate_connections_method(self, start_idx) result (sparseMat)
        class (sparse_solver), intent(inout) :: self
        integer, intent(in) :: start_idx
        real(dp), allocatable :: denseMat(:,:)
        type(COO_dp) :: sparseMat
        
        denseMat = self%generateConnections(start_idx)
        call dense2coo(denseMat, sparseMat)
    end function sparse_generate_connections_method

    ! Sparse sparseMat solver method
    subroutine sparse_solve_method(self)
        class(sparse_solver), intent(inout) :: self
        type(COO_dp) :: sparseMat
        integer :: i, startIdx, sparseMatMultiplications
        real :: T1, T2, elapsedTime, fusionFrac, fissionFrac
        real(dp), dimension(SIZE(self%matrix,1)) :: startPd, endPd
        real :: dummy(1)
        startPd = 0.0_dp
        endPd = 0.0_dp

        dummy(1) = 1.0


        do i = 1, SIZE(self%startIdxs)
            call cpu_time(T1)
            startIdx = self%startIdxs(i)
            sparseMat = self%sparseGenerateConnections(startIdx)
            startPd = self%startingGuess() !starting guess

            call self%runUntilConverged(sparseMat, startPd, endPd, sparseMatMultiplications)
            call cpu_time(T2)

            fusionFrac = self%fusionFrac(endPd)
            fissionFrac = self%fissionFrac(endPd)
            elapsedTime = T2-T1
            call self%addResult(endPd, elapsedTime, sparseMatMultiplications, dummy)
        end do
    end subroutine sparse_solve_method

    function starting_guess_method(self) result(guess) !Guesses the initial probability distribution
        class(sparse_solver), intent(inout) :: self
        real(dp), allocatable :: guess(:), prevGuess(:)
        integer :: storedResults
        allocate(guess(self%coordinateNumber()))
        allocate(prevGuess(self%coordinateNumber())) !Allocate array to correct size
        guess = 0
        storedResults = self%result%numResults

        if(storedResults > 0) then !Use previous prob distribution as initial guess
            prevGuess = self%result%getProbabilityDensity(storedResults) !Index latest guess
            guess = prevGuess
        else
            guess = 1.0_dp/self%coordinateNumber() !even distribution as first guess
        end if
    end function starting_guess_method

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
        tol = 1.0/(1e6) !Tolerance is one part in one million.
        multiplicationSteps = 20
        !tol = tol * multiplicationSteps !Dynamically change tolerance based on number of sparseMat multiplications in a row
        
        converged = .FALSE.
        prevPd = startPd
        do while (.not. converged)
            endPd = timeStepSparse(sparseMat, prevPd, multiplicationSteps) !Dont check for convergence after every step, rather take a few time steps at a time.
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
        real(dp), allocatable :: guess(:), prevGuess(:)
        integer :: storedResults
    end function starting_guess_linear_interpolation_method

end module sparse_solver_linear_interpolator_module