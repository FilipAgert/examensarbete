module solver_module
    use iso_fortran_env, only: dp=>real64, sp=>real32
    use result_module
    use denseMatrix, only: coordFromLinearIdx
    
    implicit none
    private
    public :: solver, create_solver, init_solver

    ! Define the base solver class as a derived type
    type, abstract :: solver
        real(sp), allocatable :: matrix(:,:)  ! The sparseMat is dynamically sized. sparseMat should have no connections beyond neigbours (no fusion/fission)
        type(convergedResult) :: result
        integer, allocatable :: startIdxs(:)
        real, allocatable :: energies(:)
        integer, allocatable :: fissionIdxs(:)
        integer, allocatable :: fusionIdxs(:)
        integer, allocatable :: dimSize(:)
        real, allocatable :: grid(:,:) !List of all coordinate points
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
    subroutine init_solver(self, input_matrix, grid, start_idxs, energies, fusion_idxs, fission_idxs, dim_size)
        class(solver), intent(inout) :: self
        real(sp), intent(in) :: input_matrix(:,:)
        real, intent(in) :: energies(:), grid(:,:)
        integer, intent(in) :: start_idxs(:), fusion_idxs(:), fission_idxs(:), dim_size(:)
        
        self%matrix = input_matrix  ! Store the sparseMat
        self%startIdxs = start_idxs
        self%energies = energies
        self%fusionIdxs = fusion_idxs
        self%fissionIdxs = fission_idxs
        self%dimSize = dim_size
        self%grid = grid
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
        real(sp), allocatable :: mat(:,:)
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

    subroutine add_result_method(self, pd, time, multiplications, startCoord, energy)
        class(solver), intent(inout) :: self
        real(dp), intent(in) :: pd(:)
        real, intent(in) :: time, energy
        integer, intent(in) :: multiplications
        real, intent(in) :: startCoord(:)

        call self%result%addResult(pd, startCoord, energy, time, multiplications, self%fusionFrac(pd), self%fissionFrac(pd))
    end subroutine add_result_method

end module solver_module

module dense_solver_module
    use solver_module
    use iso_fortran_env, only: dp => real64, sp=>real32
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

    subroutine dense_init_solver(self, denseMat, grid, start_idxs, energies, fusion_idxs, fission_idxs, dimSize)
        class(dense_solver), intent(inout) :: self
        real(sp), dimension(:,:), intent(inout) :: denseMat
        integer, intent(inout) :: start_idxs(:), fusion_idxs(:), fission_idxs(:), dimSize(:)
        real, intent(inout) :: grid(:,:), energies(:)
        call create_solver(self, denseMat, grid, start_idxs, energies, fusion_idxs, fission_idxs, dimSize)
    end subroutine dense_init_solver
    
    subroutine dense_solve_method(self)
        class(dense_solver), intent(inout) :: self
        real(sp), dimension(:,:), allocatable :: denseMat
        integer :: i, startIdx, matMultiplications
        real :: T1, T2, elapsedTime, fusionFrac, fissionFrac
        real(dp), allocatable :: startPd(:), endPd(:)
        allocate(startPd(self%numberOfGridPoints()))
        allocate(endPd(self%numberOfGridPoints()))
        allocate(denseMat(self%numberOfGridPoints(), self%numberOfGridPoints()))
        startPd = 0.0_dp
        endPd = 0.0_dp

        do i = 1, SIZE(self%startIdxs)
            print* , ' '
            print*, "Running calculation for energy: ", self%energies(i)
            call cpu_time(T1)
            startIdx = self%startIdxs(i)
            denseMat = self%generateConnections(startIdx) !Connects fission and fusion points to starting index
            startPd = self%startingGuess() !starting guess
            call self%runUntilConverged(denseMat, startPd, endPd, matMultiplications)
            call cpu_time(T2)

            fusionFrac = self%fusionFrac(endPd)
            fissionFrac = self%fissionFrac(endPd)
            elapsedTime = T2-T1
            print*, ' '
            print*, "Time taken: ", elapsedTime, "s"
            print*, "Number of matrix multiplications: ", matMultiplications
            call self%addResult(endPd, elapsedTime, matMultiplications, self%idxToCoords(startIdx), self%energies(i))
        end do
    end subroutine dense_solve_method

    function time_step_method(self, matrix, probabilityCoeffs, steps)result (timeStep) !This computes repeated matrix multiplication (A*(A*...*(A*V)
        class(dense_solver), intent(inout) :: self
        integer, intent(in) :: steps                     !In order to only have to store one matrix
        real(dp), dimension(:), intent(in) :: probabilityCoeffs
        real(dp), dimension(SIZE(probabilityCoeffs)) :: timeStep
        real, dimension(SIZE(probabilityCoeffs), SIZE(probabilityCoeffs)), intent(in):: matrix 
        integer :: i
        timeStep = probabilityCoeffs
        do i = 1, steps
            timeStep = MATMUL(matrix, timeStep)
        end do
        
    end function time_step_method
    subroutine run_until_converged_method(self, denseMat, startPd, endPd, numberOfMultiplications)
        class(dense_solver), intent(inout) :: self
        real(sp), dimension(:,:), intent(in) :: denseMat
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
        multiplicationSteps = 10
        tol = tol * multiplicationSteps / 10 !Dynamically change tolerance based on number of sparseMat multiplications in a row
        print*, "Tolerance: ", tol
        
        converged = .FALSE.
        prevPd = startPd
        do while (.not. converged)
            endPd = self%timeStep(denseMat, prevPd, multiplicationSteps) !Dont check for convergence after every step, rather take a few time steps at a time.
            converged = convergence(endPd, prevPd, tol, normType)
            
            if(mod(numberOfMultiplications, 5000) == 0) then

                diff = endPd - prevPd
                normVal = norm(diff, normType)
                print *, 'Matrix multiplications: ', numberOfMultiplications, " ", normType, " norm: ", normVal
            end if
            prevPd = endPd
            numberOfMultiplications = numberOfMultiplications + multiplicationSteps
        end do
    end subroutine run_until_converged_method

end module dense_solver_module

module sparse_solver_module
    use solver_module
    use fsparse
    use markovSparse, only: norm, convergence
    use sparseMatrix
    use iso_fortran_env, only: dp=>real64, sp=>real32
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
        procedure :: timeStep => time_step_dparse_method
    end type sparse_solver

    contains
    subroutine sparse_init_solver(self, input_dparseMat, grid, start_idxs, energies, fusion_idxs, fission_idxs, dim_size)
        class (sparse_solver), intent(inout) :: self
        type(COO_dp), intent(in) :: input_dparseMat
        integer, intent(in) :: start_idxs(:), fusion_idxs(:), fission_idxs(:), dim_size(:)
        real, intent(in) :: grid(:,:), energies(:)

        self%sparseMatrix = input_dparseMat
        self%startIdxs = start_idxs
        self%fusionIdxs = fusion_idxs
        self%fissionIdxs = fission_idxs
        self%dimSize = dim_size
        self%energies = energies
        self%grid = grid
    end subroutine sparse_init_solver

    function time_step_dparse_method(self, matrixSparse, probabilityCoeffs, steps) result (timeStepSparse)
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
            timeStepSparse = temp/sum(temp) ! new timeStepSparse = matrixSparse * timeStepSparse
        end do
    end function time_step_dparse_method

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
        real(sp) :: fusionFrac = 0.25_sp !Prob to fusion
        real(sp) :: fissionFrac = 0.25_sp !prob to fission

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
            call linkStatesIdx(sparseMat, fusionIdx, start_idx, real(fusionFrac,8)) !Add link between fusion coord to starting coordinate.
        end do

        do i = 1,size(self%fissionIdxs)
            fissionIdx = self%fissionIdxs(i)
            do j = 1, sparseMat%nnz
                if(sparseMat%index(2, j) == fissionIdx) then
                    neighbourIdx = sparseMat%index(1,j)
                    call sparseMat%set(sparseMat%data(j)*(1-fissionFrac), neighbourIdx , fissionIdx) !Reduce probability for all existing neighbours
                end if
            end do
            call linkStatesIdx(sparseMat, fissionIdx, start_idx, real(fissionFrac,8))
        end do
    end function sparse_generate_connections_method

    ! Sparse sparseMat solver method
    subroutine sparse_solve_method(self)
        class(sparse_solver), intent(inout) :: self
        type(COO_dp) :: sparseMat
        integer :: i, startIdx, sparseMatMultiplications
        real :: T1, T2, elapsedTime, fusionFrac, fissionFrac, T3
        real(dp), allocatable :: startPd(:), endPd(:)
        allocate(startPd(self%numberOfGridPoints()))
        allocate(endPd(self%numberOfGridPoints()))
        startPd = 0.0_dp
        endPd = 0.0_dp

        do i = 1, SIZE(self%startIdxs)
            print* , ' '
            print*, "Running calculation for energy: ", self%energies(i)
            call cpu_time(T1)
            startIdx = self%startIdxs(i)
            sparseMat = self%sparseGenerateConnections(startIdx) !Connects fission and fusion points to starting index
            call cpu_time(T2)
            startPd = self%startingGuess() !starting guess
            call self%runUntilConverged(sparseMat, startPd, endPd, sparseMatMultiplications)
            call cpu_time(T3)

            fusionFrac = self%fusionFrac(endPd)
            fissionFrac = self%fissionFrac(endPd)
            elapsedTime = T3-T1
            print*, "Time taken: ", elapsedTime, "s", " Of which ", 100*(T2-T1)/(T3-T1), "percent spent on generating matrix."
            print*, "Number of matrix multiplications: ", sparseMatMultiplications
            call self%addResult(endPd, elapsedTime, sparseMatMultiplications, self%idxToCoords(startIdx), self%energies(i))
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
        tol = 1.0/(SIZE(startPd)*1e2) !This tolerance is one part in one hundered assuming even spread of probability.
        multiplicationSteps = 20
        tol = tol * multiplicationSteps / 10 !Dynamically change tolerance based on number of sparseMat multiplications in a row
        print*, "Tolerance: ",tol
        
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
    end subroutine run_until_converged_method
end module sparse_solver_module



module sparse_solver_linear_interpolator_module
    use sparse_solver_module
    use fsparse
    use markovSparse
    use denseMatrix
    use linearInterpolation
    use iso_fortran_env, only: dp=>real64, sp=>real32

    type, extends(sparse_solver) :: sparse_linearint_solver
    contains
        procedure :: startingGuess => starting_guess_linear_interpolation_method
    end type sparse_linearint_solver

    contains
        function starting_guess_linear_interpolation_method(self) result(guess)
            class(sparse_linearint_solver), intent(inout) :: self
            real(dp), allocatable :: guess(:)
            integer :: storedResults, i, E0idx, E1idx, tempidx
            real :: energy, E0, E1, diff, potentialClosest, tempE, E1DIFF, E0DIFF
            real(dp), allocatable :: pdf0(:), pdf1(:) 
            allocate(guess(self%numberOfGridPoints()))
            guess = 0
            storedResults = self%result%numResults
            !guess = 1.0/self%numberOfGridPoints()
            !return
            if(storedResults == 0) then !Guess evenly distributed probability distribution
                guess = real(1.0_dp/real(self%numberOfGridPoints(),8),4) !Better guess than everything
                print*, "GUESS = equal probability in all states"
                return
            elseif(storedResults == 1) then
                guess = self%result%getProbabilityDensity(storedResults)
                print*, "GUESS = RE-USE 1st solution"
            else
                !guess = self%result%getProbabilityDensity(storedResults) !LINEAR INTERPOLATION SEEMS TO BE BAD.
                !return INVESTIGATE
                !THIS BLOCK OF CODE GETS INDICES AND ENERGIES OF TWO CLOSEST ENERGIES
                !E1 is closer in energy to ENERGY 
                !E0 is second closest in energy
                E0 = huge(E0)
                E1 = huge(E1) - 1

                E0idx = storedResults - 1
                E1idx = storedResults
                energy = self%energies(storedResults + 1) !get energy to guess for.

                do i = 1, storedResults !find two closest energies as they will produce best fit.
                    potentialClosest = self%result%energies(i)
                    diff = abs(energy - potentialClosest)
                    E0DIFF = abs(E0-energy)
                    E1DIFF = abs(E1-energy)
                    if(E0DIFF > E1DIFF) then 
                        if(diff < E0DIFF) then
                            E0 = potentialClosest
                            E0idx = i
                        endif
                    else
                        if(diff < E1DIFF) then
                            E1 = potentialClosest
                            E1idx = i
                        endif
                    endif
                end do

                if(abs(E0-energy) < abs(E1-energy)) then
                    tempE = E0
                    E0 = E1
                    E1 = tempE
                    tempidx = E0idx
                    E0idx = E1idx
                    E1idx = tempidx
                endif
                
            
                allocate(pdf0(self%numberOfGridPoints()), pdf1(self%numberOfGridPoints()))

                pdf0 = real(self%result%getProbabilityDensity(E0idx),4)
                pdf1 = real(self%result%getProbabilityDensity(E1idx),4)
                guess = guessPdfFast(self%grid, energy, E0, E1, pdf0, pdf1)!
                print*, "Guess = LINCOMB of IDX:", E0idx, " and ", E1idx
                
            endif
        end function starting_guess_linear_interpolation_method
end module sparse_solver_linear_interpolator_module

module sparse_solver_arnoldi_module
    use sparse_solver_linear_interpolator_module
    use fsparse
    use linearInterpolation
    use iso_fortran_env, only: dp=> real64, sp=>real32

    type, extends(sparse_linearint_solver) :: sparse_arnoldi_solver
    contains
        procedure :: runUntilConverged => run_until_converged_method_arnoldi
        
    end type sparse_arnoldi_solver

    contains
        subroutine run_until_converged_method_arnoldi(self, sparseMat, startPd, endPd, numberOfMultiplications)
            implicit none
            class(sparse_arnoldi_solver), intent(inout) :: self
            type(COO_dp) :: sparseMat
            real(dp), dimension(:), intent(in) :: startPd
            real(dp), dimension(:), intent(inout) :: endPd
            real(dp), dimension(SIZE(startPd)) :: temp
            integer, intent(inout) :: numberOfMultiplications
            

            external :: dnaupd
            external :: dneupd
            integer :: IDO
            character(len = 1) :: BMAT = 'I'
            integer :: N
            character(len = 2) ::  WHICH = 'LR'  !since multiple eigenvalues can have abs(lambda)=1, we need to choose only largest real part
            integer :: NEV!number of eigenvalues comptued
            real(dp) :: TOL 
            real(dp),allocatable :: RESID(:)
            integer :: NCV !number of vectors in calculation. computation scales as N * NCV²
            real(dp), allocatable :: V(:,:)
            integer :: IPARAM(11)
            integer :: ishfts
            integer :: maxitr !max iterations?
            integer :: mode = 1 !1 = no sigma shift. 3 = sigma shift of one
            integer :: LDV
            Real(dp), allocatable :: workd(:), workl(:), workev(:)
            integer :: LWORKL
            integer :: INFO
            integer :: INFO2 !0 = randomized initial vector, 1 = resid initial vector (use starting guess)
            logical :: converged
            integer :: IPNTR(14)
            LOGICAL :: RVEC
            character(len = 1) :: HOWMNY = 'A'
            logical, allocatable :: SELECT(:)
            real(dp), dimension(:), allocatable :: DR, DI
            real(dp), dimension(:,:), allocatable :: Z
            real(dp) :: SIGMAR, SIGMAI, TIME1, TIME2
            integer :: LDZ
            N = self%numberOfGridPoints()
            !!USER SETTINGS.
            maxitr = 50000
            NEV = 1 !number of eigenvalues calculated. test to change this
            TOL = 1./(self%numberOfGridPoints())!Tolerance is: how close to true eigenvalue can the calculated one be?
            !Tolerance means that you expect no eigenvalues to be closer to eachother than tolerance.
            !seems to only converge to true value if low enough. Should maybe be of order 1/number of grid points For 501x501
            NCV = MAX(2*NEV + 1,int((N)**(1.0/5.0))) !set at least to 2*NEV. Lower number = more matrix * vector operations.
            !!!!!!
            SIGMAR = 0 
            SIGMAI = 0
            INFO = 1 !info for dnaupd          1 = user specified vector
            INFO2 = 0 !info for dneupd         1 = user specified vector
            IDO = 0
            ishfts = 1
            LDZ = N
            allocate(Z(N, NEV + 1))
            allocate(DR(NEV + 1), DI(NEV + 1))
            allocate(SELECT(NCV))

            LWORKL = 3*NCV**2 + 6 *NCV
            allocate(workl(LWORKL))
            allocate(workd(3*N), workev(3*NCV))
            allocate(RESID(N))
            allocate(V(N,NCV))
            IPARAM = 0
            IPARAM(1) = ishfts
            IPARAM(3) = maxitr !Number of arnoldi = .FALSE.update iterations
            IPARAM (4) = 1 !needs to be 1
            IPARAM(7) = mode !specifies eigenvalue problem A*x = lambda*x

            
            LDV = N
            RESID = startPd!starting guess
            numberOfMultiplications = 0
            converged = .FALSE.
            do while(.not. converged)
                call dnaupd(IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV,IPARAM,IPNTR, workd, workl, lworkl,info)
                
                if(IDO .eq. -1 .or. IDO .eq. 1) then
                    temp = 0

                    
                    call matvec(sparseMat, workd(ipntr(1):ipntr(1) + N - 1),temp)
                    workd(ipntr(2) : ipntr(2) + N - 1) = temp
                    numberOfMultiplications = numberOfMultiplications + 1
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
                            LDV, IPARAM, IPNTR, WORKD, WORKL, LWORKL, INFO2)


                if(INFO2 .ne. 0) then
                    print *, ' '
                    print *, ' Error with _neupd, info = ', INFO2
                    print *, ' Check the documentation of _neupd. '
                    print *, ' '
                else
                    print *, 'Found eigenvector'
                    print *, 'Eigenvalue(s):', DR
                    endPd = Z(:,1)/sum(Z(:,1)) !Ensure correct phase.
                endif
            endif
            continue
        end subroutine run_until_converged_method_arnoldi
end module sparse_solver_arnoldi_module

module sparse_linear_system_solver_module
    use sparse_solver_linear_interpolator_module
    use fsparse
    use denseMatrix
    use mgmres
    use iso_fortran_env, only: dp=> real64, sp=>real32
    
    type, extends(sparse_linearint_solver) :: sparse_linear_system_solver
    contains
        procedure :: runUntilConverged => run_until_converged_lin_system
        
    end type sparse_linear_system_solver

    contains
        subroutine run_until_converged_lin_system(self, sparseMat, startPd, endPd, numberOfMultiplications)
            implicit none
            class(sparse_linear_system_solver), intent(inout) :: self
            type(COO_dp) :: sparseMat
            real(dp), dimension(:), intent(in) :: startPd
            real(dp), dimension(:), intent(inout) :: endPd
            integer, intent(inout) :: numberOfMultiplications
            type(CSR_dp) :: sparseMatCSR, sparseMatCSR2
            integer :: i, N, j
            integer, allocatable :: IA(:), JA(:)
            real(dp), allocatable :: A(:), X1(:),  RHS(:), SOL(:), a12(:), a21(:)
            integer :: ITR_MAX, MR, row, col, numchanges, NZ_NUM
            real(dp) :: TOL_ABS, TOL_REL, val, alpha
            logical :: verbose
            numberOfMultiplications = 0
            N = self%numberOfGridPoints()
            
            
            numChanges = 0
            
            do i = 1,sparseMat%nnz !iterate over all nonzero elements, if diagonal element, subtract one from it.
                row = sparseMat%index(1,i)
                col = sparseMat%index(2,i)
                if(row == col) then
                    val = sparseMat%data(i)
                    sparseMat%data(i) = val - 1
                    numChanges = numChanges + 1
                end if
            end do
            if(numChanges /= N) then !
                print*, "WARNING: MATRIX NOT SET UP CORRECTLY, DID NOT SUBTRACT -1 FROM ALL DIAGONAL ELEMENTS"
                print*, "N: ", N
                print*, "Diagonal elements found: ", numchanges
            end if
            
            !Then convert to CSR format
            call COO2ordered(sparseMat,.TRUE.)
            call coo2csr(sparseMat, sparseMatCSR)
            call remove_last_row_and_col_CSR(sparseMatCSR, sparseMatCSR2,a21, a12, alpha)
            
            !Call mgmres code
            NZ_NUM = sparseMatCSR2%nnz
            allocate(IA(N), JA(NZ_NUM), A(NZ_NUM), X1(N-1), SOL(N-1),RHS(N-1))
            RHS = - a12
            IA = sparseMatCSR2%rowptr
            JA = sparseMatCSR2%col
            A = sparseMatCSR2%data
            X1 = startPd(1:N-1)
            !MODIFY THESE
            ITR_MAX = 10
            MR = 400
            TOL_ABS = 1.0/N
            TOL_REL = 1e-3

            !MODIFY THESE
            verbose = .TRUE.
            ! print*, "IA;", size(IA)
            ! print*, "JA:", SIZE(JA)
            ! print*, "A: ", SIZE(A)
            ! print*, "NNZ: ", NZ_NUM
            ! print*, "N: ", N-1
            ! print*, "x1:", size(X1)
            ! print*, "RHS: ", size(RHS)
            ! print*, "maxinptr (should be nnz + 1) ", maxval(IA)
            ! print*, "maxval col", maxval(JA)
            
            call pmgmres_ilu_cr(N-1, NZ_NUM, IA, JA, A, X1, RHS, ITR_MAX, MR, TOL_ABS, TOL_REL,verbose)
            
            endPd(1:N-1) = X1
            ! print*, "SUM X1: ", norm(X1)
            deallocate(IA, JA, A, X1, SOL,RHS, sparseMatCSR%col, sparseMatCSR%data, sparsematcsr%rowptr)
            endPd(N) = 1
            endPd = endPd/sum(endPd)

        end subroutine run_until_converged_lin_system

        subroutine remove_last_row_and_col_CSR(CSRIN, CSROUT, rowOut, colOut, alpha)
            !!Removes last row and column from matrix. Returns row and col as vector except for last element. The last element is returned as alpha 
            !!This assumes last element exists
            implicit none
            type(CSR_dp), intent(inout) :: CSRIN, CSROUT
            real(dp), intent(out), allocatable :: rowOut(:), colOut(:)
            real(dp) :: alpha
            integer :: NNZchange, i, N, nnz, rowPtr, rowPtrNext, j, col, nonZeroIndex, colNext(CSRIN%nnz)
            integer :: newRowPtr(CSRIN%nrows)
            real(dp) :: dataNext(CSRIN%nnz)
            allocate(rowOut(CSRIN%ncols - 1), colOut(CSRIN%nrows -1))


            N = CSRIN%nrows !Assume nrows = ncols
            NNZchange = 0
            
            rowOut = 0
            colOut = 0
            !Extracts last row
            alpha = 0
            do i = CSRIN%rowptr(N), CSRIN%rowptr(N+1)-1   
                if(CSRIN%col(i) /= N) then
                    rowOut(CSRIN%col(i)) = CSRIN%data(i)
                    NNZchange = NNZchange + 1
                else
                    alpha = CSRIN%data(i)
                    
                    NNZchange = NNZchange + 1
                endif
            end do
            nnz = CSRIN%nnz
            newRowPtr = CSRIN%rowptr(1:N)
            nonZeroIndex = 0
            do i = 1, N - 1 !Iterate over all rows except last row.
                !Check if row has last column. If so, add to colout.
                rowPtr = CSRIN%rowptr(i)! This is start of row indices
                rowPtrNext = CSRIN%rowptr(i+1)
                do j = rowPtr, rowPtrNext -1
                    col = CSRIN%col(j)
                    if(col == CSRIN%ncols) then !If last column, copy data to new vector.
                        colOut(i) = CSRIN%data(j) !Copy data over.
                        newRowPtr(i+1:N) = newRowPtr(i+1:N) - 1 !Shift rowptr
                        NNZchange = NNZchange + 1 !Count number of elements copied
                    else !Else, copy to new matrix.
                        nonZeroIndex = nonZeroIndex + 1
                        dataNext(nonZeroIndex) = CSRIN%data(j)
                        colNext(nonZeroIndex) = CSRIN%col(j)
                    endif
                end do
            end do
            call CSROUT%malloc(CSRIN%nrows-1, CSRIN%ncols - 1, nonZeroIndex)
            CSROUT%data = dataNext(1:nonZeroIndex)
            CSROUT%col = colNext(1:nonZeroIndex)
            CSROUT%rowptr = newRowPtr
            
        end subroutine remove_last_row_and_col_CSR
    
end module sparse_linear_system_solver_module

module sparse_linear_system_solver_bicg_module
    !https://people.math.sc.edu/Burkardt/f_src/dlap/dlap.html
    !https://mathoverflow.net/questions/410566/matrix-free-linear-solve-for-nullspace
    use sparse_solver_linear_interpolator_module
    use fsparse
    use denseMatrix
    use dlap
    use iso_fortran_env, only: dp=> real64, sp=>real32
    
    type, extends(sparse_linearint_solver) :: sparse_linear_system_solver_bicg
    contains
        procedure :: runUntilConverged => run_until_converged_lin_system_bicg
        
    end type sparse_linear_system_solver_bicg

    contains
        subroutine run_until_converged_lin_system_bicg(self, sparseMat, startPd, endPd, numberOfMultiplications)
            implicit none
            class(sparse_linear_system_solver_bicg), intent(inout) :: self
            type(COO_dp) :: sparseMat, sparseMatFinal
            real(dp), dimension(:), intent(in) :: startPd
            real(dp), dimension(:), intent(inout) :: endPd
            integer, intent(inout) :: numberOfMultiplications
            integer :: i, N, j, newN, ISYM, ITOL
            integer, allocatable :: IA(:), JA(:), IWORK(:)
            real(dp), allocatable :: A(:), X1(:),  X(:), SOL(:), a12(:), a21(:), B(:), RWORK(:)
            integer :: row, col, numchanges, NELT, ITMAX, ITER, IERR, IUNIT, LENW, LENIW
            real(dp) :: val, alpha, TOL, ERR
            
            
            numberOfMultiplications = 0
            N = self%numberOfGridPoints()
            
            
            numChanges = 0
    
            do i = 1,sparseMat%nnz !iterate over all nonzero elements, if diagonal element, subtract one from it.
                row = sparseMat%index(1,i)
                col = sparseMat%index(2,i)
                if(row == col) then
                    val = sparseMat%data(i)
                    sparseMat%data(i) = val - 1
                    numChanges = numChanges + 1
                end if
            end do
            if(numChanges /= N) then !
                print*, "WARNING: MATRIX NOT SET UP CORRECTLY, DID NOT SUBTRACT -1 FROM ALL DIAGONAL ELEMENTS"
                print*, "N: ", N
                print*, "Diagonal elements found: ", numchanges
            end if
            
            !call COO2ordered(sparseMat,.TRUE.)
            
            call remove_last_row_and_col_COO(sparseMat, sparseMatFinal, a21, a12, alpha)
            
            ! *Arguments:
! N      :IN       Integer.
!         Order of the Matrix.
! B      :IN       Double Precision B(N).
!         Right-hand side vector.
! X      :INOUT    Double Precision X(N).
!         On input X is your initial guess for solution vector.
!         On output X is the final approximate solution.
! NELT   :IN       Integer.
!         Number of Non-Zeros stored in A.
! IA     :INOUT    Integer IA(NELT).
! JA     :INOUT    Integer JA(NELT).
! A      :INOUT    Double Precision A(NELT).
!         These arrays should hold the matrix A in either the SLAP
!         Triad format or the SLAP Column format.  See "Description",
!         below.  If the SLAP Triad format is chosen it is changed
!         internally to the SLAP Column format.
! ISYM   :IN       Integer.
!         Flag to indicate symmetric storage format.
!         If ISYM=0, all nonzero entries of the matrix are stored.
!         If ISYM=1, the matrix is symmetric, and only the upper
!         or lower triangle of the matrix is stored.
! ITOL   :IN       Integer.
!         Flag to indicate type of convergence criterion.
!         If ITOL=1, iteration stops when the 2-norm of the residual
!         divided by the 2-norm of the right-hand side is less than TOL.
!         If ITOL=2, iteration stops when the 2-norm of M-inv times the
!         residual divided by the 2-norm of M-inv times the right hand
!         side is less than TOL, where M-inv is the inverse of the
!         diagonal of A.
!         ITOL=11 is often useful for checking and comparing different
!         routines.  For this case, the user must supply the "exact"
!         solution or a very accurate approximation (one with an error
!         much less than TOL) through a common block,
!         COMMON /SOLBLK/ SOLN( )
!         if ITOL=11, iteration stops when the 2-norm of the difference
!         between the iterative approximation and the user-supplied
!         solution divided by the 2-norm of the user-supplied solution
!         is less than TOL.
! TOL    :IN       Double Precision.
!         Convergence criterion, as described above.
! ITMAX  :IN       Integer.
!         Maximum number of iterations.
! ITER   :OUT      Integer.
!         Number of iterations required to reach convergence, or
!         ITMAX+1 if convergence criterion could not be achieved in
!         ITMAX iterations.
! ERR    :OUT      Double Precision.
!         Error estimate of error in final approximate solution, as
!         defined by ITOL.
! IERR   :OUT      Integer.
!         Return error flag.
!           IERR = 0 => All went well.
!           IERR = 1 => Insufficient storage allocated
!                       for WORK or IWORK.
!           IERR = 2 => Method failed to converge in
!                       ITMAX steps.
!           IERR = 3 => Error in user input.  Check input
!                       value of N, ITOL.
!           IERR = 4 => User error tolerance set too tight.
!                       Reset to 500.0*D1MACH(3).  Iteration proceeded.
!           IERR = 5 => Preconditioning matrix, M,  is not
!                       Positive Definite.  $(r,z) < 0.0$.
!           IERR = 6 => Matrix A is not Positive Definite.
!                       $(p,Ap) < 0.0$.
!           IERR = 7 => Incomplete factorization broke down
!                       and was fudged.  Resulting preconditioning may
!                       be less than the best.
! IUNIT  :IN       Integer.
!         Unit number on which to write the error at each iteration,
!         if this is desired for monitoring convergence.  If unit
!         number is 0, no writing will occur.
! RWORK  :WORK     Double Precision RWORK(LENW).
!         Double Precision array used for workspace.  NEL is the
!         number of non-
!         zeros in the lower triangle of the matrix (including the
!         diagonal).  NU is the number of nonzeros in the upper
!         triangle of the matrix (including the diagonal).
! LENW   :IN       Integer.
!         Length of the double precision workspace, RWORK.
!         LENW >= NEL+NU+8*N.
! IWORK  :WORK     Integer IWORK(LENIW).
!         Integer array used for workspace.  NEL is the number of non-
!         zeros in the lower triangle of the matrix (including the
!         diagonal).  NU is the number of nonzeros in the upper
!         triangle of the matrix (including the diagonal).
!         Upon return the following locations of IWORK hold information
!         which may be of use to the user:
!         IWORK(9)  Amount of Integer workspace actually used.
!         IWORK(10) Amount of Double Precision workspace actually used.
! LENIW  :IN       Integer.
!         Length of the integer workspace, IWORK.
!         LENIW >= NEL+NU+4*N+12.
            newN = N - 1
            NELT = sparseMatFinal%nnz
            allocate(B(newN), X(newN), IA(NELT), JA(NELT), A(NELT))

            B = -a12
            X = startPd(1:newN)
           
            IA = sparseMatFinal%index(1,:)
            JA = sparseMatFinal%index(2,:)
            A = sparseMatFinal%data
            ISYM = 0
            ITOL = 1
            TOL = 1e-3
            ITMAX = 10000
            IUNIT = 0
            LENW = (9*newN + NELT*2)
            allocate(RWORK(LENW))
            LENIW = (12 + 2*NELT + 5*newN)
            allocate(IWORK(LENIW))
            call DSLUBC(newN, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW)
            if( IERR > 0) then
                print*, ' '
                print*, "IERR: look up: ", IERR
                print*, " "
            endif
            print*, "max in x: ", maxval(X)
            print*, "norm: x", norm2(X)
            print*, "sum: x"
            endPd(1:newN) = x
            if(1>maxval(endPd)) then
                print*, "Tolerance likely too high. Decrease tolerance."
            end if
            endPd(N) = 1
            
            endPd = endPd/sum(endPd)
            numberOfMultiplications = ITER

            print*, "ITERATIONS USED: ", ITER, ". Maximum allowed iterations: ", itmax
            print*, "IWORKSPACE USED: ", IWORK(9), ". Allocated space: ", LENIW
            print*, "DWORKSPACE USED: ", IWORK(10), ". Allocated space: ", LENW 




        end subroutine run_until_converged_lin_system_bicg

        subroutine remove_last_row_and_col_COO(COOIN, COOOUT, lastRow, lastCol, alpha)
            implicit none
            type(COO_dp), intent(inout) :: COOIN, COOOUT
            real(dp), intent(out), allocatable :: lastRow(:), lastCol(:)
            real(dp) :: alpha
            integer :: i, N, newNNz, row, col
            real(dp) :: newData(COOIN%nnz)
            integer :: newIdx(2,COOIN%nnz)
            N = COOIN%ncols
            allocate(lastRow(N-1), lastCol(N-1))
            newNNz = 0
            lastRow = 0
            lastCol = 0


            do i = 1,COOIN%nnz
                row = COOIN%index(1,i)
                col = COOIN%index(2,i)
                if (row == N .AND. col == N) then
                    alpha = COOIN%data(i)
                else if (row == N) then
                    lastRow(col) = COOIN%data(i)
                else if (col == N) then
                    lastCol(row) = COOIN%data(i)
                else
                    newNNz = newNNZ + 1
                    newIdx(1,newNNz) = row
                    Newidx(2,newNNz) = col
                    newData(newNNz) = COOIN%data(i) 
                endif
            end do
            call COOOUT%malloc(N-1,N-1, newNNz)
            COOOUT%index = newIdx(:,1:newNNZ)
            COOOUT%data = newData(1:newNNz)
        end subroutine remove_last_row_and_col_COO

end module sparse_linear_system_solver_bicg_module