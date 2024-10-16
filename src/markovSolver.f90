module markovSolver
    use potsurf
    use sparseMatrix
    use fsparse
    use result_module
    !use sparseFromPotSurf
    implicit none

    logical :: useFullMMCoordinates
    logical :: connectedToStartingCoord = .FALSE.
    integer, dimension(5) :: dimSize = [1 + MAX_II-MIN_II, 1 + MAX_JJ-MIN_JJ,&
                                        1 + MAX_KK-MIN_KK, 1 + MAX_LL-MIN_LL, 1 + MAX_MM - MIN_MM]
    integer, dimension(5,2) :: MIN_MAX_DIM = reshape([MIN_II, MIN_JJ, MIN_KK, MIN_LL, MIN_MM, &
                                        MAX_II, MAX_JJ, MAX_KK, MAX_LL, MAX_MM], [5,2])
    integer, allocatable :: fissionFusionIndices(:), fissionIdxs(:), fusionIdxs(:)
    integer :: AZ, AA
    real(kind = r_kind) :: Egs !!ground state energy.


    !!!!!these are constants we can change
    integer :: II_fusion = 3
    real(kind=r_kind) :: TOL
    integer :: NCV, num_threads
    real(kind=r_kind) :: Rneck_fission = 1.5_r_kind !Fm
    real(kind=r_kind) :: fusionChance = 0.5_r_kind
    real(kind=r_kind) :: fissionChance = 0.5_r_kind
    !!!!!!



    real(kind=r_kind), allocatable :: startingEnergies(:)
    integer, allocatable :: startingCoords(:,:)
    
    type(CSR_dp) :: markovMatCSR
    type(COO_dp) :: markovMatCOO
    
    type(convergedResult) :: results


    contains

    subroutine solveAllEnergies()
        integer :: i, coord(5)
        real(kind=r_kind) :: E, RATE, TIME1, TIME2, fusionFraction, fissionFraction, Eexc
        integer(8) :: COUNT1, COUNT2, COUNT3
        integer :: numberOfMatvecCalls
        real(kind=r_kind), allocatable :: guessPD(:), SOL(:)
        allocate(guessPD(numberOfGridPoints()), SOL(numberOfGridPoints()))

        CALL system_clock(count_rate = rate)
        
        do i = 1,SIZE(startingEnergies)
            CALL system_clock(count=COUNT1)
            coord = startingCoords(:,i)
            E = startingEnergies(i)
            Eexc = E - Egs
            print*, "Energy: ", Eexc, " MeV"
            print*, "Starting coord: ", coord
            call setupMatrix(E,coord)  !Calculate transition matrix
            CALL system_clock(count=COUNT2)
            TIME1 = (count2-count1)/rate
            print*, ' '
            print*, "Matrix setup took ", TIME1, " seconds"

            guessPD = startingGuess(E, markovMatCSR)

            print*, "------------------Starting arnoldi iteration--------------------"
            call findEigenVector(markovMatCSR, guessPD, SOL, numberOfMatvecCalls, TOL, NCV, num_threads)
            fusionFraction = getFusionFraction(SOL)
            fissionFraction = 1 - fusionFraction
            



            call system_clock(count = COUNT3)
            TIME2 = (count3-count2)/rate

            print*, "Finding eigenvector took ", TIME2, " seconds"
            print*, "Setup matrix: ", 100*TIME1/(TIME1+TIME2) , " % of total time"
            print*, "Total time: ", TIME1+TIME2, " seconds"
            call results%addResult(SOL, real(coord,r_kind), Eexc, TIME1+TIME2, numberOfMatvecCalls,fusionFraction,fissionFraction)
            call results%printResult()
            
        end do
    end subroutine solveAllEnergies

    subroutine findEigenVector(sparseMat, startPd, endPd, numberOfMatvecCalls, tol, NCV, num_threads)
            
        implicit none
        type(CSR_dp) :: sparseMat
        real(kind=r_kind), dimension(:), intent(in) :: startPd
        real(kind=r_kind), dimension(:), intent(inout) :: endPd
        real(kind=r_kind), dimension(SIZE(startPd)) :: temp
        integer, intent(inout) :: numberOfMatvecCalls
        logical :: hasNegative
        real(kind=r_kind) :: negP, posP

        external :: dnaupd
        external :: dneupd
        
        integer :: IDO, i
        character(len = 1) :: BMAT = 'I'
        integer :: N
        character(len = 2) ::  WHICH = 'LR'  !since multiple eigenvalues can have abs(lambda)=1, we need to choose only largest real part
        integer :: NEV!number of eigenvalues comptued
        real(kind=r_kind) :: TOL 
        real(kind=r_kind),allocatable :: RESID(:)
        integer :: NCV !number of vectors in calculation. computation scales as N * NCV²
        real(kind=r_kind), allocatable :: V(:,:)
        integer :: IPARAM(11)
        integer :: ishfts
        integer :: maxitr !max iterations?
        integer :: mode = 1 !1 = no sigma shift. 3 = sigma shift of one
        integer :: LDV
        Real(kind=r_kind), allocatable :: workd(:), workl(:), workev(:)
        integer :: LWORKL
        integer :: INFO
        integer :: INFO2 !0 = randomized initial vector, 1 = resid initial vector (use starting guess)
        logical :: converged
        integer :: IPNTR(14)
        LOGICAL :: RVEC
        character(len = 1) :: HOWMNY = 'A'
        logical, allocatable :: SELECT(:)
        real(kind=r_kind), dimension(:), allocatable :: DR, DI
        real(kind=r_kind), dimension(:,:), allocatable :: Z
        real(kind=r_kind) :: SIGMAR, SIGMAI
        real :: TIME1, TIME2, TIME3, RATE
        integer :: COUNT0, COUNT1, COUNT2, COUNT3
        integer :: LDZ
        integer :: num_threads
        N = numberOfGridPoints()
        !!USER SETTINGS.
        maxitr = 50000
        NEV = 1 !number of eigenvalues calculated. test to change this
        !TOL = 1./(numberOfGridPoints()*10)!Tolerance is: how close to true eigenvalue can the calculated one be?
        !Tolerance means that you expect no eigenvalues to be closer to eachother than tolerance.
        !seems to only converge to true value if low enough. Should maybe be of order 1/number of grid points For 501x501

        !NCV = MAX(2*NEV + 1,int((N)**(1.0/8))) !set at least to 2*NEV + 1. Lower number = more matrix * vector operations.
        print*, "NCV: ", NCV
        print*, "TOL:", TOL
        print*, "Tol should be: ", 1.0/(N*10)
        !!!!!!
        SIGMAR = 0 
        SIGMAI = 0
        INFO = 1 !info for dnaupd          1 = user specified vector
        INFO2 = 0 !info for dneupd         
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
        numberOfMatvecCalls = 0
        converged = .FALSE.
        CALL system_clock(count_rate=RATE)
        call system_clock(count=count0)
        do while(.not. converged)
            CALL system_clock(count = count1)
            call dnaupd(IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV,IPARAM,IPNTR, workd, workl, lworkl,info)
            
            if(IDO .eq. -1 .or. IDO .eq. 1) then
                temp = 0

                call system_clock(count=count2)

                if(num_threads > 1) then !Probably dont need this but not sure how OMP works if num_threads = 1 or num_threads < 1
                    call matvecpar(sparseMat, workd(ipntr(1):ipntr(1) + N - 1),temp, num_threads)
                else
                    call matvecpar(sparseMat, workd(ipntr(1):ipntr(1) + N - 1),temp, 1)
                endif

                workd(ipntr(2) : ipntr(2) + N - 1) = temp
                numberOfMatvecCalls = numberOfMatvecCalls + 1

                call system_clock(count=count3)
                if(mod(numberOfMatvecCalls,100) == 0) then
                    Time3 = (count3-count0)/RATE
                    Time2 = (count2-count1)/RATE
                    Time1 = (count3-count2)/rate
                    print*, " "
                    print*, "Number of matvec calls so far : ", numberOfMatvecCalls
                    print*, "Elapsed time so far: ", TIME3, " seconds."
                    print*, "Latest arnoldi iteration time: ", TIME2, " seconds. "
                    print*, "Latest matvec time: ", TIME1, " seconds. ", 100*(TIME1)/(TIME1+TIME2), "% of total"
                endif
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
                endPd = Z(:,1)/sum(Z(:,1))! !Ensure correct phase and amplitude

                !Check if eigenvector has negative and positive values. Should only have positive values.
                hasNegative = .FALSE.
                do i = 1, SIZE(endPd)
                    if(endPd(i) < 0.0_r_kind) then
                        hasNegative = .TRUE.
                        exit
                    endif
                end do
                if(hasNegative) then
                    negP = 0.0_r_kind
                    posP = 0.0_r_kind
                    do i = 1,SIZE(endPd)
                        if(endPd(i) < 0.0_r_kind) then
                            negP = negP + endPd(i)
                        else
                            posP = posP + endPd(i)
                        endif
                    end do
                    print*, " "
                    print*, "--------------------------------------------------------------------------------------------------"
                    print*, "ERROR: Eigenvector contains both negative and positive probability values."
                    print*, "Tolerance is likely too large."
                    print*, "Sum probability positive entries (should be 1): ", posP
                    print*, "Sum probability negative entries: ", negP
                    print*, "--------------------------------------------------------------------------------------------------"
                    print*, " "
                endif
                print *, 'Found eigenvector'
                print *, 'Eigenvalue:', DR(1)
                print*, "Minimum value of eigenvector: ", minval(endPd)
                print*, "Maximum value of eigenvector: ", maxval(endPd)
                
            endif
        endif
        continue
    end subroutine findEigenVector

    subroutine setupSolver(tolerance, arnoldiNCV, number_of_threads, Z, A, Rneck_fis, II_fus, Eg, startingEs, startingC, &
        useFullMM, filename_emac5D, filename_pot5D, filename_rneck5D)
        integer :: Z,A, II_fus
        real(kind=r_kind) :: Rneck_fis, startingEs(:), tolerance, Eg
        integer :: startingC(:,:), arnoldiNCV, number_of_threads
        logical :: useFullMM
        CHARACTER(100) :: filename_pot5D                   ! Total potential energy (Emac + Emic)
        CHARACTER(100) :: filename_emac5D                     ! Macroscopic potential energy Emac
        CHARACTER(100) :: filename_rneck5D
        call setMMsymmetryMode(useFullMM)
        Egs = Eg
        TOL = tolerance
        NCV = arnoldiNCV
        num_threads = number_of_threads
        AZ = Z
        AA = A
        Rneck_fission = Rneck_fis
        II_fusion = II_fus
        call read_5d_emac(AZ,AA,filename_emac5D)
        call read_5D_epot(AZ,AA,filename_pot5D)
        call read_5d_rneck(AZ,AA,filename_rneck5D) 

        fusionIdxs = getFusionIndices(II_fusion)
        fissionIdxs = getFissionIndices(Rneck_fission)
        print*, "number of fission indicies", SIZE(fissionIdxs)
        allocate(fissionFusionIndices(size(fusionIdxs) + size(fissionIdxs)))
        allocate(startingEnergies(size(startingEs)))
        allocate(startingCoords(5, size(startingC,2)))
        startingEnergies = startingEs
        startingCoords = startingC
        

    end subroutine setupSolver

    subroutine setupMatrix(Etot, startingCoord)
        real(kind=r_kind) :: Etot
        integer, dimension(5) :: startingCoord
        type(CSR_DP) :: CSR
        type(COO_DP) :: COO

        print*, "Calculate matrix from potential... "
        COO = sparseFromPotential(AZ, AA, Etot, II_fusion, fusionChance, Rneck_fission, fissionChance, dimSize, &
                                            MIN_MAX_DIM, useFullMMCoordinates)
        !Sets upp connection between fusion/fission coordinates to startnig coordinate
        call connectToStartingCoord(COO, startingCoord, connectedToStartingCoord, fusionIdxs, fissionIdxs, &
                                    fissionFusionIndices, fusionChance, fissionChance, dimSize, MIN_MAX_DIM) 
        call coo2ordered(COO, .TRUE.)
        call coo2csr(COO, CSR)
        markovMatCSR = CSR

        print*, "Calculated matrix from potential..."
    end subroutine setupMatrix


    subroutine setMMsymmetryMode(useFullMM)
        !!Set if using halved MM coordinates or full MM coordinates. Default is full.
        logical :: useFullMM !useFullMM = true => we use symmetry and cut half MM coordinates
        useFullMMCoordinates = useFullMM
        if(.NOT. useFullMM) then
            dimSize(5) = 1 + MAX_MM
            MIN_MAX_DIM(5,1) = 0
        else
            dimSize(5) = 1 + MAX_MM - MIN_MM
            MIN_MAX_DIM(5,1) = MIN_MM
        endif
    end subroutine setMMsymmetryMode

    function numberOfGridPoints() result(n)
        integer :: n
        integer :: i
        n = 1
        do i = 1,size(dimSize)
            n = n*dimSize(i)
        end do
    end function numberOfGridPoints
    

    function startingGuess(energy, CSR) result(guess)
        type(CSR_dp) :: CSR
        real(r_kind), allocatable :: guess(:)
        integer :: storedResults, i, E0idx, E1idx, tempidx
        real(r_kind) :: energy, E0, E1, diff, potentialClosest, tempE, E1DIFF, E0DIFF
        real(r_kind), allocatable :: pdf0(:), pdf1(:) 
        allocate(guess(numberOfGridPoints()))
        guess = 0
        storedResults = results%numResults
        if(storedResults == 0) then !Guess filling bathtub
            !guess = guessFromPotential(energy) 
            guess = 1.0/numberOfGridPoints()
            print*, "GUESS = even distribution"
        elseif(storedResults == 1) then
            guess = results%getProbabilityDensity(storedResults)
            print*, "GUESS = RE-USE 1st solution"
        else
            E0 = huge(E0)
            E1 = huge(E1) - 1

            E0idx = storedResults - 1
            E1idx = storedResults
            do i = 1, storedResults !find two closest energies as they will produce best fit.
                potentialClosest = results%energies(i)
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
            
        
            allocate(pdf0(numberOfGridPoints()), pdf1(numberOfGridPoints()))

            pdf0 = results%getProbabilityDensity(E0idx)
            pdf1 = results%getProbabilityDensity(E1idx)
            guess = linearInterpolation(energy, E0, E1, pdf0, pdf1)!
            print*, "Guess = LINCOMB of energies:", E0, " and ", E1, " MeV"
        endif
        guess = applyMatrix(maxval(dimSize), guess, CSR)
    end function startingGuess

    function applyMatrix(iterations, input, CSR) result(res)
        !!This applies the CSR matrix iteratively on input and returns resulting vector
        !!This is used when generating starting guess. First guess starting guess,
        !!Then apply matrix a few times to refine guess. This is because the change is quite large in the beginning
        !!And maybe arnodli will perform better with a better starting guess.
        type(CSR_dp) :: CSR
        integer :: iterations, i
        real(kind=r_kind) :: input(:)
        real(kind=r_kind), dimension(size(input)) :: res, temp
        res = input
        print*, "Size: ", sum(res)
        do i = 1,iterations
            temp = 0
            call matvecpar(CSR, res, temp, num_threads)
            res = temp
            
        end do
        
        

    end function applyMatrix

    function guessFromPotential(Etot) result(guess)
        !!Function guesses initial distribution by giving each point a higher weight if the potential energy is lower there.
        !!potsurf must have read total energy for this to work.
        !!Sort of like filling a bathtub and giving weight from depth.
        real(r_kind) :: Etot !Total energy of particles
        real(r_kind), allocatable :: guess(:)
        real(r_kind) :: Ep, height, rollingSum, maxHeight, minHeight
        integer :: idx, coord(5)
        integer :: II,JJ,KK,LL,MM
        allocate(guess(numberOfGridPoints()))
        rollingSum = 0.0_r_kind
        maxHeight = 0
        minHeight = huge(minHeight)
        do II = MIN_II, MAX_II
            do JJ = MIN_JJ, MAX_JJ
                do KK = MIN_KK, MAX_KK
                    do LL = MIN_LL, MAX_LL
                        do MM = MIN_MAX_DIM(5,1), MAX_MM
                            coord = [II,JJ,KK,LL,MM]
                            idx = linearIdxFromCoord(coord, dimSize, MIN_MAX_DIM)
                            Ep = Epot(II,JJ,KK,LL,MM)
                            height = Etot - Ep
                            if(height > 0.0_r_kind) then
                                guess(idx) = height
                                rollingSum = rollingSum + height
                                ! if(height > maxHeight) then
                                !     maxHeight = height
                                ! endif
                                ! if(height < minHeight) then
                                !     minHeight = height
                                ! endif
                            else
                                guess(idx) = 0.0_r_kind
                            endif
                        end do
                    end do
                end do
            end do
        end do
        guess = guess/rollingSum !Normalize to sum 1.
        ! print*, "Maxheight: ", maxHeight
        ! print*, "Minheight: ", minHeight
    end function guessFromPotential



    function linearInterpolation(E,E0,E1, pdf0, pdf1) result(interpolation)
        !!E is energy for distribution to be guessed.
        !!E0 and E1 correspond to energies of pdf0 and pdf1 which are previous calculations.
        real(r_kind), intent(in) ::  E, E0, E1
        real(r_kind), intent(in) :: pdf0(:), pdf1(:)
        real(r_kind) :: interpolation(size(pdf0))
        real(r_kind) :: alpha, sumInt
        alpha = (E-E0)/(E1-E0)
        interpolation = (1-alpha)*pdf0 + alpha*pdf1
        sumInt= sum(interpolation)
        interpolation = interpolation/sumInt
    end function linearInterpolation

    

    function getFusionIndices(II_fusion) result(fusionIndices) !Fusion for II less than or equal II_fusion
        integer, intent(in) :: II_fusion
        integer, dimension(:), allocatable :: fusionIndices
        integer ::  II, JJ, KK, LL ,MM, NumberOfFusionIndices, i
        integer, dimension(5) :: coord

        NumberOfFusionIndices = (MAX_MM- MIN_MAX_DIM(5,1) + 1)*(MAX_LL-MIN_LL + 1)*(MAX_KK-MIN_KK + 1)*(MAX_JJ-MIN_JJ + 1)
        NumberOfFusionIndices = NumberOfFusionIndices * (II_fusion - MIN_II + 1) !Calculates number of fusion indices

        allocate(fusionIndices(NumberOfFusionIndices))
        i = 1
        do II = MIN_II, II_fusion
            do JJ = MIN_JJ, MAX_JJ
                do KK = MIN_KK, MAX_KK
                    do LL = MIN_LL, MAX_LL
                        do MM = MIN_MAX_DIM(5,1), MAX_MM
                            coord = [II, JJ, KK, LL, MM]
                            fusionIndices(i) = linearIdxFromCoord(coord, dimSize, MIN_MAX_DIM)
                            i = i + 1
                        end do
                    end do
                end do
            end do
        end do
    end function getFusionIndices

    function getFissionIndices(Rneck_fission) result(fissionIndices) !Fusion for II less than or equal II_fusion
        integer, dimension(:), allocatable :: fissionIndices
        integer ::  II, JJ, KK, LL ,MM, NumberOfFissionIndices, i, other
        integer, dimension(5) :: coord
        real(kind = r_kind), intent(in) :: Rneck_fission
        other = 0
        NumberOfFissionIndices = 0
        do II = MIN_II, MAX_II
            do JJ = MIN_JJ, MAX_JJ
                do KK = MIN_KK, MAX_KK
                    do LL = MIN_LL, MAX_LL
                        do MM = MIN_MAX_DIM(5,1), MAX_MM
                            if(Rneck(II,JJ,KK,LL,MM) < Rneck_fission) then
                                NumberOfFissionIndices = NumberOfFissionIndices + 1
                            else
                                other = other + 1
                            endif
                        end do
                    end do
                end do
            end do
        end do
        allocate(fissionIndices(NumberOfFissionIndices))
        i = 0
        do II = MIN_II, MAX_II
            do JJ = MIN_JJ, MAX_JJ
                do KK = MIN_KK, MAX_KK
                    do LL = MIN_LL, MAX_LL
                        do MM = MIN_MAX_DIM(5,1), MAX_MM
                            if(Rneck(II,JJ,KK,LL,MM) < Rneck_fission) then
                                i = i + 1
                                coord = [II, JJ, KK, LL, MM]
                                fissionIndices(i) = linearIdxFromCoord(coord, dimSize, MIN_MAX_DIM)
                            endif
                        end do
                    end do
                end do
            end do
        end do
    end function getFissionIndices

    function getFusionFraction(pd) result (fusionFrac)
        real(r_kind), intent(in) :: pd(:)
        real(r_kind) :: fusionFrac, fusionProb, fissionProb
        integer :: i
        fusionFrac = 0
        fusionProb = 0
        fissionProb = 0
        do i = 1, SIZE(fusionIdxs)
            fusionProb = fusionProb + (pd(fusionIdxs(i)))
        end do

        do i = 1, SIZE(fissionIdxs)
            fissionProb = fissionProb + (pd(fissionIdxs(i)))
        end do

        print*, "Sum of fusionprob: ", fusionProb
        print*, "Sum of fissionprob: ", fissionProb

        fusionFrac = fusionProb/(fusionProb + fissionProb)
    end function getFusionFraction
end module markovSolver