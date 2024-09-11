module markovSparse
    use fsparse
    use denseMatrix
    use iso_fortran_env, only: dp=>real64
    implicit none
    public :: run, test
    
    
    contains
        subroutine compareSparseRunTime()
            integer , dimension(2) :: coord, dimSize, linkedRow, linkedCol
            integer i, midPointRow, midPointCol
            real(dp), dimension(:,:), allocatable :: A
            real(dp), dimension(:), allocatable :: probabilityCoeffs
            type(COO_dp) :: sparseA
            real :: T1, T2, T3

            character (len=100) :: filePath
            filePath = "A.txt"
            dimSize(1) = 41
            dimSize(2) = 41
            midPointRow = dimSize(1)/2 + 1
            midPointCol = dimSize(2)/2 + 1
            allocate(probabilityCoeffs(dimSize(1)*dimSize(2)))
            A = walkMatrix(dimSize)
            
            linkedRow = (/midPointRow,1/)
            linkedCol = (/1, midPointCol/)
            
            
            do i = 1, 3
                coord(2) = linkedRow(2)
                coord(1) = linkedRow(1) + i
                call linkStates(A, coord, getMiddleCell(dimSize), dimSize)
                coord(1) = linkedRow(1) - i
                call linkStates(A, coord, getMiddleCell(dimSize), dimSize)
            end do
            
            do i = 1, 3
                coord(1) = linkedCol(1)
                coord(2) = linkedCol(2) + i
                call linkStates(A, coord, getMiddleCell(dimSize), dimSize)
                coord(2) = linkedCol(2) - i
                call linkStates(A, coord, getMiddleCell(dimSize), dimSize)
            end do

            probabilityCoeffs = 0
            probabilityCoeffs(linearIdxFromCoord(getMiddleCell(dimSize), dimSize)) = 1.0 !Initialize middle cell to 1
            
            !timeSteps = 250
            !probabilityCoeffs = timeStep(A, probabilityCoeffs, timeSteps)
            call cpu_time(T1)
            call dense2coo (A, sparseA)!Convert A to sparse matrix
            probabilityCoeffs = timeStepUntilConvergence(sparseA, probabilityCoeffs)
            call cpu_time(T2)
            print*, "Sparse execution time: ", T2-T1 , " seconds"
            probabilityCoeffs = 0
            probabilityCoeffs(linearIdxFromCoord(getMiddleCell(dimSize), dimSize)) = 1.0 !Initialize middle cell to 1
            probabilityCoeffs = timeStepUntilConvergenceNonSparse(A, probabilityCoeffs)
            call cpu_time(T3)

            
            print*, "Classic execution time: ", T3-T2, " seconds"
        end subroutine compareSparseRunTime

        subroutine test()
            integer , dimension(2) :: coord, dimSize, linkedRow, linkedCol
            integer i, midPointRow, midPointCol
            real(dp), dimension(:,:), allocatable :: A, grid
            real(dp), dimension(:), allocatable :: probabilityCoeffs, out
            type(COO_dp) :: sparseA

            character (len=100) :: filePath
            filePath = "As.txt"

            dimSize(1) = 31
            dimSize(2) = 31
            midPointRow = dimSize(1)/2 + 1
            midPointCol = dimSize(2)/2 + 1
            allocate(probabilityCoeffs(dimSize(1)*dimSize(2)))
            allocate(out(dimSize(1)*dimSize(2)))
            A = walkMatrix(dimSize)

            
            linkedRow = (/midPointRow,1/)
            linkedCol = (/1, midPointCol/)
            
            do i = 1, 3
                coord(2) = linkedRow(2)
                coord(1) = linkedRow(1) + i
                call linkStates(A, coord, getMiddleCell(dimSize), dimSize)
                coord(1) = linkedRow(1) - i
                call linkStates(A, coord, getMiddleCell(dimSize), dimSize)
            end do
            
            do i = 1, 3
                coord(1) = linkedCol(1)
                coord(2) = linkedCol(2) + i
                call linkStates(A, coord, getMiddleCell(dimSize), dimSize)
                coord(2) = linkedCol(2) - i
                call linkStates(A, coord, getMiddleCell(dimSize), dimSize)
            end do
            probabilityCoeffs = 0
            probabilityCoeffs(1) = 1._dp !Initialize middle cell to 1
            call dense2coo(A, sparseA)
            out = 0
            out =  timeStepUntilConvergence(sparseA, probabilityCoeffs)

            grid = gridFromColumnVector(out, dimSize)

            call printMatrix(grid)
            call printMatrixToFile(filePath, grid)

        end subroutine test
        subroutine run()
            integer , dimension(2) :: coord, dimSize, linkedRow, linkedCol
            integer i, midPointRow, midPointCol
            real(dp), dimension(:,:), allocatable :: A, grid
            real(dp), dimension(:), allocatable :: probabilityCoeffs
            type(COO_dp) :: sparseA

            character (len=100) :: filePath
            filePath = "A.txt"
            dimSize(1) = 31
            dimSize(2) = 31
            midPointRow = dimSize(1)/2 + 1
            midPointCol = dimSize(2)/2 + 1
            allocate(probabilityCoeffs(dimSize(1)*dimSize(2)))
            A = walkMatrix(dimSize)
            
            linkedRow = (/midPointRow,1/)
            linkedCol = (/1, midPointCol/)
            
            
            do i = 1, 3
                coord(2) = linkedRow(2)
                coord(1) = linkedRow(1) + i
                call linkStates(A, coord, getMiddleCell(dimSize), dimSize)
                coord(1) = linkedRow(1) - i
                call linkStates(A, coord, getMiddleCell(dimSize), dimSize)
            end do
            
            do i = 1, 3
                coord(1) = linkedCol(1)
                coord(2) = linkedCol(2) + i
                call linkStates(A, coord, getMiddleCell(dimSize), dimSize)
                coord(2) = linkedCol(2) - i
                call linkStates(A, coord, getMiddleCell(dimSize), dimSize)
                
            end do
            
            
            
            probabilityCoeffs = 0
            probabilityCoeffs(linearIdxFromCoord(getMiddleCell(dimSize), dimSize)) = 1.0 !Initialize middle cell to 1
            
            !timeSteps = 250
            !probabilityCoeffs = timeStep(A, probabilityCoeffs, timeSteps)
            call dense2coo (A, sparseA)!Convert A to sparse matrix
            deallocate(A)


            probabilityCoeffs = timeStepUntilConvergence(sparseA, probabilityCoeffs)
            grid = gridFromColumnVector(probabilityCoeffs, dimSize)
            call printMatrix(grid)
            call printMatrixToFile(filePath, grid)
        end subroutine run

        function timeStepSparse(matrixSparse, probabilityCoeffs, steps)
            integer, intent(in) :: steps
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
        end function timeStepSparse
        
        
        function timeStepUntilConvergence(matrix, probabilityCoeffs)
            real(dp), dimension(:), intent(in) :: probabilityCoeffs
            type(COO_dp) :: matrix
            real(dp), dimension(SIZE(probabilityCoeffs)) :: prevProbCoeffs, timeStepUntilConvergence
            real(dp) :: tol
            logical converged
            integer :: timeSteps, stepsTaken
            
            stepsTaken = 0
            tol = 1.0/(SIZE(probabilityCoeffs) * 1e6) !Tolerance is one part in one million.
            timeSteps = 20
            tol = tol * timeSteps !Dynamically change tolerance based on number of steps taken
            
            converged = .FALSE.
            prevProbCoeffs = probabilityCoeffs
            do while (.not. converged)
                timeStepUntilConvergence = timeStepSparse(matrix, prevProbCoeffs, timeSteps) !Dont check for convergence after every step, rather take a few time steps at a time.
                converged = convergence(timeStepUntilConvergence, prevProbCoeffs, tol)
                prevProbCoeffs = timeStepUntilConvergence
                stepsTaken = stepsTaken + timeSteps
                if(mod(stepsTaken, 500) == 0) then
                    print *, 'Steps taken: ', stepsTaken
                end if
            end do
            print*, 'Steps taken: ', stepsTaken
            
        end function timeStepUntilConvergence


        function timeStep(matrix, probabilityCoeffs, steps) !This computes repeated matrix multiplication (A*(A*...*(A*V)
            integer, intent(in) :: steps                     !In order to only have to store one matrix
            double precision, dimension(:), intent(in) :: probabilityCoeffs
            double precision, dimension(SIZE(probabilityCoeffs)) :: timeStep
            double precision, dimension(SIZE(probabilityCoeffs), SIZE(probabilityCoeffs)), intent(in):: matrix 
            integer :: i
            timeStep = probabilityCoeffs
            do i = 1, steps
                timeStep = MATMUL(matrix, timeStep)
            end do
            
        end function timeStep

        function timeStepUntilConvergenceNonSparse(matrix, probabilityCoeffs)
            real(dp), dimension(:), intent(in) :: probabilityCoeffs
            real(dp), dimension(SIZE(probabilityCoeffs),SIZE(probabilityCoeffs)), intent(in):: matrix
            real(dp), dimension(SIZE(probabilityCoeffs)) :: prevProbCoeffs, timeStepUntilConvergenceNonSparse
            real(dp) :: tol
            logical converged
            integer :: timeSteps, stepsTaken
            
            stepsTaken = 0
            tol = 1.0/(SIZE(probabilityCoeffs) * 1e6) !Tolerance is one part in one million.
            timeSteps = 20
            tol = tol * timeSteps !Dynamically change tolerance based on number of steps taken
            
            converged = .FALSE.
            prevProbCoeffs = probabilityCoeffs
            do while (.not. converged)
                timeStepUntilConvergenceNonSparse = timeStep(matrix, prevProbCoeffs, timeSteps) !Dont check for convergence after every step, rather take a few time steps at a time.
                converged = convergence(timeStepUntilConvergenceNonSparse, prevProbCoeffs, tol)
                prevProbCoeffs = timeStepUntilConvergenceNonSparse
                stepsTaken = stepsTaken + timeSteps
                if(mod(stepsTaken, 500) == 0) then
                    print *, 'Steps taken: ', stepsTaken
                end if
            end do
            print*, 'Steps taken: ', stepsTaken
            
        end function timeStepUntilConvergenceNonSparse
        
        !TODO
        
        function convergence(newCoeff, oldCoeff, tol)
            real(dp), dimension(:), intent(in) :: newCoeff, oldCoeff
            real(dp), dimension(SIZE(newCoeff)) :: difference
            real(dp), intent(in) :: tol
            LOGICAL :: convergence
            integer i
            convergence = .TRUE.
            difference = newCoeff - oldCoeff
            do i = 1, SIZE(difference)
                if(difference(i) > tol) then
                    convergence = .FALSE.
                    return
                end if
            end do
        end function convergence
end module markovSparse