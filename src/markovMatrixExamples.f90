module markovMatrixExamples
    use fsparse
    use denseMatrix
    use markovSparse
    use iso_fortran_env, only: dp=>real64
    implicit none
    public compareSparseRunTime, test

    contains
        subroutine compareSparseRunTime()
            integer , dimension(2) :: dimSize
            real(dp), dimension(:,:), allocatable :: A
            real(dp), dimension(:), allocatable :: probabilityCoeffs
            type(COO_dp) :: sparseA
            real :: T1, T2, T3

            character (len=100) :: filePath
            filePath = "A.txt"
            dimSize = 41
            allocate(probabilityCoeffs(dimSize(1)*dimSize(2)))
            A = fissionFusionHole(dimSize)

            probabilityCoeffs = 0
            probabilityCoeffs(1) = 1.0 !Initialize middle cell to 1
            
            !timeSteps = 250
            !probabilityCoeffs = timeStep(A, probabilityCoeffs, timeSteps)
            call cpu_time(T1)
            call dense2coo (A, sparseA)!Convert A to sparse matrix
            probabilityCoeffs = timeStepUntilConvergenceSparse(sparseA, probabilityCoeffs)
            call cpu_time(T2)
            print*, "Sparse matrix execution time: ", T2-T1 , " seconds"
            probabilityCoeffs = 0
            probabilityCoeffs(1) = 1.0 !Initialize one cell to 1
            probabilityCoeffs = timeStepUntilConvergence(A, probabilityCoeffs)
            call cpu_time(T3)

            print*, "Sparse matrix execution time: ", T2-T1 , " seconds"
            print*, "Dense matrix execution time: ", T3-T2, " seconds"
        end subroutine compareSparseRunTime

        subroutine test()
            integer , dimension(2) :: dimSize
            real(dp), dimension(:,:), allocatable :: A, grid
            real(dp), dimension(:), allocatable :: probabilityCoeffs, out
            type(COO_dp) :: sparseA

            character (len=100) :: filePath
            filePath = "fusion-fission.txt"

            dimsize = 41

            allocate(probabilityCoeffs(dimSize(1)*dimSize(2)))
            allocate(out(dimSize(1)*dimSize(2)))
            A = fissionFusionHole(dimSize)

            probabilityCoeffs = 0
            probabilityCoeffs(1) = 1._dp !One cell is initialized to one.
            call dense2coo(A, sparseA)
            out = 0
            out =  timeStepUntilConvergenceSparse(sparseA, probabilityCoeffs)
            grid = gridFromColumnVector(out, dimSize)
            call printMatrixToFile(filePath, grid)
        end subroutine test

        function noHole(dimSize)
            !Returns a transition matrix from diSize with no holes.
            integer, intent(in), dimension(:) :: dimSize
            real(dp), allocatable, dimension(:,:) :: noHole
            integer :: i, numberOfCoordinates

            numberOfCoordinates = 1
            do i = 1, SIZE(dimSize)
                numberOfCoordinates = numberOfCoordinates * dimSize(i)
            end do

            allocate(noHole(numberOfCoordinates, numberOfCoordinates))
            noHole = walkMatrix(dimSize)
        end function noHole

        function fissionFusionHole(dimSize)
            !Creates fission, fusion hole at middle of edge at first and second dimension
            integer, intent(in), dimension(:) :: dimSize
            real(dp), allocatable, dimension(:,:) :: fissionFusionHole, A
            integer, dimension(SIZE(dimSize)):: coord1, coord2
            A = noHole(dimSize)
            coord1 = 1
            coord2 = 1

            coord1(1) = dimSize(1)/2 + 1

            coord2(1) = dimSize(1)/2 + 1
            coord2(2) = dimSize(2)
            call linkStates(A, coord1, getMiddleCell(dimSize), dimSize)
            call linkStates(A, coord2, getMiddleCell(dimSize), dimSize)
            coord1(1) = coord1(1) - 1
            coord2(1) = coord2(1) - 1
            call linkStates(A, coord1, getMiddleCell(dimSize), dimSize)
            call linkStates(A, coord2, getMiddleCell(dimSize), dimSize)
            coord1(1) = coord1(1) + 2
            coord2(1) = coord2(1) + 2
            call linkStates(A, coord1, getMiddleCell(dimSize), dimSize)
            call linkStates(A, coord2, getMiddleCell(dimSize), dimSize)

            fissionFusionHole = A
        end function fissionFusionHole
            
        
            
end module markovMatrixExamples
