module twoDimensionalMarkov
    use fsparse
    implicit none
    public :: run
    
    
    contains
        subroutine run()
            integer , dimension(2) :: coord, dimSize, linkedRow, linkedCol
            integer i, midPointRow, midPointCol
            double precision, dimension(:,:), allocatable :: A, grid
            double precision, dimension(:), allocatable :: probabilityCoeffs

            character (len=100) :: filePath
            filePath = "A.txt"
            dimSize(1) = 21
            dimSize(2) = 21
            midPointRow = dimSize(1)/2 + 1
            midPointCol = dimSize(2)/2 + 1
            allocate(probabilityCoeffs(dimSize(1)*dimSize(2)))
            A = walkMatrix(dimSize)
            
            linkedRow = (/midPointRow,1/)
            linkedCol = (/1, midPointCol/)
            
            
            do i = 1, 3
                coord(2) = linkedRow(2)
                coord(1) = linkedRow(1) + i
                call linkCells(A, coord, getMiddleCell(dimSize), dimSize)
                coord(1) = linkedRow(1) - i
                call linkCells(A, coord, getMiddleCell(dimSize), dimSize)
            end do
            
            do i = 1, 3
                coord(1) = linkedCol(1)
                coord(2) = linkedCol(2) + i
                call linkCells(A, coord, getMiddleCell(dimSize), dimSize)
                coord(2) = linkedCol(2) - i
                call linkCells(A, coord, getMiddleCell(dimSize), dimSize)
                
            end do
            
            
            
            probabilityCoeffs = 0
            probabilityCoeffs(linearIdxFromCoord(getMiddleCell(dimSize), dimSize)) = 1.0 !Initialize middle cell to 1
            
            !timeSteps = 250
            !probabilityCoeffs = timeStep(A, probabilityCoeffs, timeSteps)
            probabilityCoeffs = timeStepUntilConvergence(A, probabilityCoeffs)
            grid = gridFromColumnVector(probabilityCoeffs, dimSize)
            call printMatrix(grid)
            call printMatrixToFile(filePath, grid)
        end subroutine run

        function timeStep(matrix, probabilityCoeffs, steps) !This computes A^k then applies to V
        !MUCH SLOWER THAN TIMESTEP2
            integer, intent(in) :: steps
            double precision, dimension(:), intent(in) :: probabilityCoeffs
            double precision, dimension(SIZE(probabilityCoeffs)) :: timeStep
            double precision, dimension(SIZE(probabilityCoeffs), SIZE(probabilityCoeffs)), intent(in):: matrix 
            double precision, dimension(SIZE(probabilityCoeffs), SIZE(probabilityCoeffs)) :: tempMatrix
            integer :: i
            tempMatrix = matrix
            timeStep = 0
            do i = 2, steps
                tempMatrix = MATMUL(matrix, tempMatrix)
            end do
            timeStep = MATMUL(tempMatrix, probabilityCoeffs)
        end function timeStep
        function timeStep2(matrix, probabilityCoeffs, steps) !This computes repeated matrix multiplication (A*(A*...*(A*V)
            integer, intent(in) :: steps                     !In order to only have to store one matrix
            double precision, dimension(:), intent(in) :: probabilityCoeffs
            double precision, dimension(SIZE(probabilityCoeffs)) :: timeStep2
            double precision, dimension(SIZE(probabilityCoeffs), SIZE(probabilityCoeffs)), intent(in):: matrix 
            integer :: i
            timeStep2 = probabilityCoeffs
            do i = 1, steps
                timeStep2 = MATMUL(matrix, timeStep2)
            end do
        end function timeStep2
        
        
        function timeStepUntilConvergence(matrix, probabilityCoeffs)
            double precision, dimension(:), intent(in) :: probabilityCoeffs
            double precision, dimension(SIZE(probabilityCoeffs), SIZE(probabilityCoeffs)), intent(in) :: matrix
            double PRECISION, dimension(SIZE(probabilityCoeffs)) :: prevProbCoeffs, timeStepUntilConvergence
            double PRECISION :: tol
            logical converged
            integer :: timeSteps, stepsTaken
            
            stepsTaken = 0
            tol = 1.0/(SIZE(probabilityCoeffs) * 1e4) !Tolerance is one part in ten thousand.
            timeSteps = 50
            tol = tol * timeSteps !Dynamically change tolerance based on number of steps taken
            
            converged = .FALSE.
            prevProbCoeffs = probabilityCoeffs
            do while (.not. converged)
                timeStepUntilConvergence = timeStep2(matrix, prevProbCoeffs, timeSteps) !Dont check for convergence after every step, rather take a few time steps at a time.
                converged = convergence(timeStepUntilConvergence, prevProbCoeffs, tol)
                prevProbCoeffs = timeStepUntilConvergence
                stepsTaken = stepsTaken + timeSteps
                print *, 'Steps taken: ', stepsTaken
            end do
            
        end function timeStepUntilConvergence
        
        !TODO
        
        function convergence(newCoeff, oldCoeff, tol)
            double PRECISION, dimension(:), intent(in) :: newCoeff, oldCoeff
            double PRECISION, dimension(SIZE(newCoeff)) :: difference
            double PRECISION, intent(in) :: tol
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
    
        subroutine printMatrixToFile(fileName, A)
            double precision, dimension(:,:), intent(in) :: A
            integer :: i,j
            character (len = *), intent(in) :: fileName
            open(unit = 9, file = fileName)
            
            do i = 1,SIZE(A,1)
                do j = 1,SIZE(A,2)
                    write(9,'(E12.4)', advance='no') A(i,j)
                    
                    if (j /= SIZE(A,2)) then !Dont write semicolon after last element
                        write(9,'(A)', advance='no') '; '  ! Add space between elements except last element
                    end if
                end do
                write(9,*)
            end do
            
            close(unit = 9)
        end subroutine printMatrixToFile
        
        subroutine printMatrix(A)
            double precision, dimension(:,:), intent(in) :: A
            integer :: i, j
            print*, 'Matrix: '
            do i = 1,SIZE(A,1)
                do j = 1,SIZE(A,2)
                    write(*,'(E12.4)', advance='no') A(i,j)
                    write(*,'(A)', advance='no') '  '  ! Add space between elements
                end do
                print*
            end do
            return
        end subroutine printMatrix
        
        subroutine linkCells(A, fromCoord, toCoord, dimSize)
            double precision, dimension(:,:) :: A
            integer, dimension(:) , intent(in) :: fromCoord, toCoord, dimSize
            integer :: fromIdx, toIdx, i, idx
            integer, dimension(:,:), allocatable :: neighbours
            double precision :: prob
            
            fromIdx = linearIdxFromCoord(fromCoord, dimSize)
            toIdx = linearIdxFromCoord(toCoord, dimSize)
            neighbours = pruneNeighbours(getNeighbours(fromCoord), dimSize)
            
            do i = 1, SIZE(neighbours, 2)
                if (all(toCoord == neighbours(:, i))) then
                    return !If toCoord already in list of neighbours, no need to link.
                end if
            end do
            
            prob = 1.0/(Size(neighbours,2) + 1.0)
            A(toIdx, fromIdx) = prob
            do i = 1, SIZE(neighbours,2)
                idx = linearIdxFromCoord(neighbours(:,i), dimSize)
                A(idx, fromIdx) = prob
            end do
        end subroutine linkCells
        
        subroutine linkMultipleCells(A, fromCoords, toCoord, dimSize)
            double precision, dimension(:,:) :: A
            integer, dimension(:) , intent(in) :: toCoord, dimSize
            integer, dimension(:,:), intent(in) :: fromCoords
            integer :: i
            
            do i = 1, SIZE(fromCoords,2)
                call linkCells(A, fromCoords(:,i), toCoord, dimSize)
            end do
        end subroutine linkMultipleCells
        
        function getMiddleCell(dimSize)
            integer, intent(in), dimension(:) :: dimSize
            integer, dimension(SIZE(dimSize)) :: getMiddleCell
            integer :: i
            do i = 1 , SIZE(dimSize)
                getMiddleCell(i) = (dimSize(i)- 1)/2 + 1
            end do
        end function getMiddleCell
        
        function walkMatrix(dimSize) !Creates walk matrix
            integer, intent(in), dimension(:) :: dimSize
            double precision, dimension(:,:), allocatable :: walkMatrix
            integer :: d, i, numberOfCoords, idx
            integer, dimension(:), allocatable :: coord
            integer, dimension(:,:), allocatable :: neighbours
            real :: prob
            
            numberOfCoords = 1
            do d = 1, SIZE(dimSize)
                numberOfCoords = numberOfCoords * dimSize(d)
            end do
            
            allocate(walkMatrix(numberOfCoords, numberOfCoords))
            walkMatrix = 0
            
            do d = 1, numberOfCoords
                coord = coordFromLinearIdx(d, dimSize)
                neighbours = pruneNeighbours(getNeighbours(coord), dimSize)
                prob = 1.0/SIZE(neighbours,2)
                do i = 1, SIZE(neighbours, 2)
                    coord = neighbours(:, i)
                    idx = linearIdxFromCoord(coord, dimSize)
                    walkMatrix(idx, d) = prob
                end do
            end do
        end function walkMatrix
            
        
        function pruneNeighbours(neighbours, dimSize) !Some neighbours are outside bounds. remove them
            integer, intent(in), dimension(:) :: dimSize
            integer, intent(in), dimension(:, :) :: neighbours
            integer, dimension(:, :), allocatable :: pruneNeighbours, tempNeighbours     
            integer :: i, j, validCount, numDims, numCols
            logical :: validNeighbour
            
            validCount = 0
            numDims = size(dimSize)
            numCols = size(neighbours, 2)
            allocate(tempNeighbours(numDims, numCols))
            
            do j = 1, numCols
                validNeighbour = .TRUE.
                do i = 1, numDims
                    if (neighbours(i, j) < 1 .or. neighbours(i,j) > dimSize(i)) then
                        validNeighbour = .FALSE. !exit loop if one coordinate out of bounds
                        exit
                    end if      
                end do
                if (validNeighbour) then
                    validCount = validCount + 1
                    tempNeighbours(:, validCount) = neighbours(:,j)
                end if
            end do
            
            if (validCount > 0) then
                allocate(pruneNeighbours(numDims, validCount))
                pruneNeighbours(:, 1:validCount) = tempNeighbours(:, 1:validCount)
            else
                allocate(pruneNeighbours(numDims, 0))
            end if
        end function pruneNeighbours
            
        
        function getNeighbours(coord)
            integer, intent(in), dimension(:) :: coord
            integer, dimension(SIZE(coord), 2*SIZE(coord) + 1) :: getNeighbours !Includes self
            integer :: i, startIdx
            
            do i = 1, (SIZE(coord)*2 + 1)
                getNeighbours(:,i) = coord
            end do
             
            do i = 1, SIZE(coord)
                startIdx = 2*i
                getNeighbours(i, 2*i) = coord(i) - 1
                getNeighbours(i, 2*i + 1) = coord(i) + 1
            end do
        end function getNeighbours
        
        function coordFromLinearIdx(idx, dimSize)
        !Function converts linear idx to coordinates
            integer, intent(in) :: idx
            integer, dimension(:), intent(in) :: dimSize
            integer, dimension(SIZE(dimSize)) :: coordFromLinearIdx
            
            integer :: i, tempidx, multfactor
            
            tempidx = idx
            multfactor = 1
            
            do i = 1, size(dimSize) !Calculate multiplication factor
                multfactor = multfactor * dimSize(i)
            end do
            
            do i = SIZE(dimSize), 1, -1
                if (i == 1) then
                    coordFromLinearIdx(i) = tempidx
                else
                    multfactor = multfactor / dimSize(i)
                    coordFromLinearIdx(i) = (tempidx - 1) / multfactor + 1
                    tempidx = MOD(tempidx - 1, multfactor) + 1
                endif
            end do
        end function coordFromLinearIdx
        
        integer function linearIdxFromCoord(coord, dimSize)
        !Function gives linear index from n-dimensional coordinate system.
        !First coordinate varies first.
        !I.e. (x,y,z) with number of spots in each dim: (N_x, N_y, N_z) = (3,3,3)
        !Then (3,1,1) -> Idx: 3
        !     (1,2,1) -> Idx: 4
        !     (1,2,2) -> Idx: 13 
            implicit none
            integer, dimension(:), intent(in) :: coord, dimSize
            integer :: i, multFactor

            ! Initialize the index
            linearIdxFromCoord = coord(1)

            ! Handle more than one dimension
            if (SIZE(coord) > 1) then
                multFactor = 1
                do i = 2, SIZE(coord)
                    multFactor = multFactor * dimSize(i-1)
                    linearIdxFromCoord = linearIdxFromCoord + (coord(i) - 1) * multFactor
                end do
            end if
        end function linearIdxFromCoord
        
        function gridFromColumnVector(colVector, dimSize) !Only works for 2 dimensions
            double precision, dimension(:), intent(in) :: colVector
            integer, dimension(2), intent(in) :: dimSize
            integer :: i
            integer, dimension(2) :: coord
            double precision, dimension(dimSize(1),dimSize(2)) :: gridFromColumnVector
            
            do i = 1, SIZE(colVector)
                coord = coordFromLinearIdx(i, dimSize)
                gridfromColumnVector(coord(1), coord(2)) = colVector(i)
            end do
        end function gridFromColumnVector
        

end module twoDimensionalMarkov
    