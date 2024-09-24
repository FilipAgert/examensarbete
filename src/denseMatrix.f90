module denseMatrix
    
    use iso_fortran_env, only: dp=>real64
    implicit none
    contains    
        subroutine printMatrixToFile(fileName, A)
            !Prints two dimensional matrix A to file.
            real(dp), dimension(:,:), intent(in) :: A
            integer :: i,j
            character (len = *), intent(in) :: fileName
            character (len = 5) :: folderPath
            folderPath = 'data/' !Name cannot have prefix /
            folderPath = trim(folderPath) 

            open(unit = 9, file = trim(folderPath//trim(fileName)))
            
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
        
        function getMiddleCell(dimSize)
            integer, intent(in), dimension(:) :: dimSize
            integer, dimension(SIZE(dimSize)) :: getMiddleCell
            integer :: i
            do i = 1 , SIZE(dimSize)
                getMiddleCell(i) = (dimSize(i)- 1)/2 + 1
            end do
        end function getMiddleCell

        subroutine printMatrixD(A)
            !Prints two dimensional matrix A to console output.
            real(dp), dimension(:,:), intent(in) :: A
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
        end subroutine printMatrixD

        subroutine printMatrixS(A)
        !Prints two dimensional matrix A to console output.
            real, dimension(:,:), intent(in) :: A
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
        end subroutine printMatrixS
        
        
        
        
        function walkMatrix(dimSize)
            !Creates matrix containing transition probabilities for every state
            !Gives equal probability for each neighbour 
            integer, intent(in), dimension(:) :: dimSize
            real(dp), dimension(:,:), allocatable :: walkMatrix
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
            
        
        function pruneNeighbours(neighbours, dimSize) 
            !Removes all neighbours outside the bounds given in dimsize
            !Removes neighbours with any coordinate < 1 or larger than dimension given in dimSize
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
            !Gets all neighbours (non diagonal) of a given coordinate.
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
        subroutine linkStates(A, fromCoord, toCoord, dimSize)
            !Links two states for matrix A. Returns nothing, A is directly modified.
            !Modifies probability such that total probability to leave fromCoord is still equal 1.
            real(dp), dimension(:,:) :: A
            integer, dimension(:) , intent(in) :: fromCoord, toCoord, dimSize
            integer :: fromIdx, toIdx, i, idx
            integer, dimension(:,:), allocatable :: neighbours
            real(dp) :: prob
            
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
        end subroutine linkStates
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
        
        function gridFromColumnVector(colVector, dimSize) !Creates a grid of states from a column vector of states. ONly works for two dimensions
            real(dp), dimension(:), intent(in) :: colVector
            integer, dimension(2), intent(in) :: dimSize
            integer :: i
            integer, dimension(2) :: coord
            real(dp), dimension(dimSize(1),dimSize(2)) :: gridFromColumnVector
            
            do i = 1, SIZE(colVector)
                coord = coordFromLinearIdx(i, dimSize)
                gridfromColumnVector(coord(1), coord(2)) = colVector(i)
            end do
        end function gridFromColumnVector
    
end module denseMatrix