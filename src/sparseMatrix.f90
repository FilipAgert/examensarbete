module sparseMatrix
    use fsparse
    use iso_fortran_env, only: dp=>real64, sp=>real32
    use denseMatrix, only: linearIdxFromCoord, coordFromLinearIdx, pruneNeighbours, getNeighbours
    implicit none
    private

    public :: sparseWalkMatrix, linkStates, linkStatesIdx
    contains

    function sparseWalkMatrix(dimSize)result(sparse)
        integer, dimension(:) :: dimSize
        type(COO_sp) :: sparse
        integer :: numCoords, i, j, neighboursPerState, coord(SIZE(dimSize)), neighbourIdx, idx
        integer, dimension(:,:), allocatable :: neighbours
        real(sp) :: prob
        
        neighboursPerState = 1 + 2*SIZE(dimSize)

        numCoords = 1
        do i = 1, size(dimSize)
            numCoords = numCoords * dimSize(i) !Find total number of coordinates
        end do

        call sparse%malloc(numCoords, numCoords, neighboursPerState*numCoords) !This will allocate a bit too many non-zero elements since edges dont have as many neighbours
        idx = 1
        do i = 1,numCoords
            coord = coordFromLinearIdx(i, dimSize)

            neighbours = pruneNeighbours(getNeighbours(coord), dimSize)
            prob = 1.0/size(neighbours,2)
            do j = 1, size(neighbours,2)
                coord = neighbours(:,j)
                neighbourIdx = linearIdxFromCoord(coord,dimSize)
                call addVal(sparse, neighbourIdx, i, idx, prob)
                idx = idx + 1
            end do
        end do
        sparse%isOrdered = .FALSE.
        call changeSize(sparse, idx - 1)!Resize to smaller size.
    end function sparseWalkMatrix

    subroutine changeSize(sparseMat, new_nnz)
        type(COO_sp) :: sparseMat
        integer, intent(in) :: new_nnz
        integer :: current_nnz
        integer, allocatable :: temp_idx(:,:)
        real(sp), allocatable :: temp(:)
        current_nnz = sparseMat%nnz

        if (current_nnz < new_nnz) then
            ! Create temporary arrays for resizing
            
        
            ! Allocate the larger arrays
            allocate(temp_idx(2, new_nnz))
            allocate(temp(new_nnz))
        
            ! Copy old values to the larger arrays
            temp_idx(:, 1:current_nnz) = sparseMat%index(:, 1:current_nnz)
            temp(1:current_nnz) = sparseMat%data(1:current_nnz)
        
            ! Optionally initialize the new entries
            temp_idx(:, current_nnz+1:new_nnz) = 0
            temp(current_nnz+1:new_nnz) = 0.0_sp
        
            ! Assign the new arrays back to the COO matrix
            sparseMat%index = temp_idx
            sparseMat%data = temp
        else 
            ! Shrink the arrays before calling malloc_coo
            sparseMat%index = sparseMat%index(:, 1:new_nnz)  ! Truncate index array
            sparseMat%data  = sparseMat%data(1:new_nnz)      ! Truncate data array
        end if
        ! Now it's safe to call malloc_coo with the new nnz
        call sparseMat%malloc(sparseMat%nrows, sparseMat%ncols, new_nnz) 

    end subroutine changeSize

    subroutine addVal(sparseMat, i, j, idx, val)
        type(COO_sp):: sparseMat
        integer :: i,j, idx
        real(sp) :: val

        sparseMat%index(1,idx) = i
        sparseMat%index(2,idx) = j
        sparseMat%data(idx) = val
    end subroutine addVal
        
    subroutine linkStates(sparseMat, fromCoord, toCoord, dimSize) 
        type(COO_sp), intent(inout) :: sparseMat
        integer, dimension(:),intent(in) :: fromCoord, toCoord
        integer, dimension(:), intent(in) :: dimSize
        integer :: nnz, numRows, numCols, fromIdx, toIdx, i
        integer, dimension(:,:), allocatable :: neighbours
        real(sp) :: prob
        numRows = sparseMat%nrows
        numCols = sparseMat%ncols
        nnz = sparseMat%nnz
        
        neighbours = pruneNeighbours(getNeighbours(fromCoord),dimSize)

        do i = 1, SIZE(neighbours, 2)
            if (all(toCoord == neighbours(:, i))) then
                return !If toCoord already in list of neighbours, no need to link.
            end if
        end do
        call changeSize(sparseMat, nnz + 1)
        
        prob = 1.0/(Size(neighbours,2) + 1.0)
        fromIdx = linearIdxFromCoord(fromCoord,dimSize)
        toIdx = linearIdxFromCoord(toCoord, dimSize)

        call addVal(sparseMat,toIdx, fromIdx, nnz + 1, prob)
        do i = 1, SIZE(neighbours,2)
            toIdx = linearIdxFromCoord(neighbours(:,i), dimSize)
            call sparseMat%set(prob, toIdx, fromIdx)
        end do
    end subroutine linkStates

    subroutine linkStatesIdx(sparseMat, fromIdx, toIdx, val)
        type(COO_sp), intent(inout) :: sparseMat
        integer, intent(in) :: fromIdx, toIdx
        real(sp) :: out, val
    
        out = 0.0_sp
        call sparseMat%get(out, toIdx, fromIdx)
        if(out /= 0.0_sp) then !if already linked, do nothing
            return
        end if
        call changeSize(sparseMat, sparseMat%nnz + 1)
        call addVal(sparseMat, toIdx, fromIdx, sparseMat%nnz, val)

    end subroutine linkStatesIdx



end module sparseMatrix