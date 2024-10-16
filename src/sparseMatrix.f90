module sparseMatrix
    use fsparse
    use potsurf
    implicit none
    

    contains

    function sparseFromPotential(AZ, AA, ETOT, II_fusion, fusion_prob, &
        & Rneck_fission, fission_prob, dimSize, min_max_dim, useFullMMCoordinates) result(COO)
        !! This function generates ALL transition probabilities except from fusion indices to start index.
        type(COO_dp) :: COO
        INTEGER(kind=i_kind), intent(inout) :: AZ, AA
        INTEGER(kind = i_kind), intent(in) :: II_fusion
        real(kind=r_kind), intent(in) :: fusion_prob, fission_prob, Rneck_fission
        REAL(kind=r_kind), intent(inout) :: Etot
        integer :: coord(5)
        integer :: II, JJ, KK, LL, MM, IDX, i, neighbourIDX, NNZ, N, NNZstart, CANTEXIT
        integer(kind=8) :: INITNNZ
        integer, dimension(:,:), allocatable :: neighbours, INDEX
        real(kind=r_kind), dimension(:), allocatable :: dat
        real(kind=r_kind)  :: prob, psum
        integer, dimension(5) :: dimSize
        integer, dimension(5,2) :: MIN_MAX_DIM
        logical :: useFullMMCoordinates

        print* ,' '
        print* ,'Generating sparse matrix of transition probabilities...'
        N = 1
        do i = 1,SIZE(dimSize)
            N = N * dimSize(i) !Number of grid points.
        end do
        INITNNZ = N * 243_8!! 200_8 !11_8 !! NNZ = N*3⁵    (number of neighborus per coord: 243)
        NNZ = 0 !
        print*, "N = ", N
        print*, "NNZ guess: ", INITNNZ

        allocate(INDEX(2,INITNNZ))
        allocate(dat(INITNNZ))
        CANTEXIT = 0
        do II = MIN_II, MAX_II
            do JJ = MIN_JJ, MAX_JJ
                do KK = MIN_KK, MAX_KK
                    do LL = MIN_LL, MAX_LL
                        do MM = MIN_MAX_DIM(5,1), MAX_MM !from 0 to max_mm
                            
                            

                            coord = [II, JJ, KK, LL, MM]
                            IDX = linearIdxFromCoord(coord, dimSize, MIN_MAX_DIM) !Convert five dimensional coordinate into one dimensional index
                            neighbours = pruneNeighbours(getNeighboursDiag(coord),dimSize, MIN_MAX_DIM)
                            psum = 0.0_r_kind
                            NNZstart = NNZ

                            ! if(II .LE. II_fusion) then!IF FUSION: can only stay OR teleport to starting index.
                            !     prob = 1-fusion_prob
                            !     NNZ = NNZ + 1
                            !     index(1,NNZ) = IDX
                            !     index(2,NNZ) = IDX
                            !     dat(NNZ) = prob
                            !     exit
                            ! elseif(Rneck_fission > Rneck(II, JJ, KK, LL, MM)) then !If rneck less than limit, then we have fission
                            !     !State in fission. Can only stay or teleport to starting index.
                            !     prob = 1-fission_prob
                            !     NNZ = NNZ + 1
                            !     index(1,NNZ) = IDX
                            !     index(2,NNZ) = IDX
                            !     dat(NNZ) = prob
                            !     exit
                            !     !Here we only set probability to stay in state. probability to exit state is set later.
                            ! endif
                            ! !We dont go here if fissioned/fusioned.
                            
                            
                            do i = 1,size(neighbours, 2) !Iterate over all neighbours.
                                prob = transition_probability(AZ, AA, ETOT, II, JJ, KK, LL, MM, &
                                neighbours(1,i), neighbours(2,i), neighbours(3,i), neighbours(4,i), neighbours(5,i))


                                if(.NOT. useFullMMCoordinates) then
                                    if(MM .eq. 0 .and. neighbours(5,i) .eq. 1) Then 
                                        !If we are at MM = 0 (mass symmetry), then to account for 
                                        !removing all indices MM < 0, make double probability to increase MM (this would then )
                                        prob = prob*2
                                    endif
                                endif
                                
                                psum = psum + prob
                                if(prob > 0.0_r_kind) then
                                    
                                    neighbourIDX = linearIdxFromCoord(neighbours(:,i),dimSize, MIN_MAX_DIM)
                                    NNZ = NNZ + 1
                                    index(1,NNZ) = neighbourIDX
                                    index(2,NNZ) = IDX
                                    dat(NNZ) = prob
                                endif
                            end do
                            if(NNZstart /= NNZ) then
                                dat(NNZstart + 1:NNZ) = dat(NNZstart + 1 : NNZ) / psum !Normalize to total prob 1.
                                
                                if(II .LE. II_fusion) then !if II less than or equal to, then its fusion
                                    dat(NNZstart + 1:NNZ) = dat(NNZstart + 1:NNZ)*(1-fusion_prob)
                                elseif(Rneck_fission > Rneck(II, JJ, KK, LL, MM)) then !If rneck less than limit, then we have fission
                                    dat(NNZstart + 1:NNZ) = dat(NNZstart + 1:NNZ)*(1-fission_prob)
                                endif
                            else
                                CANTEXIT = CANTEXIT + 1
                            endif
                        end do
                    end do
                end do
            end do
            !print*, "Percentage complete: ", 100*II/dimSize(1), "%"

        end do
        call COO%malloc(N,N, NNZ)
        COO%data(1:NNZ) = dat(1:NNZ)
        COO%index(:,1:NNZ) = index(:,1:NNZ)
        print*, "Generated sparse matrix."
        ! print*, "NNZ guess per grid point: ", real(INITNNZ)/N
        ! print*, "Actual nnz per grid point: ", real(NNZ)/N
        ! print*, "Number of columns without values: ", CANTEXIT
        deallocate(dat, index)
    end function sparseFromPotential

    subroutine connectToStartingCoord(COO, start_c, connectedToStartingIdx, fusionIdxs, fissionIdxs, fissionFusionIndices, &
                                    fusionChance, fissionChance, dimSize, min_max_dim) 
        !modifies existing matrix to connect fusion and fission indices to starting index
        !!TODO: DOES NOT WORK AS INTENDED IF FUSION OR FISSION INDEX NEIGHBOUR TO START_IDX
        type(COO_dp) :: COO
        integer, intent(in) :: start_c(5)
        integer :: i, j, fusionIdx, fissionIdx, startIdx
        integer :: connectionCounter, connectionIndex
        integer :: fissionFusionIndices(:), fissionIdxs(:), fusionIdxs(:)
        logical :: connectedToStartingIdx
        real(kind=r_kind) :: fissionChance, fusionChance
        integer, dimension(5) :: dimSize
        integer, dimension(5,2) :: MIN_MAX_DIM

        startIdx = linearIdxFromCoord(start_c, dimSize, min_max_dim)
        connectionCounter = 0
        connectionIndex = COO%NNZ
        ! print*, "size fission: ", size(fissionIdxs)
        ! print*, "size fusion: ", size(fusionIdxs)
        if(.NOT. connectedToStartingIdx) then  !Generate new connections
            !Extend size of sparsematrix data and index arrays.
            call changeSize(COO, COO%nnz + size(fusionIdxs) + size(fissionIdxs))
            do i = 1, size(fusionIdxs) 
                fusionIdx = fusionIdxs(i)
                !!!!GENERATE CONNECTION TO START_IDX
                connectionCounter = connectionCounter + 1
                connectionIndex = connectionIndex + 1
                fissionFusionIndices(connectionCounter) = connectionIndex
                COO%index(1,connectionIndex) = startIdx
                COO%index(2,connectionIndex) = fusionIdx
                COO%data(connectionIndex) = fusionChance
            end do
            do i = 1,size(fissionIdxs)
                fissionIdx = fissionIdxs(i)
                !!!!GENERATE CONNECTION TO START_IDX
                connectionCounter = connectionCounter + 1
                connectionIndex = connectionIndex + 1
                fissionFusionIndices(connectionCounter) = connectionIndex
                COO%index(1,connectionIndex) = startIdx
                COO%index(2,connectionIndex) = fissionIdx
                COO%data(connectionIndex) = fissionChance
            end do
        else!modify 
            do i = 1,size(fissionFusionIndices)
                COO%index(1, fissionFusionIndices(i)) = startIdx
            end do
        endif
    end subroutine connectToStartingCoord

    function getNeighbours(coord)
        !Gets all neighbours (non diagonal) of a given coordinate. Includes invalid neighbours.
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

    function getNeighboursDiag(coord) result(neighbours) !Gets all neighbours including diagonals and self.
        integer, intent(in), dimension(5) :: coord
        integer, dimension(5, 3**5) :: neighbours !Includes self
        integer :: i,j,k,l,m, startIdx

        startIdx = 1

        do i = coord(1)-1, coord(1) + 1
            do j = coord(2) - 1, coord(2) + 1
                do k = coord(3) - 1, coord(3) + 1
                    do l = coord(4) - 1, coord(4) + 1
                        do m = coord(5) - 1, coord(5) + 1
                            neighbours(:, startIdx) = [i,j,k,l,m]
                            startIdx = startIdx + 1
                        end do
                    end do
                end do
            end do
        end do
    end function getNeighboursDiag

    function pruneNeighbours(neighbours, dimSize, MIN_MAX_DIM)  result(pruned)
        !Removes all neighbours outside the bounds given by MIN_MAX_DIM
        integer, intent(in), dimension(:, :) :: neighbours
        integer, dimension(:, :), allocatable :: pruned, tempNeighbours     
        integer, dimension(5) :: dimSize
        integer, dimension(5,2) :: MIN_MAX_DIM
        integer :: i, j, validCount, numDims, numCols
        logical :: validNeighbour


        
        validCount = 0
        numDims = size(dimSize)
        numCols = size(neighbours, 2)
        allocate(tempNeighbours(numDims, numCols))
        
        do j = 1, numCols !loops over list of all neighbours
            validNeighbour = .TRUE.
            do i = 1, numDims !Checks for each dimension if coordinate out of bounds.
                if (neighbours(i, j) < MIN_MAX_DIM(i,1) .or. neighbours(i,j) > MIN_MAX_DIM(i,2)) then
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
            allocate(pruned(numDims, validCount))
            pruned(:, 1:validCount) = tempNeighbours(:, 1:validCount)
        else
            allocate(pruned(numDims, 0))
        end if
    end function pruneNeighbours

    subroutine matvecpar(matrix, vec_x, vec_y, num_threads)!Modified verson of fsparses matvec to include parallel computing.
        implicit none
        type(CSR_dp), intent(in) :: matrix
        real(r_kind), intent(in)    :: vec_x(:)
        real(r_kind), intent(inout) :: vec_y(:)
        ! integer :: rowptr(matrix%nrows + 1), col(matrix%nnz)
        integer :: i, j, num_threads
        call omp_set_num_threads(num_threads)
        ! dat = matrix%data
        ! col = matrix%col
        ! rowptr = matrix%rowptr
        associate(data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, nnz => matrix%nnz, nrows => matrix%nrows, &
            &ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                !$omp parallel shared(matrix, vec_y, vec_x) private(i,j)
                !$omp do schedule(static)
                do i = 1, nrows
                    do j = matrix%rowptr(i), matrix%rowptr(i+1)-1
                        vec_y(i) = vec_y(i) + matrix%data(j) * vec_x(matrix%col(j))
                    end do
                end do
                !$omp end do
                !$omp end parallel
            endif
        end associate
    end subroutine matvecpar



    function coordFromLinearIdx(idx, dimSize, min_max_dim) result (coord)
    !Function converts linear idx to coordinates I,J,K,L,M
        integer, intent(in) :: idx !Linear index (position in probability denstiy vector.)
        integer, dimension(5) :: dimSize
        integer, dimension(5,2) :: MIN_MAX_DIM
        integer, dimension(SIZE(dimSize)) :: coord
        
        integer :: i, tempidx, multfactor
        
        tempidx = idx
        multfactor = 1
        
        do i = 1, size(dimSize) !Calculate multiplication factor
            multfactor = multfactor * dimSize(i)
        end do
        
        do i = SIZE(dimSize), 1, -1
            if (i == 1) then
                coord(i) = tempidx
            else
                multfactor = multfactor / dimSize(i)
                coord(i) = (tempidx - 1) / multfactor + 1
                tempidx = MOD(tempidx - 1, multfactor) + 1
            endif
        end do

        do i = 1, SIZE(dimSize)
            coord(i) = coord(i) - 1 + MIN_MAX_DIM(i,1) !This converts coord from range (1,16) to (-5, 10)
        end do
    end function coordFromLinearIdx
    
    function linearIdxFromCoord(coord,dimSize, min_max_dim) result(linearIdx)
    !Function gives linear index from n-dimensional coordinate system.
    !First coordinate varies first.
    !I.e. (x,y,z) with number of spots in each dim: (N_x, N_y, N_z) = (3,3,3)
    !Then (3,1,1) -> Idx: 3
    !     (1,2,1) -> Idx: 4
    !     (1,2,2) -> Idx: 13 
        implicit none
        integer, dimension(:), intent(in) :: coord
        integer, dimension(5) :: dimSize
        integer, dimension(5,2) :: MIN_MAX_DIM
        integer, dimension(SIZE(coord)) :: coordTransformed
        integer :: i, multFactor
        integer :: linearIdx

        do i = 1,SIZE(dimSize)
            coordTransformed(i) = coord(i) - MIN_MAX_DIM(i,1) + 1 !This converts coordinates from range e.g. (-5, 10) to (1, 16)
        end do

        ! Initialize the index
        linearIdx = coordTransformed(1)

        ! Handle more than one dimension
        if (SIZE(coord) > 1) then
            multFactor = 1
            do i = 2, SIZE(coord)
                multFactor = multFactor * dimSize(i-1)
                linearIdx = linearIdx + (coordTransformed(i) - 1) * multFactor
            end do
        end if
    end function linearIdxFromCoord


    subroutine changeSize(sparseMat, new_nnz)
        type(COO_dp) :: sparseMat
        integer, intent(in) :: new_nnz
        integer :: current_nnz
        integer, allocatable :: temp_idx(:,:)
        real(r_kind), allocatable :: temp(:)
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
            temp(current_nnz+1:new_nnz) = 0.0_r_kind
        
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
end module sparseMatrix