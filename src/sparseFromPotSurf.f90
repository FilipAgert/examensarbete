module sparseFromPotSurf
    use potsurf
    use fsparse
    implicit none
    integer, dimension(5,2) :: MIN_MAX_DIM = reshape([MIN_II, MIN_JJ, MIN_KK, MIN_LL, 0, & !MIN_MM, & REMOVE COMMENT IF INCLUDE SYMMETRY
                                                     MAX_II, MAX_JJ, MAX_KK, MAX_LL, MAX_MM], [5,2]) !Array containing min and max of coordinates.

    ! integer, dimension(5,2) :: MIN_MAX_DIM = reshape([MIN_II, MIN_JJ, MIN_KK, MIN_LL, MIN_MM, & !REMOVE COMMENT IF INCLUDE SYMMETRY
    !                                                   MAX_II, MAX_JJ, MAX_KK, MAX_LL, MAX_MM], [5,2]) !Array containing min and max of coordinates.

    integer, dimension(5) :: dimSize = [1 + MAX_II-MIN_II, 1 + MAX_JJ-MIN_JJ,&
                                        1 + MAX_KK-MIN_KK, 1 + MAX_LL-MIN_LL, 1 + MAX_MM]!-MIN_MM] remove comment if including symmetry
    ! integer, dimension(5) :: dimSize = [1 + MAX_II-MIN_II, 1 + MAX_JJ-MIN_JJ,&
    !                                      1 + MAX_KK-MIN_KK, 1 + MAX_LL-MIN_LL, 1 + MAX_MM - MIN_MM]!-MIN_MM] remove comment if including symmetry


contains
    function COOfromPotSurf(AZ, AA, ETOT, II_fusion, fusion_prob, Rneck_fission, fission_prob) result(COO)
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

        

        print* ,' '
        print* ,'Generating sparse matrix of transition probabilities...'
        N = 1
        do i = 1,SIZE(dimSize)
            N = N * dimSize(i) !Number of grid points.
        end do
        INITNNZ = N * 11!! 200_8 !11_8 !! NNZ = N*3⁵    (number of neighborus per coord: 243)
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
                            IDX = linearIdxFromCoordShifted(coord, dimSize) !Convert five dimensional coordinate into one dimensional index
                            neighbours = pruneNeighbours5D(getNeighbours(coord))
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

                                if(MM .eq. 0 .and. neighbours(5,i) .eq. 1) Then 
                                    !If we are at MM = 0 (mass symmetry), then to account for 
                                    !removing all indices MM < 0, make double probability to increase MM (this would then )
                                    prob = prob*2
                                    !print*, "Doubled prob!"
                                endif
                                psum = psum + prob
                                if(prob > 0.0_r_kind) then
                                    
                                    neighbourIDX = linearIdxFromCoordShifted(neighbours(:,i), dimSize)
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
            print*, "Percentage complete: ", 100*II/dimSize(1), "%"
        end do
        
        call COO%malloc(N,N, NNZ)
        COO%data(1:NNZ) = dat(1:NNZ)
        COO%index(:,1:NNZ) = index(:,1:NNZ)
        deallocate(dat, index)
        print*, "Generated sparse matrix."
        print*, "NNZ guess per grid point: ", INITNNZ/N
        print*, "Actual nnz per grid point: ", NNZ/N
        print*, "Number of columns without values: ", CANTEXIT
    end function COOfromPotSurf

    function grid()  !returns matrix with all grid points indexed.
        integer :: N, i, II, JJ, KK ,LL ,MM
        integer, allocatable :: grid(:,:)

        N = 1
        do i = 1,SIZE(dimSize)
            N = N * dimSize(i) !Number of grid points.
        end do
        allocate(grid(5,N))
        i = 0


        do II = MIN_II, MAX_II
            do JJ = MIN_JJ, MAX_JJ
                do KK = MIN_KK, MAX_KK
                    do LL = MIN_LL, MAX_LL
                        do MM = MIN_MAX_DIM(5,1), MAX_MM
                            i = i + 1
                            grid(:,i) = [II, JJ, KK, LL, MM]
                        end do
                    end do
                end do
            end do
        end do


    end function grid

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

    function getNeighboursDiag5D(coord) result(neighbours) !Gets all neighbours including diagonals and self.
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
    end function getNeighboursDiag5D

    function pruneNeighbours5D(neighbours)  result(pruned)
        !Removes all neighbours outside the bounds given by MIN_MAX_DIM
        integer, intent(in), dimension(:, :) :: neighbours
        integer, dimension(:, :), allocatable :: pruned, tempNeighbours     
        
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
    end function pruneNeighbours5D


    function coordFromLinearIdxShifted(idx, dimSize)
    !Function converts linear idx to coordinates
        integer, intent(in) :: idx
        integer, dimension(:), intent(in) :: dimSize
        integer, dimension(SIZE(dimSize)) :: coordFromLinearIdxShifted
        
        integer :: i, tempidx, multfactor
        
        tempidx = idx
        multfactor = 1
        
        do i = 1, size(dimSize) !Calculate multiplication factor
            multfactor = multfactor * dimSize(i)
        end do
        
        do i = SIZE(dimSize), 1, -1
            if (i == 1) then
                coordFromLinearIdxShifted(i) = tempidx
            else
                multfactor = multfactor / dimSize(i)
                coordFromLinearIdxShifted(i) = (tempidx - 1) / multfactor + 1
                tempidx = MOD(tempidx - 1, multfactor) + 1
            endif
        end do

        do i = 1, SIZE(dimSize)
            coordFromLinearIdxShifted(i) = coordFromLinearIdxShifted(i) - 1 + MIN_MAX_DIM(i,1) !This converts coord from range (1,16) to (-5, 10)
        end do
    end function coordFromLinearIdxShifted
    
    function linearIdxFromCoordShifted(coord, dimSize)
    !Function gives linear index from n-dimensional coordinate system.
    !First coordinate varies first.
    !I.e. (x,y,z) with number of spots in each dim: (N_x, N_y, N_z) = (3,3,3)
    !Then (3,1,1) -> Idx: 3
    !     (1,2,1) -> Idx: 4
    !     (1,2,2) -> Idx: 13 
        implicit none
        integer, dimension(:), intent(in) :: coord, dimSize
        integer, dimension(SIZE(coord)) :: coordTransformed
        integer :: i, multFactor
        integer :: linearIdxFromCoordShifted

        do i = 1,SIZE(dimSize)
            coordTransformed(i) = coord(i) - MIN_MAX_DIM(i,1) + 1 !This converts coordinates from range e.g. (-5, 10) to (1, 16)
        end do

        ! Initialize the index
        linearIdxFromCoordShifted = coordTransformed(1)

        ! Handle more than one dimension
        if (SIZE(coord) > 1) then
            multFactor = 1
            do i = 2, SIZE(coord)
                multFactor = multFactor * dimSize(i-1)
                linearIdxFromCoordShifted = linearIdxFromCoordShifted + (coordTransformed(i) - 1) * multFactor
            end do
        end if
    end function linearIdxFromCoordShifted

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
                            fusionIndices(i) = linearIdxFromCoordShifted(coord, dimSize)
                            i = i + 1
                        end do
                    end do
                end do
            end do
        end do
    end function getFusionIndices

    function getFissionIndices(Rneck_fission) result(fissionIndices) !Fusion for II less than or equal II_fusion
        integer, dimension(:), allocatable :: fissionIndices
        integer ::  II, JJ, KK, LL ,MM, NumberOfFissionIndices, i
        integer, dimension(5) :: coord
        real(kind = r_kind), intent(in) :: Rneck_fission

        NumberOfFissionIndices = 0
        do II = MIN_II, MAX_II
            do JJ = MIN_JJ, MAX_JJ
                do KK = MIN_KK, MAX_KK
                    do LL = MIN_LL, MAX_LL
                        do MM = MIN_MAX_DIM(5,1), MAX_MM
                            if(Rneck(II,JJ,KK,LL,MM) < Rneck_fission) then
                                NumberOfFissionIndices = NumberOfFissionIndices + 1
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
                                fissionIndices(i) = linearIdxFromCoordShifted(coord, dimSize)
                            endif
                        end do
                    end do
                end do
            end do
        end do
    end function getFissionIndices
end module sparseFromPotSurf