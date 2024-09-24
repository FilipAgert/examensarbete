program testSolver
    use iso_fortran_env, only: dp=>real64
    use sparse_solver_linear_interpolator_module
    use sparse_solver_arnoldi_module
    use dense_solver_module
    use sparseMatrix
    use markovSparse, only: linearIdxFromCoord, gridFromColumnVector
    use denseMatrix, only: gridFromColumnVector, printMatrixToFile, coordFromLinearIdx
    implicit none

    integer :: dimSize(2), fusionCoord(2), fissionCoord(2), fissionIdx(1), fusionIdx(1)
    real(dp), allocatable :: pd(:,:)
    real, allocatable :: grid(:,:)
    integer, allocatable :: startCoord(:,:), startIdxs(:)
    type(sparse_arnoldi_solver) :: solverLin
    type(COO_dp) :: sparseMat
    integer :: i, startNumber
    real, allocatable :: Es(:)
    character(len=8) :: filetext
    startNumber = 15
    allocate(startCoord(2, startNumber))
    allocate(startIdxs(startNumber))
    allocate(Es(startNumber))
    dimSize = [25,25]
    allocate(pd(dimSize(2)*dimSize(1), startNumber))
    allocate(grid(SIZE(dimSize),dimSize(1)*dimSize(2)))
    sparseMat = sparseWalkMatrix(dimSize)
    fusionCoord(1) = dimSize(1) / 2 + 1
    fissionCoord(1) = dimSize(1) / 2 + 1
    fusionCoord(2) = 1
    fissionCoord(2) = dimSize(2)
    
    startCoord(1,:) = dimSize(1)/2 + 1
    do i = 1,startNumber
        Es(i) = real(i)/(real(startNumber + 1))
        startCoord(2,i) = dimSize(2)*i/(startNumber + 1)
        startIdxs(i) = linearIdxFromCoord(startCoord(:,i),dimSize)
    end do

    do i = 1, dimSize(1)*dimSize(2)
        grid(:,i) = coordFromLinearIdx(i, dimSize)
    end do


    fissionIdx(1) = linearIdxFromCoord(fusionCoord,dimSize)
    fusionIdx(1) = linearIdxFromCoord(fissionCoord, dimSize)

    call solverLin%init(sparseMat, grid, startIdxs, Es, fusionIdx, fissionIdx, dimSize)
    call solverLin%solve()
    call solverLin%printResult()
    
    do i = 1,solverLin%result%numResults
        grid = gridFromColumnVector(solverLin%result%getProbabilityDensity(i),dimSize)
        filetext = "file"
        write(filetext, '(I1)')i
        call printMatrixToFile(filetext, real(grid,8))
    end do

    




    

end program testSolver