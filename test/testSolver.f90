program testSolver
    use iso_fortran_env, only: dp=>real64
    use sparse_solver_module
    use dense_solver_module
    use sparseMatrix
    use markovSparse, only: linearIdxFromCoord, gridFromColumnVector
    use denseMatrix, only: gridFromColumnVector, printMatrixToFile
    implicit none

    integer :: dimSize(2), fusionCoord(2), fissionCoord(2), fissionIdx(1), fusionIdx(1)
    real(dp), allocatable :: pd(:,:)
    integer, allocatable :: startCoord(:,:), startIdxs(:)
    type(sparse_solver) :: sparseSolver
    type(COO_dp) :: sparseMat
    integer :: i, startNumber
    character(len=5) :: integerString
    startNumber = 2
    allocate(startCoord(2, startNumber))
    allocate(startIdxs(startNumber))
    dimSize = [11,11]
    allocate(pd(dimSize(2)*dimSize(1), startNumber))
    sparseMat = sparseWalkMatrix(dimSize)
    fusionCoord(1) = dimSize(1) / 2 + 1
    fissionCoord(1) = dimSize(1) / 2 + 1
    fusionCoord(2) = 1
    fissionCoord(2) = dimSize(2)
    
    startCoord(1,:) = dimSize(1)/2 + 1
    do i = 1,startNumber
        startCoord(2,i) = dimSize(2)*i/(startNumber + 1)
        startIdxs(i) = linearIdxFromCoord(startCoord(:,i),dimSize)
    end do


    fissionIdx(1) = linearIdxFromCoord(fusionCoord,dimSize)
    fusionIdx(1) = linearIdxFromCoord(fissionCoord, dimSize)

    call sparseSolver%init(sparseMat, startIdxs, fusionIdx, fissionIdx, dimSize)
    call sparseSolver%solve()
    call sparseSolver%printResult()





    

end program testSolver