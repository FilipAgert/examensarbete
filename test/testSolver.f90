program testSolver
    use iso_fortran_env, only: dp=>real64
    use sparse_solver_module
    use denseMatrix
    implicit none

    integer :: dimSize(2), fusionCoord(2), fissionCoord(2), fissionIdx(1), fusionIdx(1)
    real(dp), allocatable :: mat(:,:)
    real(dp), allocatable :: pd(:,:)
    integer, allocatable :: startCoord(:,:), startIdxs(:)
    type(sparse_solver) :: sparseSolver
    integer :: i, startNumber
    character(len=5) :: integerString
    startNumber = 10
    allocate(startCoord(2, startNumber))
    allocate(startIdxs(startNumber))
    dimSize = [101, 101]
    allocate(mat(dimSize(2)*dimSize(1), dimSize(2)*dimSize(1)))
    allocate(pd(dimSize(2)*dimSize(1), startNumber))
    mat = walkMatrix(dimSize)
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

    call sparseSolver%init(mat, startIdxs, fusionIdx, fissionIdx)
    call sparseSolver%solve()
    call sparseSolver%printResult()
    do i = 1, sparseSolver%result%numResults
        write(integerString, '(I0)') i
        pd(:,i) = sparseSolver%result%getProbabilityDensity(i)
        call printMatrixToFile("RUN-"//trim(integerString),gridFromColumnVector(pd(:,i),dimSize))
        integerString = ""
    end do





    

end program testSolver