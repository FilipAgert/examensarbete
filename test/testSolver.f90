program testSolver
    use iso_fortran_env, only: dp=>real64, sp=>real32
    use sparse_solver_linear_interpolator_module
    use sparse_solver_arnoldi_module
    use sparse_solver_arnoldi_shift_module
    use sparse_linear_system_solver_module
    use sparse_linear_system_solver_bicg_module
    use dense_solver_module
    use sparseMatrix
    use markovSparse, only: linearIdxFromCoord, gridFromColumnVector
    use denseMatrix, only: gridFromColumnVector, printMatrixToFile, coordFromLinearIdx
    implicit none

    integer :: dimSize(2)
    integer, allocatable :: fusionCoord(:,:), fissionCoord(:,:), fissionIdx(:), fusionIdx(:)
    real, allocatable :: grid(:,:)
    integer, allocatable :: startCoord(:,:), startIdxs(:)
    type(sparse_arnoldi_shift_solver) :: solverLin
    type(sparse_arnoldi_solver) :: solverArnoldi
    type(COO_dp) :: sparseMat
    integer :: i, startNumber
    real, allocatable :: Es(:)
    character(len=30) :: filetext
    integer :: fusionPoints, fissionPoints
    
    
    startNumber = 10
    allocate(startCoord(2, startNumber))
    allocate(startIdxs(startNumber))
    allocate(Es(startNumber))
    dimSize = [100,100]
    fusionPoints = dimSize(1)/10
    fissionPoints = dimSize(1)/10
    allocate(fusionCoord(2,2*fusionPoints - 1), fissionCoord(2,2*fissionPoints - 1))
    allocate(fusionIdx(2*fusionPoints-1), fissionIdx(2*fissionPoints-1))
    allocate(grid(SIZE(dimSize),dimSize(1)*dimSize(2)))
    sparseMat = sparseWalkMatrix(dimSize)

    fusionCoord(1,1) = dimSize(1) / 2 + 1
    fissionCoord(1,1) = dimSize(1) / 2 + 1
    fusionCoord(2,:) = 1
    fissionCoord(2,:) = dimSize(2)
    startCoord(1,:) = dimSize(1)/2 + 1
    
    do i = 1,(fissionPoints - 1)
        fissionCoord(1,2*i) = dimSize(1)/2 + 1 + i
        fissionCoord(1,2*i+1) = dimSize(1)/2 + 1 - i
    end do

    do i = 1,(fusionPoints - 1)
        fusionCoord(1,2*i) = dimSize(1)/2 + 1 + i
        fusionCoord(1,2*i+1) = dimSize(1)/2 + 1 - i
    end do

    do i = 1,SIZE(fusionIdx)
        fusionIdx(i) = linearIdxFromCoord(fusionCoord(:,i),dimSize)
    end do

    do i = 1,SIZE(fissionIdx)
        fissionIdx(i) = linearIdxFromCoord(fissionCoord(:,i),dimSize)
    end do

    do i = 1,startNumber
        Es(i) = real(i)/(real(startNumber + 1))
        startCoord(2,i) = dimSize(2)*i/(startNumber + 1) + 1
        startIdxs(i) = linearIdxFromCoord(startCoord(:,i),dimSize)
    end do

    do i = 1, dimSize(1)*dimSize(2)
        grid(:,i) = coordFromLinearIdx(i, dimSize)
    end do
    

    call solverLin%init(sparseMat, grid, startIdxs, Es, fusionIdx, fissionIdx, dimSize)
    call solverLin%solve()
    call solverLin%printResult()
    print*, ' '
    print*, 'Arnoldi: '
    call solverArnoldi%init(sparseMat, grid, startIdxs, Es, fusionIdx, fissionIdx, dimSize)
    call solverArnoldi%solve()
    call solverArnoldi%printResult()

    
do i = 1,solverLin%result%numResults
    grid = gridFromColumnVector(solverLin%result%getProbabilityDensity(i),dimSize)
    filetext = "file"
    write(filetext, '(A,I0)')"bicgBIG",i
    fileText = trim(fileText)
    call printMatrixToFile(filetext, grid)
end do

do i = 1,solverArnoldi%result%numResults
    grid = gridFromColumnVector(solverArnoldi%result%getProbabilityDensity(i),dimSize)
    filetext = "file"
    write(filetext, '(A,I0)')"arnoldiBIG-",i
    fileText = trim(fileText)
    call printMatrixToFile(filetext, grid)
end do

    




    

end program testSolver