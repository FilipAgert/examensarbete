program testmgmres
    use iso_fortran_env, only: dp=>real64, sp=>real32
    use sparse_solver_linear_interpolator_module
    use sparse_solver_arnoldi_module
    use mgmres
    use dense_solver_module
    use sparseMatrix
    use markovSparse, only: linearIdxFromCoord, gridFromColumnVector
    use denseMatrix, only: gridFromColumnVector, printMatrixToFile, coordFromLinearIdx
    implicit none

    integer :: dimSize(2)
    integer, allocatable :: fusionCoord(:,:), fissionCoord(:,:), fissionIdx(:), fusionIdx(:)
    real, allocatable :: grid(:,:)
    integer, allocatable :: startCoord(:,:), startIdxs(:)
    type(sparse_solver) :: solverLin
    type(sparse_arnoldi_solver) ::solverArnoldi
    type(COO_dp) :: sparseMat
    type(CSR_dp) :: sparseMatCSR
    integer :: i, startNumber, N
    integer, allocatable :: IA(:), JA(:)
    real(dp), allocatable :: A(:), X1(:), X2(:), RHS(:), SOL(:)
    integer :: ITR_MAX, MR, row, col, numchanges, NZ_NUM
    real(dp) :: TOL_ABS, TOL_REL, val
    logical :: verbose
    

    real, allocatable :: Es(:)
    character(len=30) :: filetext
    integer :: fusionPoints, fissionPoints
    
    
    startNumber = 1
    allocate(startCoord(2, startNumber))
    allocate(startIdxs(startNumber))
    allocate(Es(startNumber))
    dimSize = [101,101]
    N = dimSize(1)*dimSize(2)
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
    numChanges = 0
    call solverLin%init(sparseMat, grid, startIdxs, Es, fusionIdx, fissionIdx, dimSize)
    sparseMat =  solverLin%sparseGenerateConnections(startIdxs(1))

    
    do i = 1,sparseMat%nnz !iterate over all nonzero elements, if diagonal element, subtract one from it.
        row = sparseMat%index(1,i)
        col = sparseMat%index(2,i)
        if(row == col) then
            val = sparseMat%data(i)
            sparseMat%data(i) = val - 1
            numChanges = numChanges + 1
        end if
    end do
    if(numChanges /= N) then !
        print*, "WARNING: MATRIX NOT SET UP CORRECTLY, DID NOT SUBTRACT -1 FROM ALL DIAGONAL ELEMENTS"
    end if
    !Then convert to CSR format
    call COO2ordered(sparseMat,.TRUE.)
    call coo2csr(sparseMat, sparseMatCSR)
    !We can deallocate COO matrix.

    allocate(IA(N+1), JA(NZ_NUM), A(NZ_NUM), X1(N), X2(N), SOL(N),RHS(N))
    
    !Call mgmres code
    NZ_NUM = sparseMatCSR%nnz
    IA = sparseMatCSR%rowptr
    JA = sparseMatCSR%col
    A = sparseMatCSR%data
    X1 = 1.0/N
    X2 = 1.0/N
    RHS = 0
    ITR_MAX = 1
    MR = 40
    TOL_ABS = 1
    TOL_REL = 1/(N*1e3)
    verbose = .FALSE.

    RHS(1) = 1
    call pmgmres_ilu_cr(N, NZ_NUM, IA, JA, A, X1, RHS, ITR_MAX, MR, TOL_ABS, TOL_REL,verbose)
    RHS(1) = 0
    RHS(2) = 1
    call pmgmres_ilu_cr(N, NZ_NUM, IA, JA, A, X2, RHS, ITR_MAX, MR, TOL_ABS, TOL_REL,verbose)
    SOL = X2-X1
    SOL = SOL/sum(SOL)
    


    grid = gridFromColumnVector(SOL,dimSize)
    filetext = "file"
    write(filetext, '(A,I0)')"mgmres"
    fileText = trim(fileText)
    call printMatrixToFile(filetext, grid)

    




    

end program testmgmres