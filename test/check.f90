program check
    use iso_fortran_env, only: dp=>real64
    use sparseMatrix
    use denseMatrix, only: printMatrix, walkMatrix
    use fsparse
    implicit none

    type(COO_dp) :: sparse
    integer, dimension(2) :: dimSize, fromCoord, toCoord
    real(dp), dimension(:,:), allocatable :: dense
    dimSize = 3
    fromCoord = 1
    toCoord = 3
    allocate(dense(dimSize(1)*dimSize(2), dimSize(1)*dimSize(2)))
    ! Initialize arrays
    sparse = sparseWalkMatrix(dimSize)
    call linkStates(sparse,fromCoord, toCoord, dimSize)
    call coo2dense(sparse, dense)
end program check

