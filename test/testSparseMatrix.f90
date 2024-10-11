program testSparseMatrix
    use iso_fortran_env, only: dp=>real64
    use denseMatrix, only: coordFromLinearIdx, linearIdxFromCoord
    implicit none

    integer :: dimSize(5), coord(5)
    integer :: index
    index = 60
    dimSize = [53, 15, 15, 15, 41]
    coord = [7,2,1,1,1]
    print*, coordFromLinearIdx(index, dimSize)
    print*, linearIdxFromCoord(coord, dimSize)
    

end program testSparseMatrix