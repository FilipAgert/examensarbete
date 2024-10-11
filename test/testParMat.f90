program testParMat
    use fsparse
    implicit none

    type(CSR_dp) :: SparseMat
    type(COO_dp) :: COO
    real(8), allocatable :: denseMat(:,:), vec(:), tempVec(:), data(:)
    integer, allocatable :: index(:,:)
    real :: time1, time2
    integer :: i, j, matSize, nonzerocount
    matSize = 300*300

    !allocate(denseMat(matSize, matSize))
    allocate(vec(matSize), tempVec(matSize))
    allocate(data(matSize*700), index(2,matSize*700))
    tempVec = 0
    vec = 0
    nonzerocount = 1

    do i = 1,matSize
        do j = 1,matSize
            if (abs(i-j) < 350) then
                data(nonzerocount) = 1 - abs(i-j)/350
                index(1,nonzerocount) = i
                index(2,nonzerocount) = j
                nonzerocount = nonzerocount + 1
            end if
        end do
        vec(i) = 1-1/i
    end do

    call COO%malloc(matSize, matSize, nonzerocount - 1)
    COO%index = index(:,1:nonzerocount - 1)
    COO%data = data(1:nonzerocount - 1)
    !call dense2coo(denseMat, COO)
    call coo2csr(COO, SparseMat)

    call CPU_TIME(time1)
    do i = 1,10
        tempVec = 0
        call matVec(SparseMat, vec, tempVec)
    end do
    call CPU_TIME(time2)
    print*, "2 norm: ", norm2(tempVec)
    print*, "TIME: ", time2 - time1
end program testParMat