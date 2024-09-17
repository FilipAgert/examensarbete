program testLinearInterpolation
    use linearInterpolation, only: linear1D

    implicit none
    real :: x0, x1, y0, y1, prevPoint(2),x, y, prevValue(2)
    x0 = 1
    y0 = 2
    x1 = 2
    y1= 4
    prevPoint(1) = x0
    prevValue(1) = y0
    prevPoint(2) = x1
    prevValue(2) = y1
    x = 4
    y = linear1D(x, prevPoint, prevValue)
    print*, "Interpolation: ", y, "Expected result: 8"
end program testLinearInterpolation