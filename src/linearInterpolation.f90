module linearInterpolation


    implicit none
    public linear1D
contains
    function linear(coord, prevPoints, prevValues) result (interpolation)
        real, intent(in) :: coord(:), prevPoints(:,:), prevValues(:) !Prevpoints is ordered as first number rotates coordinates, second number rotates point number.
        real :: interpolation
        integer :: dim, numberOfPoints
        dim = SIZE(coord)
        numberOfPoints = size(prevValues)

        






    end function linear
    function linear1D(coord, prevPoints, prevValues) result (interpolation)
        real, intent(in) :: coord, prevPoints(2), prevValues(:,:)
        real :: x1,x2, x
        real, allocatable :: interpolation(:)
        allocate(interpolation(SIZE(prevValues,2)))
        x1 = prevPoints(1)
        x2 = prevPoints(2)
        x = coord
        interpolation = prevValues(1,:)*((x2-x)/(x2-x1)) + prevValues(2,:) *((x-x1)/(x2-x1)) 
    end function linear1D

end module linearInterpolation