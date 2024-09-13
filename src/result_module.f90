module result_module
    use iso_fortran_env, only: dp=>real64
    implicit none
    private
    public :: convergedResult

    type convergedResult
        real(dp), allocatable :: probabilityDensity(:,:)
        real, allocatable :: solveTime(:)
        integer, allocatable :: matrixMultiplications(:)
        real, allocatable :: fusionFraction(:)
        real, allocatable :: fissionFraction(:)
        real, allocatable :: startCoordinate(:,:)
        integer :: numResults = 0
    contains 
        procedure, public :: addResult => add_result_method
        procedure, public :: printResult => print_result_method
        procedure, public :: getProbabilityDensity => get_probability_density_method
        procedure, public :: getStartCoord => get_coordinate_method
        procedure, public :: hasResult => has_result_method
    end type convergedResult

contains
    function has_result_method(self)result(res) !Checks if there is available result
        class(convergedResult), intent(in) :: self
        logical :: res
        res = self%numResults > 0
    end function has_result_method

    subroutine add_result_method(self, pd, startCoord, time, multiplications, fusionFrac, fissionFrac)
        class(convergedResult), intent(inout) :: self
        real(dp), intent(in) :: pd(:)
        real, intent(in) :: time
        integer, intent(in) :: multiplications
        real, intent(in) :: fusionFrac
        real, intent(in) :: fissionFrac
        real, intent(in) :: startCoord(:)

        real(dp), allocatable :: tempPd(:,:)
        real, allocatable :: tempReal(:)
        integer, allocatable :: tempInt(:)
        real, allocatable :: tempCoord(:,:)


        integer :: n !number of stored results
        n = self%numResults

        if (n == 0) then !allocate arrays
            allocate(self%probabilityDensity(size(pd), 1))
            allocate(self%solveTime(1))
            allocate(self%matrixMultiplications(1))
            allocate(self%fusionFraction(1))
            allocate(self%fissionFraction(1))
            allocate(self%startCoordinate(size(startCoord), 1))
            
        else 
            ! Resize probabilityDensity with a temporary array
            allocate(tempPd(size(pd), n))
            tempPd = self%probabilityDensity
            deallocate(self%probabilityDensity)
            allocate(self%probabilityDensity(size(pd), n+1))
            self%probabilityDensity(:, 1:n) = tempPd
            deallocate(tempPd)

            ! Resize solveTime with a temporary array
            allocate(tempReal(n))
            tempReal = self%solveTime
            deallocate(self%solveTime)
            allocate(self%solveTime(n+1))
            self%solveTime(1:n) = tempReal
            deallocate(tempReal)

            ! Resize matrixMultiplications with a temporary array
            allocate(tempInt(n))
            tempInt = self%matrixMultiplications
            deallocate(self%matrixMultiplications)
            allocate(self%matrixMultiplications(n+1))
            self%matrixMultiplications(1:n) = tempInt
            deallocate(tempInt)

            ! Resize fusionFraction with a temporary array
            allocate(tempReal(n))
            tempReal = self%fusionFraction
            deallocate(self%fusionFraction)
            allocate(self%fusionFraction(n+1))
            self%fusionFraction(1:n) = tempReal
            deallocate(tempReal)

            ! Resize fissionFraction with a temporary array
            allocate(tempReal(n))
            tempReal = self%fissionFraction
            deallocate(self%fissionFraction)
            allocate(self%fissionFraction(n+1))
            self%fissionFraction(1:n) = tempReal
            deallocate(tempReal)

            ! Resize startCoordinate with a temporary array
            allocate(tempCoord(size(startCoord), n))
            tempCoord = self%startCoordinate
            deallocate(self%startCoordinate)
            allocate(self%startCoordinate(size(startCoord), n+1))
            self%startCoordinate(:, 1:n) = tempCoord
            deallocate(tempCoord)
        end if

        ! Add the new result data
        n = n + 1
        self%numResults = n

        self%probabilityDensity(:, n) = pd
        self%solveTime(n) = time
        self%matrixMultiplications(n) = multiplications
        self%fusionFraction(n) = fusionFrac
        self%fissionFraction(n) = fissionFrac
        self%startCoordinate(:, n) = startCoord
    end subroutine add_result_method

    subroutine print_result_method(self)
        class(convergedResult), intent(in) :: self
        integer :: i, j
        integer :: coordDim
        integer :: coordWidth
    
        if (self%numResults == 0) then
            print *, "No results available"
            return
        end if
    
        ! Determine the dimension size for startCoordinate
        coordDim = size(self%startCoordinate, 1)
        coordWidth = 2  ! Adjust this based on expected coordinate width
    
        ! Print the table headers
        print *, "Results: "
        print *, "--------------------------------------------------------------"
        print *, "Index | ", &
            'Coords (', coordDim, 'D)', &
            " | Fusion Fraction | Fission Fraction | Solve Time [s] | Matrix Mults"
        print *, "--------------------------------------------------------------"
        
        ! Print each result
        do i = 1, self%numResults
            ! Print the index and vertical bar
            write(*, '(I5, 1X, A)', advance="no") i, " |"
    
            ! Print the coordinates with proper alignment
            write(*, '(A)', advance="no") " "
            do j = 1, coordDim
                if (j > 1) then
                    write(*, '(A)', advance="no") " |"
                end if
                write(*, '(F5.2, A)', advance="no") self%startCoordinate(j, i), " "
            end do
            write(*, '(A)', advance = "no") " "
    
            ! Print the remaining fields with vertical bars separating columns
            write(*, '(A, F4.2, 15X, F4.2, 15X, F5.2, 12X, I5)', advance="no") &
                " |  ", self%fusionFraction(i), self%fissionFraction(i), &
                self%solveTime(i), self%matrixMultiplications(i)
            
            print *, "" ! Move to the next line
        end do
    
        print *, "--------------------------------------------------------------"
    end subroutine print_result_method

    function get_probability_density_method(self, index) result(pd)
        class(convergedResult), intent(in) :: self
        real(dp), allocatable :: pd(:)
        integer :: index

        if (index < 1 .or. index > self%numResults) then
            print *, "Error: Index out of range."
            error stop
        end if

        pd = self%probabilityDensity(:, index)
    end function get_probability_density_method

    function get_coordinate_method(self, index) result(coord)
        class(convergedResult), intent(in) :: self
        real, dimension(:), allocatable :: coord(:)
        integer :: index

        if (index < 1 .or. index > self%numResults) then
            print *, "Error: Index out of range."
            error stop
        end if

        coord = self%startCoordinate(:, index)
    end function get_coordinate_method

end module result_module