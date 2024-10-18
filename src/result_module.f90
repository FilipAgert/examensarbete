module result_module
    use iso_fortran_env, only: r_kind=>real64
    implicit none
    private
    public :: convergedResult

    type convergedResult
        real(kind = r_kind), allocatable :: probabilityDensity(:,:)
        real(kind=r_kind), allocatable :: solveTime(:)
        integer, allocatable :: matrixMultiplications(:)
        real(kind=r_kind), allocatable :: fusionFraction(:)
        real(kind=r_kind), allocatable :: fissionMassDistribution(:,:)
        integer, allocatable :: startCoordinate(:,:)
        real(kind=r_kind), allocatable :: energies(:)
        integer :: numResults = 0
    contains 
        procedure, public :: addResult => add_result_method
        procedure, public :: printResult => print_result_method
        procedure, public :: getProbabilityDensity => get_probability_density_method
        procedure, public :: getStartCoord => get_coordinate_method
        procedure, public :: getEnergy => get_energy_method
        procedure, public :: hasResult => has_result_method
        procedure, public :: printResultToFile => print_results_to_file
        procedure, public :: printMassDistribution => print_mass_distribution
    end type convergedResult

contains
    function has_result_method(self)result(res) !Checks if there is available result
        class(convergedResult), intent(in) :: self
        logical :: res
        res = self%numResults > 0
    end function has_result_method

    subroutine add_result_method(self, pd, startCoord, energy, time, multiplications, fusionFrac, fissionMassDistribution)
        class(convergedResult), intent(inout) :: self
        real(kind=r_kind), intent(in) :: pd(:)
        real(kind=r_kind), intent(in) :: time, energy
        integer, intent(in) :: multiplications
        real(kind=r_kind), intent(in) :: fusionFrac
        real(kind=r_kind), dimension(:), intent(in) :: fissionMassDistribution
        integer, intent(in) :: startCoord(:)

        real(kind=r_kind), allocatable :: tempPd(:,:)
        real(kind=r_kind), allocatable :: tempReal(:)
        integer, allocatable :: tempInt(:)
        integer, allocatable :: tempCoord(:,:)


        integer :: n !number of stored results
        n = self%numResults

        if (n == 0) then !allocate arrays
            allocate(self%probabilityDensity(size(pd), 1))
            allocate(self%solveTime(1))
            allocate(self%matrixMultiplications(1))
            allocate(self%energies(1))
            allocate(self%fusionFraction(1))
            allocate(self%fissionMassDistribution(size(fissionMassDistribution),1))
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

            allocate(tempPd(size(fissionMassDistribution), n))
            tempPd = self%fissionMassDistribution
            deallocate(self%fissionMassDistribution)
            allocate(self%fissionMassDistribution(size(fissionMassDistribution), n+1))
            self%fissionMassDistribution(:, 1:n) = tempPd
            deallocate(tempPd)

            allocate(tempReal(n))
            tempReal = self%energies
            deallocate(self%energies)
            allocate(self%energies(n+1))
            self%energies(1:n) = tempReal
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
        self%fissionMassDistribution(:, n) = fissionMassDistribution
        self%startCoordinate(:, n) = startCoord
        self%energies(n) = energy
    end subroutine add_result_method

    subroutine print_result_method(self)
        class(convergedResult), intent(in) :: self
        integer :: i
        integer :: totMults
        real :: totS
        totS = 0
        totMults = 0
    
        if (self%numResults == 0) then
            print *, "No results available"
            return
        end if
    
        ! Print the table headers
        print *, "Results: "
        print *, "---------------------------------------------------------------------------------------------"
        print *, "Index | I | J | K | L | M | Energy [MeV] | Fusion Fraction |  Solve Time [s] | Matrix Mults |"
        print *, "---------------------------------------------------------------------------------------------"
        
        ! Print each result
        do i = 1, self%numResults
            ! Print the index and vertical bar
            write(*, '(I5, 1X, A)', advance="no") i, " |"
    
            !write(*, '(A)', advance = "no") " "
    
            ! Print the remaining fields with vertical bars separating columns
            write(*, '(I3,A, I3,A,I3,A,I3,A,I3,A, F7.1, 6X,A, F14.10, 9X, A,F0.2, 11X,A, I5)', advance="no") &
                self%startCoordinate(1,i), "|",self%startCoordinate(2,i), "|",self%startCoordinate(3,i), "|",&
                self%startCoordinate(4,i),"|",self%startCoordinate(5,i),"|", self%energies(i), "|", self%fusionFraction(i), &
                "|",self%solveTime(i), "|",self%matrixMultiplications(i)
            totS = totS + self%solveTime(i)
            totMults = totMults + self%matrixMultiplications(i)
            print *, "" ! Move to the next line
        end do
        print *, "---------------------------------------------------------------------------------------------"
        print *, "TOTAL TIME [s]: ", totS
        print *, "TOTAL MATRIX MULTIPLICATIONS: ", totMults
        print *, "---------------------------------------------------------------------------------------------"
    end subroutine print_result_method

    subroutine print_mass_distribution(self,idx)
        class(convergedResult), intent(in) :: self
        integer :: idx
        print*, "Fission mass distribution: "
        print*, self%fissionMassDistribution(:,idx)


    end subroutine

    subroutine print_results_to_file(self)
        class(convergedResult), intent(in) :: self
        integer :: i, E
        character(len=100) :: filename
        real(kind=r_kind) :: printPd(size(self%probabilityDensity, 1),1)

        do i = 1,self%numResults
            printPd(:,1) = self%getProbabilityDensity(i)
            WRITE(filename, '(A, F4.1)') "PD-OLD_POT-", self%energies(i)
            print*, filename
            call printMatrixToFile(filename, printPd)
        end do
        
    end subroutine print_results_to_file

    subroutine printMatrixToFile(fileName, A)
        !Prints two dimensional matrix A to file.
        real(r_kind), dimension(:,:), intent(in) :: A
        integer :: i,j
        character (len = *), intent(in) :: fileName
        character (len = 8) :: folderPath
        character (len = 100) :: fullName
        folderPath = '../data/' !Name cannot have prefix /
        folderPath = trim(folderPath) 
        fullName = trim(folderPath//trim(fileName))

        open(unit = 9, file = fullName)
        
        do i = 1,SIZE(A,1)
            do j = 1,SIZE(A,2)
                write(9,'(E12.4)', advance='no') A(i,j)
                
                if (j /= SIZE(A,2)) then !Dont write semicolon after last element
                    write(9,'(A)', advance='no') '; '  ! Add space between elements except last element
                end if
            end do
            write(9,*)
        end do
        
        close(unit = 9)
    end subroutine printMatrixToFile

    function get_probability_density_method(self, index) result(pd)
        class(convergedResult), intent(in) :: self
        real(kind=r_kind), allocatable :: pd(:)
        integer :: index

        if (index < 1 .or. index > self%numResults) then
            print *, "Error: Index out of range."
            error stop
        end if

        pd = self%probabilityDensity(:, index)
    end function get_probability_density_method

    function get_coordinate_method(self, index) result(coord)
        class(convergedResult), intent(in) :: self
        integer ,dimension(:), allocatable :: coord(:)
        integer :: index

        if (index < 1 .or. index > self%numResults) then
            print *, "Error: Index out of range."
            error stop
        end if

        coord = self%startCoordinate(:, index)
    end function get_coordinate_method

    function get_energy_method(self,index) result(energy)
        class(convergedResult), intent(in) :: self
        real(kind=r_kind) :: energy
        integer :: index

        if (index < 1 .or. index > self%numResults) then
            print *, "Error: Index out of range."
            error stop
        end if

        energy = self%energies(index)
    end function get_energy_method

end module result_module