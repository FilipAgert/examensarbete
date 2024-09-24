module fileUtils
    use iso_fortran_env, only: dp=>real64
    implicit none
    contains
    subroutine printMatrixToFile(fileName, A)
        !Prints two dimensional matrix A to file.
        real(dp), dimension(:,:), intent(in) :: A
        integer :: i,j
        character (len = *), intent(in) :: fileName
        character (len = 5) :: folderPath
        folderPath = 'data/' !Name cannot have prefix /
        folderPath = trim(folderPath) 

        open(unit = 9, file = trim(folderPath//trim(fileName)))
        
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

    subroutine readMatrixFromFile(fileName, matrix)
        ! Reads a matrix from a file where values are separated by semicolons
        character(len=*), intent(in) :: fileName
        real(dp), allocatable, dimension(:,:) :: matrix
        character(len=1000) :: line
        character(len=5) :: folderPath
        character(len=20), dimension(:), allocatable :: elements
        integer :: i, j, rowCount, colCount
        real(dp) :: value
    
        folderPath = 'data/'
        folderPath = trim(folderPath)
    
        ! Open the file for reading
        open(unit=9, file=trim(folderPath//trim(fileName)))
    
        ! First, count rows and columns to allocate the matrix
        rowCount = 0
        colCount = 0
    
        ! Loop over the file to count rows and columns
        do
            read(9,'(A)', iostat=i) line
            if (i /= 0) exit  ! Exit loop if end of file or error
    
            ! Count columns (number of semicolons + 1)
            call splitLine(line, elements)
            if (rowCount == 0) then
                colCount = size(elements)
            end if
            rowCount = rowCount + 1
        end do
    
        ! Allocate the matrix based on rowCount and colCount
        allocate(matrix(rowCount, colCount))
    
        ! Rewind the file to start reading values
        rewind(9)
    
        ! Read the matrix values from the file
        i = 1
        do
            read(9,'(A)', iostat=j) line
            if (j /= 0) exit  ! Exit loop if end of file or error
    
            call splitLine(line, elements)
            do j = 1, size(elements)
                read(elements(j), *) value
                matrix(i,j) = value
            end do
            i = i + 1
        end do
    
        close(9)
    end subroutine readMatrixFromFile

    subroutine splitLine(line, elements)
        ! Splits a line of text based on semicolons
        character(len=*), intent(in) :: line
        character(len=20), dimension(:), allocatable, intent(out) :: elements
        integer :: i, startPos, endPos, count
    
        allocate(elements(1))
        count = 0
        startPos = 1
    
        do i = 1, len_trim(line)
            if (line(i:i) == ';') then
                endPos = i - 1
                count = count + 1
                if (count > size(elements)) then
                    allocate(elements(count))
                end if
                elements(count) = adjustl(line(startPos:endPos))
                startPos = i + 2  ! Skip semicolon and space
            end if
        end do
    
        ! Capture the last element after the last semicolon
        if (startPos <= len_trim(line)) then
            count = count + 1
            allocate(elements(count))
            elements(count) = adjustl(line(startPos:len_trim(line)))
        end if
    end subroutine splitLine


end module fileUtils