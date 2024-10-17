module sparse
!This module is a copy of https://github.com/jalvesz/FSPARSE from github
!only change is to change integer(4) -> integer(8) in indexing of data in order to handle
!sparse matrices with nnz > 2^32 / 2 (signed int )= 2e9
    implicit none
    INTEGER(4), PARAMETER :: wpr=8 ! Number of bytes for real FP numbers
    INTEGER(4), PARAMETER :: wpi=8 ! Number of bytes for integer numbers

    private
    type, public, abstract :: sparse_m
      integer :: nrows = 0 !! number of rows
      integer :: ncols = 0 !! number of columns
      integer(kind = wpi) :: nnz   = 0 !! number of non-zero values
    end type

    type, public, extends(sparse_m) :: COO_dp
        integer, allocatable :: index(:,:)
        real(kind=wpr), allocatable :: data(:)
    contains 
        procedure :: malloc => malloc_coo
    end type COO_dp

    type, public, extends(sparse_m) :: CSR_dp
        integer, allocatable :: col(:)
        integer, allocatable :: rowptr(:)
        real(kind = wpr), allocatable :: data(:)
    contains
        procedure :: malloc => malloc_csr
    end type CSR_dp

    public :: coo2ordered, coo2csr

    contains
        subroutine malloc_coo(self, num_rows, num_cols, nnz)
            class(COO_dp) :: self
            integer :: num_rows, num_cols
            integer(kind=wpi) :: nnz
            real(kind=wpr), allocatable :: temp(:)
            integer,  allocatable :: temp_idx(:,:)
            self%nrows = num_rows
            self%ncols = num_cols
            self%nnz = nnz

            if(.not.allocated(self%index)) then
                allocate(temp_idx(2,nnz) , source = 0 )
            else
                allocate(temp_idx(2,nnz) , source = self%index )
            end if
            call move_alloc(from=temp_idx,to=self%index)
            
            if(.not.allocated(self%data)) then
                allocate(temp(nnz)); temp = 0.0_wpr
            else
                allocate(temp(nnz) , source = self%data )
            end if
            call move_alloc(from=temp,to=self%data)
        end subroutine malloc_coo

        subroutine malloc_csr(self,num_rows,num_cols,nnz)
            class(CSR_dp) :: self
            integer, intent(in) :: num_rows
            integer, intent(in) :: num_cols
            integer(kind=wpi), intent(in) :: nnz
            real(kind=wpr), allocatable :: temp(:)    
            integer,  allocatable :: temp_idx(:)
            !-----------------------------------------------------
    
            self%nrows = num_rows
            self%ncols = num_cols
            self%nnz   = nnz
    
            if(.not.allocated(self%col)) then
                allocate(temp_idx(nnz) , source = 0 )
            else
                allocate(temp_idx(nnz) , source = self%col )
            end if
            call move_alloc(from=temp_idx,to=self%col)
    
            if(.not.allocated(self%rowptr)) then
                allocate(temp_idx(num_rows+1) , source = 0 )
            else
                allocate(temp_idx(num_rows+1) , source = self%rowptr )
            end if
            call move_alloc(from=temp_idx,to=self%rowptr)

            
            if(.not.allocated(self%data)) then
                allocate(temp(nnz)); temp = 0.0_wpr
            else
                allocate(temp(nnz) , source = self%data )
            end if
            call move_alloc(from=temp,to=self%data)
        end subroutine malloc_csr

        recursive subroutine quicksort_i(a, first, last)
            ! ref: https://gist.github.com/t-nissie/479f0f16966925fa29ea
            integer, intent(inout) :: a(*)
            integer(wpi), intent(in)    :: first, last
            integer(wpi) :: i, j, x, t

            x = a( (first+last) / 2 )
            i = first
            j = last
            do
                do while (a(i) < x)
                    i=i+1
                end do
                do while (x < a(j))
                    j=j-1
                end do
                if (i >= j) exit
                t = a(i);  a(i) = a(j);  a(j) = t
                i=i+1
                j=j-1
            end do
            if (first < i-1) call quicksort_i(a, first, i-1)
            if (j+1 < last)  call quicksort_i(a, j+1, last)
        end subroutine

        recursive subroutine quicksort_i_dp(a, b, first, last)
            integer, parameter :: wp = wpi
            integer, intent(inout)  :: a(*) !! reference table to sort
            real(kind=wpr), intent(inout) :: b(*) !! secondary real data to sort w.r.t. a(:)
            integer(kind=wpi), intent(in)     :: first, last
            integer(kind=wpi)  :: i, j, x, t
            real(kind=wpr) :: d

            x = a( (first+last) / 2 )
            i = first
            j = last
            do
                do while (a(i) < x)
                    i=i+1
                end do
                do while (x < a(j))
                    j=j-1
                end do
                if (i >= j) exit
                t = a(i);  a(i) = a(j);  a(j) = t
                d = b(i);  b(i) = b(j);  b(j) = d
                i=i+1
                j=j-1
            end do
            if (first < i-1) call quicksort_i_dp(a, b, first, i-1)
            if (j+1 < last)  call quicksort_i_dp(a, b, j+1, last)
        end subroutine 


        subroutine sort_coo_unique_dp( a, data, n, num_rows, num_cols )
            !! Sort a 2d array in increasing order first by index 1 and then by index 2
            integer, parameter :: wp = wpr
            real(wpr), intent(inout) :: data(*)
            integer, intent(inout) :: a(2,*)
            integer(wpi), intent(inout) :: n
            integer, intent(in) :: num_rows
            integer, intent(in) :: num_cols
    
            integer(wpi) :: stride, adr0, adr1, dd
            integer(wpi) :: n_i, pos, ed
            integer, allocatable :: count_i(:), count_i_aux(:), rows_(:), cols_(:)
            real(wpr), allocatable :: temp(:)
            !---------------------------------------------------------
            ! Sort a first time with respect to first index using Count sort
            !How this works is: we count number of values containing each index
            !So we get a list containing number of values in each row
            allocate( count_i( 0:num_rows ) , source = 0 )
            do ed = 1, n
                count_i( a(1,ed) ) = count_i( a(1,ed) ) + 1
            end do
            !Then, we iteratively add each value of the index below it to the index above it.
            !This way, count_i represents pointers to where in data vector the corresponding row starts its data.
            do n_i = 2, num_rows
                count_i(n_i) = count_i(n_i) + count_i(n_i-1)
            end do
            allocate( count_i_aux( 0:num_rows ) , source = count_i )
    
            allocate( rows_(n), cols_(n), temp(n) )
            do ed = n, 1, -1
                n_i = a(1,ed)
                pos = count_i(n_i)
                rows_(pos) = a(1,ed)
                cols_(pos) = a(2,ed)
                temp(pos)  = data(ed)
                count_i(n_i) = count_i(n_i) - 1
            end do
            !---------------------------------------------------------
            ! Sort with respect to second colum using a quicksort
            do n_i = 1, num_rows
                adr0 = count_i_aux(n_i-1)+1
                adr1 = count_i_aux(n_i)
                dd = adr1-adr0+1
                if(dd>0) call quicksort_i_dp(cols_(adr0),temp(adr0),1_wpi,dd)
            end do
            !---------------------------------------------------------
            ! Remove duplicates
            do ed = 1,n
                a(1:2,ed) = [rows_(ed),cols_(ed)]
            end do
            data(1:n) = temp(1:n)
            stride = 0
            do ed = 2, n
                if( a(1,ed) == a(1,ed-1) .and. a(2,ed) == a(2,ed-1) ) then
                    data(ed-1-stride) = data(ed-1-stride) + data(ed) ; data(ed) = data(ed-1-stride)
                    stride = stride + 1
                else
                    a(1:2,ed-stride) = a(1:2,ed)
                    data(ed-stride) = data(ed)
                end if
            end do
            n = n - stride
        end subroutine

        subroutine coo2ordered(COO)
            class(COO_dp), intent(inout) :: COO
            integer, allocatable :: itemp(:,:)
            real(wpr), allocatable :: temp(:)

            call sort_coo_unique_dp(COO%index, COO%data, COO%nnz, COO%nrows, COO%ncols)
            allocate( temp(COO%nnz) , source=COO%data(1:COO%nnz) )
            call move_alloc( temp , COO%data )
            allocate( itemp(2,COO%nnz) , source=COO%index(1:2,1:COO%nnz) )
            call move_alloc( itemp , COO%index )

        end subroutine

        subroutine coo2csr(COO,CSR)
            !! coo2csr: This function enables transfering data from a COO matrix to a CSR matrix
            !! under the hypothesis that the COO is already ordered.
            type(COO_dp), intent(in)    :: COO
            type(CSR_dp), intent(inout) :: CSR
            integer(wpi) :: i
    
            CSR%NNZ = COO%nnz; 
            CSR%nrows = COO%nrows; 
            CSR%ncols = COO%ncols

            if( allocated(CSR%col) ) then
                CSR%col(1:COO%nnz)  = COO%index(2,1:COO%nnz)
                CSR%rowptr(1:COO%nrows) = 0
                CSR%data(1:COO%nnz) = COO%data(1:COO%nnz)
            else 
                allocate( CSR%col(COO%nnz)  , source = COO%index(2,1:COO%nnz) )
                allocate( CSR%rowptr(COO%nrows+1) , source = 0 )
                allocate( CSR%data(COO%nnz) , source = COO%data(1:COO%nnz) )
            end if
    
            CSR%rowptr(1) = 1
            do i = 1, COO%nnz
                CSR%rowptr( COO%index(1,i)+1 ) = CSR%rowptr( COO%index(1,i)+1 ) + 1
            end do
            do i = 1, COO%nrows
                CSR%rowptr( i+1 ) = CSR%rowptr( i+1 ) + CSR%rowptr( i )
            end do
        end subroutine
        
end module sparse