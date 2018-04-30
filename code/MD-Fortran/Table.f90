!A module for function tabulation 
module unit_Table
    use GlobalParameters
    implicit none
    type :: Table                                           !A class representing a table of a function and its derivative
        real*8, public :: min_val, max_val                  !minimum and maximum values of arguments
        integer, public :: n_points                         !number of points in the table
        real*8, public, allocatable :: data(:), d_data(:)   !an array with values of the function and its derivative
        character(len = :), public, allocatable :: name     !a name of the table
        logical, public :: energy_units                     !table should be converted into program units by dividing by ENUNIT
    contains
        procedure :: init                                   !allocate arrays according to a given table size
        procedure :: get                                    !return value of the function for a given arg
        procedure :: get_d                                  !return value of the function derivative for a given arg
        procedure :: get_tot                                !return both the value of function and its derivative for a given arg
    endtype 
    
    public :: load_table                                    !load table from a file
    public :: save_table                                    !save table to a file
    public :: print_table                                   !print table to a file for plotting
    
contains
    subroutine init(this, n)
        class(Table), intent(inout) :: this
        integer, intent(in) :: n
        
        if (n < 1) return
        this%n_points = n
        if(allocated(this%data)) deallocate(this%data)
        if(allocated(this%d_data)) deallocate(this%d_data)
        allocate(this%data(n), this%d_data(n))
    endsubroutine

    real*8 pure function get(this, arg)                     !return value of the function for a given arg
        class(Table), intent(in) :: this
        real*8, intent(in) :: arg
        integer :: idx
        real*8 :: correction
        
        if(arg < this%min_val) then                         !out of bound check, return a value based on the derivative value
            get = this%data(1) + this%d_data(1)*(arg - this%min_val)
            return
        elseif(arg >= this%max_val) then
            get = this%data(this%n_points) + this%d_data(this%n_points)*(arg - this%max_val)
            return
        else                                                !linear interpolation of tabulated values
            correction = (this%n_points - 1)*(arg - this%min_val)/(this%max_val - this%min_val)
            idx = correction
            correction = correction - idx
            idx = idx + 1                                   !array indexes start from 1
            get = this%data(idx)*(1.0d0 - correction) + this%data(idx + 1)*correction
            return
        endif
    endfunction
    
    real*8 pure function get_d(this, arg)                   !return value of the function derivative for a given arg
        class(Table), intent(in) :: this
        real*8, intent(in) :: arg
        integer :: idx
        real*8 :: correction
        
        if(arg < this%min_val) then                         !out of bound check
            get_d = this%d_data(1)
            return
        elseif(arg >= this%max_val) then
            get_d = this%d_data(this%n_points)
            return
        else                                                !linear interpolation of tabulated values
            correction = (this%n_points - 1)*(arg - this%min_val)/(this%max_val - this%min_val)
            idx = correction
            correction = correction - idx
            idx = idx + 1                                   !array indexes start from 1
            get_d = this%d_data(idx)*(1.0d0 - correction) + this%d_data(idx + 1)*correction
            return
        endif
    endfunction
    
    pure subroutine get_tot(this, arg, get_data, get_d_data)!return both the value of function and its derivative for a given arg
        class(Table), intent(in) :: this
        real*8, intent(in) :: arg
        real*8, intent(out) :: get_data, get_d_data
        integer :: idx
        real*8 :: correction
        
        if(arg < this%min_val) then                         !out of bound check, return a value based on the derivative value
            get_data = this%data(1) + this%d_data(1)*(arg - this%min_val)
            get_d_data = this%d_data(1)
            return
        elseif(arg >= this%max_val) then
            get_data = this%data(this%n_points) + this%d_data(this%n_points)*(arg - this%max_val)
            get_d_data = this%d_data(this%n_points)
            return
        else                                                !linear interpolation of tabulated values
            correction = (this%n_points - 1)*(arg - this%min_val)/(this%max_val - this%min_val)
            idx = correction
            correction = correction - idx
            idx = idx + 1                                   !array indexes start from 1
            get_data = this%data(idx)*(1.0d0 - correction) + this%data(idx + 1)*correction
            get_d_data = this%d_data(idx)*(1.0d0 - correction) + this%d_data(idx + 1)*correction
            return
        endif
    endsubroutine
	
    subroutine load_table(file_id, l_table)                 !load table from a file
        integer, intent(in) :: file_id
        type(Table), intent(inout) :: l_table
        character*80 :: table_name
        integer :: i
        real*8 :: arg

        read(file_id,*) table_name, l_table%n_points, l_table%min_val, l_table%max_val
        l_table%name = trim(table_name)
        write(*,"('loading ',A,' table: Np = ',I5,', arg_min = ',F6.2,', arg_max = ',F6.2)") &
            l_table%name, l_table%n_points, l_table%min_val, l_table%max_val
        call l_table%init(l_table%n_points)
        do i = 1, l_table%n_points
            read(file_id,*) arg, l_table%data(i), l_table%d_data(i)
        enddo
        if(l_table%energy_units) then
            l_table%data(:) = l_table%data(:)/ENUNIT
            l_table%d_data(:) = l_table%d_data(:)/ENUNIT
        endif
    endsubroutine
    
    subroutine save_table(file_id, s_table)                 !save table to a file
        integer, intent(in) :: file_id
        type(Table), intent(in) :: s_table
        integer :: i
        real*8 :: arg, unit_conversion
        
        if(s_table%energy_units) then
            unit_conversion = ENUNIT
        else
            unit_conversion = 1.0d0
        endif
        write(*,*) "saving ", s_table%name, " table"
        write(file_id,*) s_table%name, s_table%n_points, s_table%min_val, s_table%max_val
        do i = 1, s_table%n_points
            arg = s_table%min_val + (s_table%max_val - s_table%min_val)*(i - 1)/(s_table%n_points - 1)
            write(file_id,*) arg, s_table%data(i)*unit_conversion, s_table%d_data(i)*unit_conversion
        enddo
    endsubroutine
    
    subroutine print_table(s_table)                         !print table to a file for plotting
        type(Table), intent(in) :: s_table
        integer :: i, file_id
        real*8 :: arg, unit_conversion
        
        if(s_table%energy_units) then
            unit_conversion = ENUNIT
        else
            unit_conversion = 1.0d0
        endif
        write(*,*) "printing ", s_table%name, " table"
        open(newunit = file_id, file = s_table%name//".pot")
        do i = 1, s_table%n_points
            arg = s_table%min_val + (s_table%max_val - s_table%min_val)*(i - 1)/(s_table%n_points - 1)
            write(file_id,*) arg, s_table%data(i)*unit_conversion, s_table%d_data(i)*unit_conversion
        enddo
        close(file_id)
    endsubroutine
endmodule