!   loading interatomic potential for EAM interaction
!
!   The format of *.tab files:
!   number of particle types (N)
!   -----> N blocks for each type
!   mass, atomic_volume, material_id
!   (atomic_volume and material_id are not used in the current code and introduced for backward compatibility of the table format)
!   U: table name, number of points (np), arg_min, arg_max
!   U_data(np): r, U, dU/dr
!   f: table name, number of points (np), arg_min, arg_max
!   f_data(np): r, f, df/dr
!   F: table name, number of points (np), arg_min, arg_max
!   F_data(np): ro, F, dF/dro
!   -----
!   -----> N*(N - 1)/2 blocks for cross interaction potential (t1-t2, t1-t3, ... t1-tN, t2-t3, ...)
!   U: table name, number of points (np), arg_min, arg_max
!   U_data(np): r, U, dU/dr
!   -----

subroutine Load_EAM(filename)
    use unit_Table
    use GlobalVars
    implicit none
    character(*), intent(in) :: filename
    integer :: file_id, ios, i, j, material_id
    real*8 :: atomic_volume

    table_U(:)%energy_units = .true.    !mark if tables should be converted into program units
    table_F(:)%energy_units = .true.
    table_fe(:)%energy_units = .false.
    
    write(*,*) "loading ", trim(filename)
    open(newunit = file_id, file = trim(filename), status = "old", iostat = ios)
    if(ios.ne.0) then
        write(*,*) "cannot open the file ", trim(filename)
        stop
    endif
    read(file_id,*) NTYPE
    if(NTYPE <= 0) then
        write(*,*) "corrupted file: ", trim(filename)
        stop
    elseif(NTYPE > 3) then
        write(*,*) "too many materials specified in ", trim(filename)
        stop
    endif

    write(*,*) "loading tables:"
    do i = 1, NTYPE
        !atomic_volume and material_id are not used in the current code and introduced for backward compatibility of the table format
        read(file_id,*) mass(i), atomic_volume, material_id
        call load_table(file_id, table_U(IJINDEX(i,i)))
        call load_table(file_id, table_fe(i))
        call load_table(file_id, table_F(i))
    enddo
    do i = 1, NTYPE
        do j = i + 1, NTYPE
            call load_table(file_id, table_U(IJINDEX(i,j)))
        enddo
    enddo
    FullList = .false.

    return
endsubroutine
