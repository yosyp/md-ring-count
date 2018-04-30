!   THIS SUBROUTINE OPENS MOST OF THE INPUT/OUTPUT FILES
!   MSE 6270, Leonid Zhigilei
subroutine OpenFiles(info_file)
    use GlobalVars
    implicit none
    character(*), intent(in) :: info_file
    integer, parameter :: input_id = 10
    integer :: stat, iunit
    character*80 :: info, name, mode

    open(unit = input_id, file = info_file, iostat = stat)
    if(stat.ne.0) then
        write(*,*) "Cannot open ", info_file
        stop
    endif
    do
        read(input_id, *, iostat = stat) info, iunit, name, mode
        if(is_iostat_end(stat)) exit
        if(info(1:1) == "#") then
            cycle
        elseif(info(1:1) == "*") then
            OutDir = trim(name)
            if(len_trim(name) == 0) then
                write(*,*) "Wrong name of a directory ", OutDir
                stop
            endif
        else
            if(len_trim(name) == 0) then
                write(*,*) "Wrong name of a file ", trim(name)
                stop
            endif
            open(unit = iunit, file = trim(name), status = trim(mode), iostat = stat)
            if(stat.ne.0) then
                write(*,*) "Cannot open ", trim(name)
                stop
            endif
        endif
    enddo
    close(input_id)
    
    return
endsubroutine
