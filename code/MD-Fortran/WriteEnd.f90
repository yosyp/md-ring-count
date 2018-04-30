!   THIS SUBROUTINE WRITES DATA AT THE END OF THE RUN
!   MSE 6270, Leonid Zhigilei

subroutine WriteEnd(info_file)
    use GlobalVars
    implicit none
    character(*), intent(in) :: info_file
    integer, parameter :: input_id = 10
    integer :: stat, iunit
    character*80 :: info, name, mode
   
!   Write output coordinate file
    call WriteRestart()
    
!   Close the files    
    open(unit = input_id, file = info_file, iostat = stat)
    if(stat == 0) then
        do
            read(input_id, *, iostat = stat) info, iunit, name, mode
            if(is_iostat_end(stat)) exit
            if((info(1:1) == "#").or.(info(1:1) == "*")) cycle
            close(iunit)
        enddo
    endif
    close(input_id)

    return
endsubroutine

















