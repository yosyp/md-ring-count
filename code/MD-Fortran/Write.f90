!   Write output files for further analysis
!   MSE 6270, Leonid Zhigilei

subroutine WriteSnapshot()
    use GlobalVars
    implicit none
    integer :: itime, file_id, i
    character*5 :: chtime
    real*8 :: T, E_tot
    
    itime = anint(time)
    if(Itime < 1000) Then
        write(chtime, FMT='(I3.3)') itime
    elseif(Itime < 10000) Then
        write(chtime, FMT='(I4.4)') itime
    elseif(Itime < 100000) Then
        write(chtime, FMT='(I5.5)') itime
    else
        write(*,*) "Simulation is too long... time = ", time, " ps."
        stop
    endif
    
    open(newunit = file_id, FILE='./'//trim(OutDir)//'/time'//trim(chtime)//'.d')
    do i = 1, Natoms
        T = Ek(i)*ENUNIT*2.0d0/(BK*dim)
        E_tot = (Ep(i) + Ek(i))*ENUNIT
        write(file_id, "(I7,1X,I2,3(1X,E12.5),7(1X,E10.3),1X,I2)") &
            i, KHIST(i), XD(1:3,i), Q1D(1:3,i)/dt, Ep(i)*ENUNIT, Ek(i)*ENUNIT, T, E_tot, KTYPE(i)
    enddo
    close(file_id)
    
    return
endsubroutine

subroutine WriteRestart()
    use GlobalVars
    integer :: i
    
    rewind 16
      write(16,*) Natoms, time
      write(16,*) CellSize(1:3), CellCenter(1:3)
      write(16,*) (KTYPE(i),i = 1, Natoms)
      write(16,*) (KHIST(i), XD(1:3,i), i = 1, Natoms)
      write(16,*) (Q1D(1:3,i)/dt, i = 1, Natoms)
    return
endsubroutine

