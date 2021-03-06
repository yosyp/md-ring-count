!     Stretching of the material in the z-direction by scaling
!     the axis slowly, similar to HEATING()
subroutine STRETCHING()
    use GlobalVars
    implicit none
    real*8, parameter :: SKF = 1.500d0
    real*8 :: dzeta, ZDC, DZ, myCellCenter
    logical, save :: done = .false.

    ! We want to scale the Z direction by SKF percent over the entire simulation.

    ! ZDI is previous defined when the input .data file is read using:
    !     ZDI = MAXVAL(XD(3,:)) - MINVAL(XD(3,:)) in ReadFiles.f90
    !
    ! DZ is the incremental change in Z at every timestep:
    DZ = (( SKF*ZDI) - ZDI) / (NSTEP/50);

    ! ZDC is the _current_ length of the system
    ZDC = MAXVAL(XD(3,:)) - MINVAL(XD(3,:))

    ! write(*,*) "delta ZDC=", ZDC
    myCellCenter = 0.5*(MAXVAL(XD(3,:)) + MINVAL(XD(3,:)))

    ! dzeta is the amount to multiply all z-coordinates to achieve incremental 
    !   stretching at the current timestep
    dzeta = (ZDC + DZ) / ZDC
    
    XD(3,1:Natoms) = (XD(3,1:Natoms) - myCellCenter)*dzeta + myCellCenter
    ! XD(3,1:Natoms) = (XD(3,1:Natoms) - CellCenter(3))*dzeta + CellCenter(3)
    ! XD(3,1:Natoms) = XD(3,1:Natoms)*dzeta

    return
endsubroutine
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            