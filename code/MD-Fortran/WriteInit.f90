!    HIS SUBROUTINE WRITES INITIAL INFORMATIOIN FOR THE RUN
!   MSE 6270, Leonid Zhigilei

subroutine WriteInit()
    use GlobalVars
    implicit none
    integer :: KTT(NTYPE), i
    real*8 :: DENSM, DENSN, M_tot, volume

!   Counting particles of different types
    KTT(1:NTYPE) = 0
    do i = 1, Natoms
        if((KTYPE(i) > 0).and.(KTYPE(i) <= NTYPE)) KTT(KTYPE(i)) = KTT(KTYPE(i)) + 1
    enddo
    if(sum(KTT).ne.Natoms) then
        write(*,*)  "KTYPE(I), NTYPE, and Natoms do not match, program stops"
        stop
    endif
    M_tot = dot_product(KTT(1:NTYPE), mass(1:NTYPE))    ! Total mass in amu
    volume = product(CellSize(1:dim))*(1.0d-8)**dim     ! Volume in cm3
    DENSM = M_tot*AMUTOKG*1.0D+03/volume                ! gr/cm3
    DENSN = Natoms/volume                               ! particles/cm3

    write(17,*) "   -----> MD RUN:"
    write(17,"('NSTEP = ',I7,' (Number of steps)')") NSTEP
    write(17,"('dt = ',F7.5,' (Integration time step in [psec])')") dt
    write(17,"('KFLAG = ',I1,' (1-Quench, 2-Velocity distribution, 5-Constant Temperature)')") KFLAG
    write(17,"('LFLAG = ',I1,' (1-Constant Pressure)')") LFLAG
    write(17,"('KEYBS = ',I1,' (0-pair potential, 1-Stillinger Weber, 2-EAM, 3-Tersoff)')") KEYBS
    write(17,"('LIDX = ',I1,' (1-periodicity in X direction,0-free boundaries)')") Boundary(1)
    write(17,"('LIDY = ',I1,' (1-periodicity in Y direction,0-free boundaries)')") Boundary(2)
    write(17,"('LIDZ = ',I1,' (1-periodicity in Z direction,0-free boundaries)')") Boundary(3)
    write(17,"('KBOUND = ',I2,' (0-Free,1-rigid)')") KBOUND
    write(17,"('dim = ',I2,' (3 - 3D simulation,2- 2D simulation)')") dim
    write(17,"('T0 = ',F6.1,' (Temperature distributed by VEL,[K])')") T0
    write(17,"('P0 = ',F6.1,' (Desired pressure in con. pressure control,[GPa])')") P0
    write(17,"('NEWTAB = ',I3,' (Step of neighbors list renewal)')") NEWTAB
    write(17,"('NEPRT = ',I4,' (Step of printing output information)')") NEPRT
    write(17,"('NWRITE = ',I5,' (Step of writing output information)')") NWRITE
    write(17,"('Ngather = ',I4,' (Step of gathering molecules to the comp. cell)')") Ngather
    write(17,"('RSkin = ',F6.3,' (Skin depth in neighbor list calculation [A])')") RSkin
    write(17,*)

    write(17,*) "   -----> MATERIAL:"
    write(17,"('XL = ',D11.5,' (X size of the computational cell [A])')") CellSize(1)
    write(17,"('YL = ',D11.5,' (Y size of the computational cell [A])')") CellSize(2)
    write(17,"('ZL = ',D11.5,' (Z size of the computational cell [A])')") CellSize(3)
    if(dim == 2) then
        write(17,"('DENSM = ',D11.5,' (Density of the material [gm/cm2])')") DENSM
        write(17,"('DENSN = ',D11.5,' (Density of the material [molec./cm2])')") DENSN
    else
        write(17,"('DENSM = ',D11.5,' (Density of the material [gm/cm3])')") DENSM
        write(17,"('DENSN = ',D11.5,' (Density of the material [molec./cm3])')") DENSN
    endif
    write(17,"('Natoms = ',I6,' (Number of particles in the computational cell)')") Natoms
    write(17,"('NTYPE = ',I2,' (Number of particle types)')") NTYPE
    write(17,*)
    
    write(17,*) "   -----> TYPES OF PARTICLES:"
    write(17,"('Type ',I2,': ',I5,' particles.  mass = ',D11.5,' [amu].')") (i, KTT(i), mass(i), i = 1, NTYPE)
    write(17,*)
    write(17,*) " Step     time     Energy      Kinetic    Potential  Temperature    Pressure"
    
    return
endsubroutine




