!     ******************************************************
!     ** Molecular Dynamic Code usee in MSE 6270 course   **
!     ******************************************************
!     ** Version 1.6 by Leonid V. Zhigilei,               **
!     ** Supports multi-component systems,                **
!     ** periodic, free, and rigid boundary conditions,   **
!     ** Lennard-Jones potential for Ar and a few metals  **
!     ** Stillinger-Weber potential for Si                **
!     ******************************************************
      
!**********************************************************************
!  UNITS:
!    LENGTH      - 1A
!    time        - 1psec
!    MASS        - 1amu = 1.66057D-27 Kg
!    TEMPERATURE - K
!    ENERGY      - 1.0364381D-04 eV = 1.66057D-23 J
!
!**********************************************************************
! All input/output files that you use should be listed in file md.rc **
! UNIT 14 | INPUT DATA (THE MATERIAL AND RUN PARAMETERS)  |READ
! UNIT 15 | INPUT COORDINATES AND VELOCITIES              |READ
! UNIT 16 | OUTPUT COORDINATES AND VELOCITIES             |       WRITE
! UNIT 17 | OUTPUT FILE (ENERGY PER ATOM,TEMP.,time, etc.)|       WRITE
!**********************************************************************
!  Natoms is total number of particles
!  Ntype is the number of particle types
!  Npots=Ntype*(Ntype+1)/2 is number of types of particle pairs
!  For example, for Ntype=3, Npots=6 and for molecules I of type Ktype(I)
!  and J of type Ktype(J) IJindex is defined as
!
!           KTYPE(J) = 1:Ntype
!
!          | 1  2  3
!      K  ----------
!      T  1| 1  3  6
!      Y   |           IJINDEX(KTYPE(I),KTYPE(J)) = 1:Npots
!      P  2| 3  2  5
!      E   |
!     (I) 3| 6  5  4
!
!**********************************************************************
!  HISTORY VARIABLE DEFINES PROPERTIES OF THE PARTICLE OTHER THAN ITS TYPE
!
!     KHIST(J)=1 - full dynamics
!     KHIST(J)=2 - temperature control
!     KHIST(J)=3 - rigid
!     KHIST(J)=4 - analytic constrain - dynamic boundary
!
!**********************************************************************
!  KFLAG DEFINES TEMPERATURE CONTROL
!  KFLAG = 1-Quench, 2-Vel, 3-Heating
!**********************************************************************
!  KEYBS = 0 -pair potential, 1 - Stillinger-Weber

program MD
    use GlobalVars
    implicit none
    integer :: step
    real*8 :: timestart, timeend
    real*8 :: E_p, E_k, E_tot, T, P_tot, P(3)
    character(*), parameter :: file_list = "md.rc"

    call CPU_time(timestart)
    call OpenFiles(file_list)       ! open and read the namelist datafile
    call ReadFiles()                ! reading input files
    
!   Create energies & forces tables
!   Pair potentials and forces are tabulated
    if(KEYBS == 0) then
        NTYPE = 1                   ! Number of particle types
        call EF1LJ(1,1,1)           ! Ar-Ar
    elseif(KEYBS == 1) then
        NTYPE = 2                   ! Number of particle types
        call EF1_SW(1,1,1)          ! Si
        call EF1_SW(2,2,2)          ! Ge
        call EF1_SW(3,1,2)          ! Si-Ge
    elseif(KEYBS == 2) then
        NTYPE = 2                   ! Number of particle types
        call Load_EAM("potentials/Mishin_CuAg.tab")
    elseif(KEYBS == 3) then
        NTYPE = 2                   ! Number of particle types
        call EF1_TF(1,3,3)          ! Cc
!        call EF1_TF(1,1,1)          ! Si
!        call EF1_TF(2,2,2)          ! Ge
!        call EF1_TF(3,1,2)          ! Si-Ge
    else
        write(*,*) "KEYBS = ", KEYBS, " is not defined. Program will stop."
        stop
    endif

    call SetInit()                  ! defining some initial parameters
    call WriteInit()                ! writing initial output before MD

! DEBUG
    time = 10*NWRITE-1
    call WriteSnapshot()
    call WriteRestart()    
    time = 0
    write(*,*) "Wrote inital snapshots and time is reset now to t=", TIME
! END DEBUG

    write(*,*) "NWRITE=", NWRITE
    write(*,*) "NEPRT=", NEPRT


    if(KFLAG == 2) call Vel()

!   Beginning of MD loop
    do step = 0, NSTEP
!       Updating Neighbor List
        if(mod(step,NEWTAB) == 0) call NbList()
        call Nord5()                ! Integration
!       Gathering all particles to the initial cell
        if(mod(step,Ngather) == 0) call Gather()
!       Quenching
        if(KFLAG == 1) call Quench()
!       Slow heating of the material
        if(KFLAG == 3) call Heating()
!       Using the Berendsen thermostat in the stochastic region
        if(KFLAG == 5) call Bere_T()
!       Stretching using custom subrouting by scaling z-direction
        if(KFLAG == 6) call Stretching()
!       Berendsen barostat for constant pressure simulation
        if(LFLAG == 1) call Bere_P()

!       Printing output information
        if(mod(step,NEPRT) == 0) then
            call Temper(E_tot,E_k,E_p,T)
            call Pressure(P_tot, P)
            write(17,"(I6,6(1X,E14.8))") step, time, E_tot, E_k, E_p, T, P_tot
            flush(17)
        endif

!       Writing data for future analysis
        if(mod(step,NWRITE) == 0) then
            call WriteSnapshot()
            call WriteRestart()
            write(*,*) "WROTE TO FILE! time = ", time
        endif

        time = time + dt
    enddo
!   End of MD loop
    
    call WriteEnd(file_list)
    call CPU_time(timeend)
    write(*,*) "The total CPU time = ", timeend - timestart
endprogram
