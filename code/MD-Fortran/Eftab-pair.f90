!   Making the Energy and Force tables
!   MSE 6270, Leonid Zhigilei
!   Number of pair potentials defined here should be
!   Npots=Ntype*(Ntype+1)/2,  Ntype is number of particle types
!   Lennard - Jones is defined here so far.

subroutine EF1LJ(KTEF, KE1, KE2)
!   KTEF - index of the potential
!   KE1, KE2 - elements from the list Element(N)
!   LJeps is in eV, LJsig is in A, LJmass is in amu
!   [D.V.Matyushov and R. Schmid, J.Chem.Phys. 104, 8627 (1996)] - Ar
!   [T.Halicioglu and G.M.Pound, Phys.Stat.Sol.A 30, 619 (1975)] - metals
!   Cross-species parameters are fixed using the Lorentz-Berthelot rules
!   [see, e.g. J. S. Rowlinson, Liquid and liquid mixtures
!   (Butterworth Scientific, London, 1982]
    
    use GlobalVars
    implicit none
    integer, intent(in) :: KTEF, KE1, KE2
    real*8 :: LJCeps, LJCsig, r, sr6
    real*8 :: min_r, max_r
    integer :: i
    
    integer, parameter :: Np = 8192
    integer, parameter :: Nmaterials = 19
    character*2, parameter :: Element(Nmaterials) = &
    (/'Ar','Al','Ca','Au','Pb','Ni','Pd','Pt','Ag','Cu', &
      'Cr','Fe','Li','Mo','W ','Na','K ','Si','Ge'/)

    real*8, parameter :: LJeps(Nmaterials) = &
    (/0.0103d0, 0.3922d0, 0.2152d0, 0.4415d0, 0.2364d0, 0.5197d0, &
      0.4267d0, 0.6817d0, 0.3448d0, 0.4094d0, 0.5020d0, 0.5264d0, &
      0.2054d0, 0.8379d0, 1.0678d0, 0.1379d0, 0.1144d0, 3.3900d0, &
      2.7400d0/)

    real*8, parameter :: LJsig(Nmaterials) = &
    (/3.405d0,  2.62d0,   3.60d0,   2.637d0,  3.197d0,  2.282d0, &
      2.52d0,   2.542d0,  2.644d0,  2.338d0,  2.336d0,  2.321d0, &
      2.839d0,  2.551d0,  2.562d0,  3.475d0,  4.285d0,  2.0936d0, &
      2.1827d0/)

    real*8, parameter :: mass_LJ(Nmaterials) = &
    (/40.0d0,  27.0d0,   40.1d0,   197.0d0,  207.2d0,  58.7d0, &
      106.4d0,   195.1d0,  107.9d0,  63.6d0,  52.0d0,  55.9d0, &
      6.94d0,  95.9d0,  183.85d0,  23.0d0,  39.1d0,  28.09d0, &
      72.59d0/)

    write(17,*) "Creating tables for ", Element(KE1),"-", Element(KE2)," potential"

    if(KE1 == KE2) then
        LJCeps=LJeps(KE1)
        LJCsig=LJsig(KE1)
        mass(IJINDEX(KTEF, KTEF)) = mass_LJ(KE1)   ! Mass of the particle in pr.unit [aem]
    else
        LJCeps=SQRT(LJeps(KE1)*LJeps(KE2))
        LJCsig=(LJsig(KE1)+LJsig(KE2))/2.0d0
    endif
    min_r = 0.6*LJCsig
    max_r = 2.5*LJCsig                  ! Cutoff distance for the interaction potential [A]

    table_U(KTEF)%energy_units = .true. !mark if tables should be converted into program units
    call table_U(KTEF)%init(Np)
    table_U(KTEF)%name = "LJ-"//Element(KE1)//"-"//Element(KE2)
    table_U(KTEF)%min_val = min_r
    table_U(KTEF)%max_val = max_r
    
    do i = 1, Np
        r = min_r + (max_r - min_r)*(i - 1)/(Np - 1)
        sr6 = (LJCsig/r)**6                                                 ! (sigma/r)**6
        table_U(KTEF)%data(i) = 4.0d0 * LJCeps * sr6 * (sr6 - 1.0d0)        ! potential
        table_U(KTEF)%d_data(i) = -LJCeps * 48.0d0/r * sr6 * (sr6 - 0.5d0)  ! potential derivative
    enddo
    table_U(KTEF)%data(1:Np) = table_U(KTEF)%data(1:Np)/ENUNIT              ! convert to program units
    table_U(KTEF)%d_data(1:Np) = table_U(KTEF)%d_data(1:Np)/ENUNIT
    FullList = .false.
    
!W  Uncomment this line if you want to plot the potential
!    call print_table(table_U(KTEF))
      
    return
endsubroutine











