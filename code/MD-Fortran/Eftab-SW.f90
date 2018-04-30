!   Two-body part of Stillinger-Weber potential for Silicon
!   MSE 6270, Leonid Zhigilei and Avinash Dongare
!   [F. H. Stillinger and T. A. Weber, Phys. Rev. B 31, 5262-5271 (1985)]
      
subroutine EF1_SW(KTEF,KE1,KE2)
!   KTEF - type of the potential
!   KE1, KE2 - elements from the list Element(N)
!   Units: XT in Angs, UT in eV, FT in eV/Angs, mass in amu
!   Parameters for Stillinger-Weber potential for Silicon
!     [F. H. Stillinger and T. A. Weber, Phys. Rev. B 31, 5262-5271 (1985)]
!   There are two subsequent papers that give Ge and Si-Ge interactions
!   Paramerets for Stillinger-Weber potential for Ge:
!     [K. Ding and H. Andersen. Phys. Rev. B 34, 6987 (1986)]
!   Paramerets for Stillinger-Weber potential for Si-Ge:
!     [M. Laradji, D.P. Landau and B. Dunweg. Phys. Rev. B 51 4894 (1995)]
!
!   Stillinger and Weber use reduced units.  We use real ones.
    
    use GlobalVars
    implicit none
    integer, intent(in) :: KTEF, KE1, KE2
    real*8 :: A, B, P, Q, ee, rp, rpp, rq, rqq, SIG
    integer :: i
    real*8 :: r, min_r, max_r

    integer, parameter :: Np = 8192
    character*2, parameter :: Element(2) = (/'Si','Ge'/)
    real*8, parameter :: mass_SW(2) = (/28.0855D+00, 72.64D+00/)
    real*8, parameter :: TSW_eps(3) = (/2.167222D+00, 1.93D+00, 2.0427D+00/)
!   SW length units in A
    real*8, parameter :: TSW_sig(3) = (/2.0951D+00, 2.181D+00, 2.13805D+00/)

!   Parameters for 2-body term
!   q and p are integers but we define them
!   as reals to be ready for Si-F potential
    real*8, parameter :: SW_P(3) = (/4.0D+00, 4.0D+00, 4.0D+00/)
    real*8, parameter :: SW_Q(3) = (/0.0D+00, 0.0D+00, 0.0D+00/)
!   A and B in non-reduced units
    real*8, parameter :: SW_A(3) = (/15.2779D+00, 13.6056D+00, 14.40013D+00/)
    real*8, parameter :: SW_B(3) = (/0.60222D+00, 0.60222D+00, 0.60222D+00/)

!   Parameters for 3-body term - Lambda and Gamma in non-reduced units
    real*8, parameter :: TSW_LAM(2) = (/21.0D+00, 31.0D+00/)
    real*8, parameter :: TSW_GAM(3) = (/1.2D+00, 1.2D+00, 1.2D+00/)
    real*8, parameter :: TSW_RM(3) = (/3.7712D+00, 3.9258D+00, 3.8437D+00/)

!   For Si S-W give SW_eps = 3.4723D-19 J = 2.167222 eV
!   Thijsse (following Balamane) proposes to multiply the energy scale
!   by 1.0676 so that the cohesive energy would be 4.63 eV instead of 4.34 eV
    
    write(17,*) "Creating tables for ", Element(KE1),"-", Element(KE2)," Stillinger-Weber potential"

    if(KE1 == KE2) then
        min_r = 0.6*TSW_sig(KE1)
        max_r = TSW_RM(KE1)
        SW_lam(KTEF) = TSW_lam(KE1)
        SW_gam(KTEF) = TSW_gam(KE1)
        SW_sig(KTEF) = TSW_sig(KE1)
        SW_eps(KTEF) = TSW_eps(KE1)
        A = SW_A(KE1)
        P = SW_P(KE1)
        B = SW_B(KE1)
        Q = SW_Q(KE1)
        mass(KTEF) = mass_SW(KE1)
    else
        min_r = 0.6*TSW_sig(3)
        max_r = TSW_RM(3)
        SW_gam(KTEF) = TSW_gam(3)
        SW_sig(KTEF) = TSW_sig(3)
        SW_eps(KTEF) = TSW_eps(3)
        A = SW_A(3)
        P = SW_P(3)
        B = SW_B(3)
        Q = SW_Q(3)
    endif

    table_U(KTEF)%energy_units = .true.    !mark if tables should be converted into program units
    call table_U(KTEF)%init(Np)
    table_U(KTEF)%name = "SW-"//Element(KE1)//"-"//Element(KE2)
    table_U(KTEF)%min_val = min_r
    table_U(KTEF)%max_val = max_r

    SIG = SW_SIG(KTEF)
    do i = 1, Np
        r = min_r + (max_r - min_r)*(i - 1)/(Np - 1)
        ee = exp(SIG/(r - max_r))
        rp = (r/SIG)**(-P)
        rq = (r/SIG)**(-Q)
        rpp = (r/SIG)**(-P-1.0d0)
        rqq = (r/SIG)**(-Q-1.0d0)
        table_U(KTEF)%data(i) = A*(B*rp-rq)*ee
        table_U(KTEF)%d_data(i) = A*(-B*P*rpp+Q*rqq)*ee/SIG &
            - table_U(KTEF)%data(i)*SIG/(r - max_r)**2
    enddo
    table_U(KTEF)%data(Np) = 0.0d0  !set the last point manually, since there is dividing by zero
    table_U(KTEF)%d_data(Np) = 0.0d0
    table_U(KTEF)%data(1:Np) = table_U(KTEF)%data(1:Np)/ENUNIT  ! convert to program units
    table_U(KTEF)%d_data(1:Np) = table_U(KTEF)%d_data(1:Np)/ENUNIT
    FullList = .true.
    
!W  Uncomment this line if you want to plot the potential
!    call print_table(table_U(KTEF))

    return
endsubroutine
