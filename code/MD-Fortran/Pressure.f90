!     Calculation of the total pressure 
!     MSE 6270, Leonid Zhigilei
subroutine Pressure(P_tot, P)
    use GlobalVars
    implicit none
    real*8, intent(out) :: P_tot, P(3)
    real*8 :: volume, Pv, Pt
    integer :: i, l

    P_tot = 0.0d0
    P(1:3) = 0.0d0
    do i = 1, Natoms
!   Virial part
        Pv = 0.0d0
        do l = 1, dim
            Pv = Pv + STEN(l,l,i)
        enddo
        Pv = 0.5d0*Pv/dim
!   Thermal part
        Pt = Ek(i)*2.0/dim
        P_tot = P_tot + Pv + Pt
        P(1:3) = P(1:3) + Pv + Pt
    enddo
    
    volume = product(CellSize(1:dim))*(1.0d-10)**dim
!   convert units [Pr.En.units/Ang^3] --> [J/m^3] == [Pa]
    P_tot = P_tot*ENUNIT*EVTOJOU/volume
    P = P*ENUNIT*EVTOJOU/volume
    
    return
endsubroutine
