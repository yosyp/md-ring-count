!   Implementing Berendsen barostat for keeping constant P
!   At each time step coordinates are scaled by the factor Dzeta.
!   P0 is external pressure in GPa
!   Berendsen et al., J. Chem. Phys., 81 (1984), p. 3684.
!   MSE 6270, Leonid Zhigilei
subroutine BERE_P 
    use GlobalVars
    real*8, parameter :: pbeta = 0.1d-01!pbeta = 0.005d0
    real*8 :: dzeta(3), P_tot, P(3) 
    integer :: i
    
    call Pressure(P_tot, P)
    if(IPCON == 1) then
        dzeta(1:3) = 1.0d0 - pbeta*dt*(P0-P_tot*1.0d-09)
    elseif(IPCON == 2) then
        dzeta(1:2) = 1.0d0 - pbeta*dt*(P0 - (P(1) + P(2))*0.5d-09)
        dzeta(3) = 1.0d0 - pbeta*dt*(P0-P(3)*1.0d-09) 
    else  
        dzeta(1:3) = 1.0d0 - pbeta*dt*(P0-P(1:3)*1.0d-09)
    endif

    do i = 1, 3
        if(boundary(i) == 0) cycle
        CellSize(i) = CellSize(i)*dzeta(i)
        XD(i,1:Natoms) = (XD(i,1:Natoms) - CellCenter(i))*dzeta(i) + CellCenter(i)
    enddo
    
    return
endsubroutine

!   Implementing van Gunstern-Berendsen thermostat for keeping constant T
!   At each time step velocities are scaled by the factor Vscale.
!   Berendsen, J. Chem. Phys., 81 (1984), p. 3684.
!   Although this method does not reproduce canonical ensemble, it is 
!   widely used, and usually gives the same results as more rigorous methods
!   such as GLE. It reproduces the correct average energy,
!   but the distribution is wrong. Therefore, averages are usually correct,
!   but fluctuations are not. 
!   MSE 6270, Leonid Zhigilei
subroutine BERE_T
    use GlobalVars
    implicit none
!   TAU is the time constant that defines the strength of coupling to 
!   the thermal bath.  It is usually chosen from the range of 0.1 - 2 ps.
    real*8, parameter :: tau = 2.0d0
    real*8 :: Vscale, E_p, E_k, E_tot, T
    
    call Temper(E_tot, E_k, E_p, T)
    Vscale = SQRT(1.0d0 + (T0/T - 1.0d0)*dt/tau)  
    Q1D(1:3,1:Natoms) = Q1D(1:3,1:Natoms)*Vscale
    
    return
endsubroutine
