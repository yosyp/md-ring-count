!     Heating of the material in a less violent manner then
!     just by velocities distribution (as in Vel.f [KFLAG=2])
!     MSE 6270, Leonid Zhigilei

subroutine HEATING()
    use GlobalVars
    implicit none
    real*8, parameter :: SKF = 0.001d0
    real*8 :: E_tot, E_k, E_p, T, Vscale
    logical, save :: done = .false.            !keep the value between function calls
    
    if(done) return
    call TEMPER(E_tot, E_k, E_p, T)
    if(T >= T0) then
        done =.true.
        return
    else
        Vscale = sqrt(1.0d0 + SKF*(T0-T)/T0)
        Q1D(1:3,1:Natoms) = Q1D(1:3,1:Natoms)*Vscale
    endif

    return
endsubroutine
      