!     Stretching of the material in the z-direction by scaling
!     the axis slowly, similar to HEATING()
!     Yosyp Schwab 

subroutine STRETCHING()
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
      