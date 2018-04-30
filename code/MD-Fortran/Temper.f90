!     Calculation of temperature and averaged energies
!     Energies in eV/particle, Temperature in K.
!     You can modify this to calculate temperature
!     for atoms with KHIST(J)=2 only.
!     MSE 6270, Leonid Zhigilei

subroutine TEMPER(E_tot, E_k, E_p, T)
    use GlobalVars
    implicit none
    real*8, intent(out) :: E_tot, E_k, E_p, T
      
    E_p = sum(Ep(1:Natoms))*ENUNIT
    E_k = sum(Ek(1:Natoms))*ENUNIT
    E_tot = E_p + E_k
      
    if(NRIGID >= Natoms) then
        T = 0.0d0
    else
        T = E_k*2.0d0/(BK*dim*(Natoms - NRIGID))
    endif

    return
endsubroutine
      