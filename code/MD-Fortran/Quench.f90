!     This is the fastest way to stop motions in the system and
!     to quench it down to the minimum (local or global) of the
!     potential surface
!     MSE 6270, Leonid Zhigilei

subroutine Quench()
    use GlobalVars
    implicit none
    integer :: i
!   We stop atom I when its kinetic energy Ek(I) starts to
!   decrease. Ek_prev(I) is the kinetic energy of the atom I from
!   the previous step.

    do i = 1, Natoms
        if(Ek(i) > Ek_prev(i)) then
            Ek_prev(i) = Ek(i)
        else
            Ek_prev(i) = 0.0d0
            Q1D(1:3,i) = 0.0d0
            Ek(i) = 0.0d0
        endif
    enddo

    return
endsubroutine





