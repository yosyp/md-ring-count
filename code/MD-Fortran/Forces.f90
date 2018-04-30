!   Force calculation
!   MSE 6270, Leonid Zhigilei

subroutine Forces()
    use GlobalVars
    implicit none

    FD(1:3,1:Natoms) = 0.0d0
    Ep(1:Natoms) = 0.0d0
    STEN(1:3,1:3,1:Natoms) = 0.0d0

    if(KEYBS == 0) then
        call F_pair()
    elseif (KEYBS == 1) then
        call F_SW()
    elseif (KEYBS == 2) then
        call F_EAM()
    elseif (KEYBS == 3) then
        call F_TF()
    endif
      
    return
endsubroutine





