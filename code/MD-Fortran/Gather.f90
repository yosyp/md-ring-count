!     Gathering molecules to the computational cell
!     MSE 6270, Leonid Zhigilei
!     It is necessary to gather particles back to the
!     computational cell from time to time for systems
!     with active diffusion.  We have to avoid the
!     situation when the particle diffuses farther than
!     the size of the computational cell.
      
subroutine Gather()
    use GlobalVars
    implicit none
    integer :: i, l
    real*8 :: d(3)
      
    do i = 1, Natoms
        d(1:3) = XD(1:3,i) - CellCenter(1:3)
        do l = 1, 3
            if(boundary(l) == 0) cycle
            if(d(l) > 0.5d0*CellSize(l))  XD(l,i) = XD(l,i) - CellSize(l)
            if(d(l) < -0.5d0*CellSize(l)) XD(l,i) = XD(l,i) + CellSize(l)
        enddo
    enddo
      
    return
endsubroutine
