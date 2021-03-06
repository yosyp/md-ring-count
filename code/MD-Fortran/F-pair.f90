!     Evaluation of forces for pair potentials with neighbor list
!     Multiple potentials can be used
!     Potential energy is calculated
!     Frenkel & Smit, p.59
!     MSE 6270, Leonid Zhigilei

subroutine F_pair()
    use GlobalVars
    implicit none
    real*8 :: d_ij(3), r2, r, r_cut2(KPMX), E_p, force_ij
    integer :: i, j, k, l, ij
      
    r_cut2(1:KPMX) = (table_U(1:KPMX)%max_val)**2
    do i = 1, Natoms
        do k = 1, NNG(i)
            j = NNNG(k,i)
            d_ij(1:3) = XD(1:3,j) - XD(1:3,i)
          
            do l = 1, 3   !Periodicity
                if(boundary(l) == 0) cycle
                if(d_ij(l) > 0.5d0*CellSize(l))  d_ij(l) = d_ij(l) - CellSize(l)
                if(d_ij(l) < -0.5d0*CellSize(l)) d_ij(l) = d_ij(l) + CellSize(l)
            enddo

            ij = ijindex(KTYPE(i),KTYPE(j))
            r2 = dot_product(d_ij,d_ij)
            if(r2 > r_cut2(ij)) cycle
            r = sqrt(r2)

            call table_U(ij)%get_tot(r, E_p, force_ij)

!           Potential energy (give half the potential to each particle)
            Ep(i) = Ep(i) + 0.5d0*E_p
            Ep(j) = Ep(j) + 0.5d0*E_p
!           force
            force_ij = -force_ij/r
            FD(1:3,i) = FD(1:3,i) - force_ij*d_ij(1:3)
            FD(1:3,j) = FD(1:3,j) + force_ij*d_ij(1:3)
!           Static Portion of Stress Tensor
            do l = 1, 3
                STEN(l,l,i) = STEN(l,l,i) + force_ij*d_ij(l)**2
                STEN(l,l,j) = STEN(l,l,j) + force_ij*d_ij(l)**2
            enddo
        enddo
    enddo
    
    return
endsubroutine









