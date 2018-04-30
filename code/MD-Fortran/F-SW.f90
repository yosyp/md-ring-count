!   Evaluation of forces for the Stillinger Weber potential for Silicon
!   [F. H. Stillinger and T. A. Weber, Phys. Rev. B 31, 5262-5271 (1985)]
!   MSE 6270, Leonid Zhigilei and Avinash Dongare

subroutine F_SW()
    use GlobalVars
    implicit none
    real*8 :: d_ij(3), d_ik(3), d_jk(3), r_ij, r2_ij, r_ik, r2_ik, r_jk, r2_jk
    real*8 :: f_ij, f_ik, f_jk, f_tmp, E_p
    real*8 :: r_cut(KPMX), r_cut2(KPMX)
    real*8 :: costh, dcosdrij, dcosdrik, dcosdrjk, deij, deik, eij, eik, eps, lam
    integer :: i, j, k, l, ij, ik, jk, kj, kk 

    r_cut(1:KPMX) = table_U(1:KPMX)%max_val
    r_cut2(1:KPMX) = (r_cut(1:KPMX))**2
    do i = 1, Natoms                            ! loop over ALL particles
        do kj = 1, NNG(i)                       ! loop over known possible neighbours of i
            j = NNNG(kj,i)                      ! extract the atom number
            d_ij(1:3) = XD(1:3,j) - XD(1:3,i)   ! distances dXij, dYij and dZij
!           Periodicity
            do l = 1, 3
                if(boundary(l) == 0) cycle
                if(d_ij(l) > 0.5d0*CellSize(l))  d_ij(l) = d_ij(l) - CellSize(l)
                if(d_ij(l) < -0.5d0*CellSize(l)) d_ij(l) = d_ij(l) + CellSize(l)
            enddo

            ij = ijindex(KTYPE(i), KTYPE(j))    ! potential type for the pair
            r2_ij = dot_product(d_ij,d_ij)
            if (r2_ij > r_cut2(ij)) cycle
            r_ij = sqrt(r2_ij)                  ! Distance Rij

!           We are using full neighbour list. Therefore, only do pair potential if j > i.
            if(j > i) then
                call table_U(ij)%get_tot(r_ij, E_p, f_ij)
                Ep(i) = Ep(i) + 0.5*E_p         ! Potential energy
                Ep(j) = Ep(j) + 0.5*E_p
                
                f_ij = -f_ij/r_ij               ! Force
                FD(1:3,i) = FD(1:3,i) - f_ij*d_ij(1:3)
                FD(1:3,j) = FD(1:3,j) + f_ij*d_ij(1:3)
!               Static Portion of Stress Tensor
                do l = 1, 3
                    STEN(l,l,i) = STEN(l,l,i) + f_ij*d_ij(l)**2
                    STEN(l,l,j) = STEN(l,l,j) + f_ij*d_ij(l)**2
                enddo
            endif

!           Now we add the three body part due to the neighbours of i
            do kk = kj + 1, NNG(i)             ! Loop over all neighbours of atom I
                k = NNNG(kk,i)
                d_ik(1:3) = XD(1:3,k) - XD(1:3,i) ! distances dXik, dYik and dZik
!               Periodicity
                do l = 1, 3
                    if(boundary(l) == 0) cycle
                    if(d_ik(l) > 0.5d0*CellSize(l))  d_ik(l) = d_ik(l) - CellSize(l)
                    if(d_ik(l) < -0.5d0*CellSize(l)) d_ik(l) = d_ik(l) + CellSize(l)
                enddo
                ik = ijindex(KTYPE(i), KTYPE(k))
                r2_ik = dot_product(d_ik,d_ik)
                
                if(r2_ik > r_cut2(ik)) cycle
                r_ik = sqrt(r2_ik)
!               Calculations of the distance between j and k atoms
                d_jk(1:3) = d_ik(1:3) - d_ij(1:3)
                r2_jk = dot_product(d_jk,d_jk)
                r_jk = sqrt(r2_jk)
                jk = ijindex(KTYPE(j), KTYPE(k))

!               This now forms the triplet of I, J and K atoms
                eij = SW_SIG(ij)/(r_ij - r_cut(ij))
                eik = SW_SIG(ik)/(r_ik - r_cut(ik))
                deij = -SW_SIG(ij)/(r_ij - r_cut(IJ))**2
                deik = -SW_SIG(ik)/(r_ik - r_cut(IK))**2
!               Use the cosine rule to calculate cos(theta)
                costh = (r2_ij + r2_ik - r2_jk)/(2.0d0*r_ij*r_ik)

                dcosdrij = 1.0d0/r_ik - costh/r_ij
                dcosdrik = 1.0d0/r_ij - costh/r_ik
                dcosdrjk = -r_jk/(r_ij*r_ik)

                eps = sqrt(SW_EPS(IJ) * SW_EPS(IK))
                lam = (SW_LAM(KTYPE(j))*((SW_LAM(KTYPE(i)))**2.0d0)*SW_LAM(KTYPE(k)))**0.25d0

                f_tmp = eps*lam*exp(SW_GAM(ij)*eij + SW_GAM(ik)*eik)
                E_p = f_tmp*(costh + 1.0d0/3.0d0)**2
                f_tmp = f_tmp*2.0d0*(costh + 1.0d0/3.0d0)
!               Give all the three body energy to the central atom
                Ep(i) = Ep(i) + E_p/ENUNIT
!               Calculate the Forces
!               Rij derivatives:
                f_ij = -(E_p*SW_GAM(ij)*deij + f_tmp * dcosdrij)/(r_ij*ENUNIT)
                FD(1:3,i) = FD(1:3,i) - f_ij*d_ij(1:3)
                FD(1:3,j) = FD(1:3,j) + f_ij*d_ij(1:3)
!               Rik derivatives:
                f_ik = -(E_p*SW_GAM(ik)*deik + f_tmp * dcosdrik)/(r_ik*ENUNIT)
                FD(1:3,i) = FD(1:3,i) - f_ik*d_ik(1:3)
                FD(1:3,k) = FD(1:3,k) + f_ik*d_ik(1:3)
!               Rjk derivatives:
                f_jk = -f_tmp*dcosdrjk/(r_jk*ENUNIT)
                FD(1:3,j) = FD(1:3,j) - f_jk*d_jk(1:3)
                FD(1:3,k) = FD(1:3,k) + f_jk*d_jk(1:3)
!               Static portion of stress tensor - 3 body
                do l = 1, 3
                    STEN(l,l,i) = STEN(l,l,i) + f_ij*d_ij(l)**2 + f_ik*d_ik(l)**2
                    STEN(l,l,j) = STEN(l,l,j) + f_ij*d_ij(l)**2 + f_jk*d_jk(l)**2
                    STEN(l,l,k) = STEN(l,l,k) + f_ik*d_ik(l)**2 + f_jk*d_jk(l)**2
                enddo
            enddo
        enddo
    enddo
    
    return
endsubroutine
