!   Evaluation of forces for the Tersoff  potential for Silicon
!   Avinash M. Dongare and leonid V. Zhigilei [March 2003]
!
!   Copyright (c) 2000-2004 Computational Materials Group, UVa
!   Leonid Zhigilei, http://www.faculty.virginia.edu/CompMat/
!***************************************************************

subroutine F_TF()
    use GlobalParameters
    use GlobalVars
    implicit none

    real*8 :: d_ij(3), d_ik(3), d_jk(3), r_ij, r2_ij, r_ik, r2_ik, r_jk, r2_jk
    real*8 :: f_ij, f_ik, f_jk, dVdBij, E_p, c, dc, tn, beta, sig_ij, tmp, b_ij, db_ij
    real*8 :: costh, dcosdrij, dcosdrik, dcosdrjk, g, dg
    real*8 :: r_cut(KPMX), r_cut2(KPMX)
    integer :: i, j, k, l, ii, ij, ik, jk, kj, kk 

    
    r_cut(1:KPMX) = table_U(1:KPMX)%max_val
    r_cut2(1:KPMX) = (r_cut(1:KPMX))**2
    
    do i = 1, Natoms                            ! loop over ALL particles
        ii = ijindex(KTYPE(i),KTYPE(i))
        tn = TF_N(ii)
        beta = TF_Beta(ii)
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
            
!           Calculates sig_ij in order to calculate bij = (1 + beta^ni*sig_ij^ni)^(-1/2ni)
            sig_ij = 0.0d0
            do kk = 1, NNG(i)                     ! Loop over all neighbours of atom I
                k = NNNG(kk,i)
                if(k == j) cycle
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

!               Use the cosine rule to calculate cos(theta)
                costh = (r2_ij + r2_ik - r2_jk)/(2.0d0*r_ij*r_ik)
                sig_ij = sig_ij + table_C(ik)%get(r_ik)*table_g(ii)%get(costh)
            enddo
            
            tmp = (1.0d0 + (beta*sig_ij)**tn)
            b_ij = tmp**(-0.5d0/tn)
            db_ij = -0.5d0*b_ij*beta*(beta*sig_ij)**(tn - 1.0)/tmp
            Ep(i) = Ep(i) + 0.5d0*table_UA(ij)%get(r_ij)*b_ij
            
!           pair forces
            call table_U(ij)%get_tot(r_ij, E_p, f_ij)
            Ep(i) = Ep(i) + 0.5d0*E_p             ! Potential energy   
            f_ij = -0.5d0*(f_ij + table_UA(ij)%get_d(r_ij)*b_ij)/r_ij  ! Force
            FD(1:3,i) = FD(1:3,i) - f_ij*d_ij(1:3)
            FD(1:3,j) = FD(1:3,j) + f_ij*d_ij(1:3)
!           Static Portion of Stress Tensor
            do l = 1, 3
                STEN(l,l,i) = STEN(l,l,i) + f_ij*d_ij(l)**2
                STEN(l,l,j) = STEN(l,l,j) + f_ij*d_ij(l)**2
            enddo
            
!           3 body forces
            dVdBij = table_UA(ij)%get(r_ij)
            do kk = 1, NNG(i)                     ! Loop over all neighbours of atom I
                k = NNNG(kk,i)
                if(k == j) cycle
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
                
                costh = (r2_ij + r2_ik - r2_jk)/(2.0d0*r_ij*r_ik)
                dcosdrij = 1.0d0/r_ik - costh/r_ij
                dcosdrik = 1.0d0/r_ij - costh/r_ik
                dcosdrjk = -r_jk/(r_ij*r_ik)
                
                call table_g(ii)%get_tot(costh, g, dg)
                call table_C(ik)%get_tot(r_ik, c, dc)
!               Calculate the Forces                
                f_ij = -0.5d0*dVdBij*db_ij*c*dg*dcosdrij/r_ij
                f_ik = -0.5d0*dVdBij*db_ij*(c*dg*dcosdrik + dc*g)/r_ik
                f_jk = -0.5d0*dVdBij*db_ij*c*dg*dcosdrjk/r_jk
                  
                FD(1:3,i) = FD(1:3,i) - f_ij*d_ij(1:3) - f_ik*d_ik(1:3)
                FD(1:3,j) = FD(1:3,j) + f_ij*d_ij(1:3) - f_jk*d_jk(1:3)
                FD(1:3,k) = FD(1:3,k) + f_ik*d_ik(1:3) + f_jk*d_jk(1:3)
                
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
