!     Neighbour list for all atoms
!     There are several schemes for building the Neighbour list.
!     The one used here is described in Frenkel & Smit pp. 363-367
!     It is NOT memory efficient but easy to understand and better
!     suited for parallel algorithms.
!
!     In the method described in Frenkel & Smit pp. 363-363 they fill
!     all the neighbours of both atoms. For pair potentials, half the
!     list is sufficient.  The array NNNG, however, is less than half
!     full in this case.  For some many-body potentials, such as S-W
!     for Si, and in Monte Carlo calculations, the full neighbour list
!     is needed.  In this code logical variable FullList is used to
!     choose between half list and full list calculation.
!
!     Box filtering is used.  We evaluate dz first and see if it is
!     larger than RList.  Then dx and dy are checked and only then r^2
!     is checked against RList2.  Box filtering is more efficient when
!     the system size is large or one or two dimensions are larger
!     than the others (here we assume that z dimension is the largest.
!
!     Other schemes:
!     Method with a one-dimensional array.  It is more memory efficient
!     but probably less transparent AND, it cannot be made parallel.
!
!   Linked cell method, Frenkel & Smit pp. 368-371, is better/faster
!   for really large systems.
!
!   MSE 6270, Leonid Zhigilei
!**********************************************************************

subroutine NbList()
    use GlobalVars
    implicit none
    real*8 :: d_ij(3), r2, r_cut(KPMX), r_cut2(KPMX)
    integer :: l, i, ij, j, ktypei

    NNG(1:Natoms) = 0  !the number of neighbours of particle i
    r_cut(1:KPMX) = table_U(1:KPMX)%max_val + Rskin
    r_cut2(1:KPMX) = r_cut(1:KPMX)**2
    do i = 1, Natoms
        ktypei = KTYPE(i)
        inner_loop: do j = i + 1, Natoms
            ij = ijindex(KTYPEI, KTYPE(J))
            do l = 1, 3
                d_ij(l) = XD(l,j) - XD(l,i)
!               Periodicity
                if(boundary(l) .ne. 0) then
                    if(d_ij(l) > 0.5d0*CellSize(l))  d_ij(l) = d_ij(l) - CellSize(l)
                    if(d_ij(l) < -0.5d0*CellSize(l)) d_ij(l) = d_ij(l) + CellSize(l)
                endif
                if(abs(d_ij(l)) > r_cut(ij)) cycle inner_loop
            enddo
          
            r2 = dot_product(d_ij,d_ij)
            if(r2 > r_cut2(ij)) cycle
                
!           We've found a neighbour close enough to be on the list
            NNG(i) = NNG(i) + 1     ! increase the count of neighbours
            NNNG(NNG(i),i) = j      ! add the neighbour of atom I

!           Many-body potentials (e.g. SW,EAM) need the complete list
            if(FullList) then
                NNG(j) = NNG(j) + 1
                NNNG(NNG(j),j) = i
            endif
        enddo inner_loop
        
        if(NNG(i) > MAXNNB) then
            write(*,*) "Too many neighbours: ", NNG(I), " > ", MAXNNB
            stop
        endif
    enddo

    return
endsubroutine
