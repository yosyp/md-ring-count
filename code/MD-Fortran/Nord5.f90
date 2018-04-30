!     Nordsieck Integrator (fifth-order predictor-corrector algorithm)
!     [Nordsieck A. Math.Comput.,v.16,p.22,1962]
!     [Allen & Tildesley, p.340]
!     [Lecture notes for MSE 6270]
!     Kinetic energy is calculated here as well.
!     MSE 6270, Leonid Zhigilei

subroutine NORD5()
    use GlobalVars
    implicit none
    real*8, parameter :: C(6) = (/ 3.d0/20.d0, 251.d0/360.d0, 1.0d0, 11.d0/18.d0, 1.d0/6.d0, 1.d0/60.d0 /)
    real*8 :: correction(1:dim), m
    integer :: i
    
    do i = 1, Natoms
        if(KRIGID(i) == 1) cycle
        XD(1:dim,i) = XD(1:dim,i) + Q1D(1:dim,i) + Q2D(1:dim,i) + Q3D(1:dim,i) + Q4D(1:dim,i) + Q5D(1:dim,i)
        Q1D(1:dim,i) = Q1D(1:dim,i) + 2.0*Q2D(1:dim,i) + 3.0*Q3D(1:dim,i)  + 4.0*Q4D(1:dim,i) + 5.0*Q5D(1:dim,i)
        Q2D(1:dim,i) = Q2D(1:dim,i) + 3.0*Q3D(1:dim,i) + 6.0*Q4D(1:dim,i)  + 10.0*Q5D(1:dim,i)
        Q3D(1:dim,i) = Q3D(1:dim,i) + 4.0*Q4D(1:dim,i) + 10.0*Q5D(1:dim,i)
        Q4D(1:dim,i) = Q4D(1:dim,i) + 5.0*Q5D(1:dim,i)
    enddo

    call Forces()
    
    do i = 1, Natoms
        if(KRIGID(i) == 1) cycle
        m = mass(KTYPE(i))
        correction(1:dim) = FD(1:dim, i)*0.5d0*dt**2 / m - Q2D(1:dim,i)
        Ek(i) = dot_product(Q1D(1:dim,i),Q1D(1:dim,i))*0.5d0*m/dt**2
        XD(1:dim,i) = XD(1:dim,i) + C(1)*correction(1:dim)
        Q1D(1:dim,i) = Q1D(1:dim,i) + C(2)*correction(1:dim)
        Q2D(1:dim,i) = Q2D(1:dim,i) + C(3)*correction(1:dim)
        Q3D(1:dim,i) = Q3D(1:dim,i) + C(4)*correction(1:dim)
        Q4D(1:dim,i) = Q4D(1:dim,i) + C(5)*correction(1:dim)
        Q5D(1:dim,i) = Q5D(1:dim,i) + C(6)*correction(1:dim)
    enddo

    return
endsubroutine
