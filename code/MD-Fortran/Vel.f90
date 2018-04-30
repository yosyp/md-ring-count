!     This subroutine distributes velocities to Natoms particles
!     in order to make even distribution of kinetic energies
!     from 0 to 6kT. Half of this energy will be transferred
!     to the potential energy. (I Assume that the initial
!     configuration is quenched [POT.EN.=0])
!     More about initialization of velocities in Frenkel & Smit p 56-57
!     MSE 6270, Leonid Zhigilei

subroutine Vel()
    use GlobalVars
    implicit none
    real*8 :: r_number(3), p_avr(3), Ek_avr, Ek_tmp, T, correction
    integer :: i

    do i = 1, Natoms
        call random_number(r_number)    !triplet of random numbers
        Q1D(1:3,i) = 2.0d0*dt*sqrt(T0*BK*r_number(1:3)/(ENUNIT*mass(KTYPE(i))))
        
        call random_number(r_number)
        r_number(1:3) = r_number(1:3) - 0.5
        Q1D(1:3,i) = sign(Q1D(1:3,i), r_number(1:3))
    enddo
    
!   average momentum
    p_avr = 0.0d0
    do i = 1, Natoms
        p_avr(1:3) = p_avr(1:3) + mass(KTYPE(i))*Q1D(1:3,i)
    enddo
    if (Natoms > 0) p_avr(1:3) = p_avr(1:3)/Natoms
!   momentum correction
    do i = 1, Natoms
        Q1D(1:3,i) = Q1D(1:3,i) - p_avr(1:3)/mass(KTYPE(i))
    enddo
    
!   kinetic energy and temperature
    Ek_avr = 0.0d0
    do i = 1, Natoms
        Ek_tmp = dot_product(Q1D(1:dim,i),Q1D(1:dim,i))*0.5d0*mass(KTYPE(i))/dt**2
        Ek_avr = Ek_avr + Ek_tmp
    enddo
    if (Natoms > 0) Ek_avr = Ek_avr/Natoms
    T = Ek_avr*ENUNIT/(3.0*BK)
    write(*,"('VEL------> T=',F6.1, ', Tvel=',F6.1, ' K')") T0, T
    correction = sqrt(T0/T)
!   temperature correction
    Q1D(1:3,1:Natoms) = Q1D(1:3,1:Natoms)*correction
    Ek_avr = 0.0d0
    do i = 1, Natoms
        Ek_tmp = dot_product(Q1D(1:dim,i),Q1D(1:dim,i))*0.5d0*mass(KTYPE(i))/dt**2
        Ek_avr = Ek_avr + Ek_tmp
    enddo
    if (Natoms > 0) Ek_avr = Ek_avr/Natoms
    T = Ek_avr*ENUNIT/(3.0*BK)
    write(*,"('VEL------> T=',F6.1, ', Corrected Tvel=',F6.1, ' K')") T0, T
    
    if(dim == 2) Q1D(3,i) = 0.0d0
    do i = 1, Natoms
        if(KRIGID(i) == 1) Q1D(1:3,i) = 0.0d0
    enddo
 
    return
endsubroutine
