!   Energy tables for the TERSOFF potential.
!   Parameters obtained from 
!   "J. Tersoff, Modelling solid-state chemistry: Interatomic potentials 
!   for multicomponent systems, Phys. Rev. B 39, 5566 (1989)".
!   Avinash M. Dongare and Leonid V. Zhigilei [July 2003]
!
!   Copyright (c) 2000-2004 Computational Materials Group, UVa
!   Leonid Zhigilei, http://www.faculty.virginia.edu/CompMat/
!***************************************************************

subroutine EF1_TF(KTEF,KE1,KE2)
!   KTEF - type of the potential
!   KE1, KE2 - elements from the list Element(N)   
    use GlobalParameters
    use GlobalVars
    implicit none
    integer, intent(in) :: KTEF, KE1, KE2

    integer, parameter :: Np = 8192
    character*2, parameter :: Element(2) = (/'Si','Ge'/)
    real*8, parameter :: mass_TF(2) = (/28.0855d0, 72.64d0/)
    real*8, parameter :: TSF_A(2) = (/1830.8d0, 1769.0d0/)      !Parameter A
    real*8, parameter :: TSF_B(2) = (/471.18d0, 419.23d0/)      !Parameter B
!   Constants in /A
    real*8, parameter :: TSF_lamda1(2) = (/2.4799d0, 2.4451d0/) !Parameter lamda
    real*8, parameter :: TSF_lamda2(2) = (/1.7322d0, 1.7047d0/) !Parameter mu 
!   Parameters: 
    real*8, parameter :: TSF_KSI(3) = (/1.0d0, 1.0d0, 1.00061d0/) !Parameter Kai   
    real*8, parameter :: TSF_WIJ(3) = (/1.0d0, 1.0d0, 1.0d0/)   !Parameter w
    real*8, parameter :: TSF_Bet(2) = (/1.1d-06, 9.0166d-07/)   !Parameter beta
    real*8, parameter :: TSF_C(2) = (/1.0039d5, 1.0643d5/)      !Parameter c
    real*8, parameter :: TSF_D(2) = (/16.217d0, 15.652d0/)      !Parameter d
    real*8, parameter :: TSF_H(2) = (/-0.59825d0, -0.43884d0/)  !Parameter h
    real*8, parameter :: TSF_N(2) = (/0.78734d0, 0.75627d0/)    !Parameter n
!   Cutoff distance
!   Req + one eighth the distance between second neighbor
    real*8, parameter :: TSF_RC(2) = (/2.7d0, 2.8d0/) 
!   Req + seven eighth the distance between second neighbor
    real*8, parameter :: TSF_SC(2) = (/3.0d0, 3.1d0/)

    integer :: i
    real*8 :: A, B, RC, SC, TF_C, TF_D, TF_H, TF_KSI, TF_WIJ
    real*8 :: r, min_r, max_r, cos_th, f_a, df_a, f_r, df_r, f_c, df_c, lamda1, lamda2
             
    write(17,*) 'Creating tables for ', Element(KE1),'-', Element(KE2),' Tersoff potential' 
    if(KE1 == KE2) then
        min_r = 0.6*TSF_SC(KE1)
        max_r = TSF_SC(KE1)
        TF_C = TSF_C(KE1)
        TF_D = TSF_D(KE1)
        TF_H = TSF_H(KE1)
        TF_N(KTEF) = TSF_N(KE1)
        TF_Beta(KTEF) = TSF_Bet(KE1)
        lamda1 = TSF_lamda1(KE1)
        lamda2 = TSF_lamda2(KE1)
        RC = TSF_RC(KE1) 
        SC = TSF_SC(KE1) 
        A = TSF_A(KE1)
        B = TSF_B(KE1) 
        TF_KSI = TSF_KSI(KE1) 
        TF_WIJ = TSF_WIJ(KE1)
        mass(KTEF) = mass_TF(KE1)
    else
        min_r = 0.6*sqrt(TSF_SC(KE1)*TSF_SC(KE2))
        max_r = sqrt(TSF_SC(KE1)*TSF_SC(KE2))
        lamda1 = 0.5*(TSF_lamda1(KE1) + TSF_lamda1(KE2))
        lamda2 = 0.5*(TSF_lamda2(KE1) + TSF_lamda2(KE2))
        RC = sqrt(TSF_RC(KE1)*TSF_RC(KE2)) 
        SC = sqrt(TSF_SC(KE1)*TSF_SC(KE2)) 
        A = sqrt(TSF_A(KE1)*TSF_A(KE2))
        B = sqrt(TSF_B(KE1)*TSF_B(KE2))
        TF_KSI = TSF_KSI(3)
        TF_WIJ = TSF_WIJ(3)
    endif

    table_U(KTEF)%energy_units = .true.    !mark if tables should be converted into program units
    call table_U(KTEF)%init(Np)
    table_U(KTEF)%name = "TF-"//Element(KE1)//"-"//Element(KE2)
    table_U(KTEF)%min_val = min_r
    table_U(KTEF)%max_val = max_r
    
    table_UA(KTEF)%energy_units = .true.   !mark if tables should be converted into program units
    call table_UA(KTEF)%init(Np)
    table_UA(KTEF)%name = "TF_AC-"//Element(KE1)//"-"//Element(KE2)
    table_UA(KTEF)%min_val = min_r
    table_UA(KTEF)%max_val = max_r
    
    call table_C(KTEF)%init(Np)
    table_C(KTEF)%name = "TF_C-"//Element(KE1)//"-"//Element(KE2)
    table_C(KTEF)%min_val = min_r
    table_C(KTEF)%max_val = max_r
    
    call table_g(KTEF)%init(Np)
    table_g(KTEF)%name = "TF_g-"//Element(KE1)//"-"//Element(KE2)
    table_g(KTEF)%min_val = -1.0
    table_g(KTEF)%max_val = 1.0

    do i = 1, Np
        r = min_r + (max_r - min_r)*(i - 1)/(Np - 1)
        
        f_a = -B*exp(-lamda2*r)
        df_a = -lamda2*f_a
        f_r = A*exp(-lamda1*r)
        df_r = -lamda1*f_r
        f_c = 0.5 + 0.5*cos(PI*(r - RC)/(max_r - RC))
        df_c = -0.5*PI/(max_r - RC)*sin(PI*(r - RC)/(max_r - RC))
        if(r < RC) then
            f_c = 1.0d0
            df_c = 0.0d0
        endif
        
        table_U(KTEF)%data(i) = f_c*f_r
        table_U(KTEF)%d_data(i) = df_c*f_r + f_c*df_r
        
        table_UA(KTEF)%data(i) = TF_KSI*f_c*f_a
        table_UA(KTEF)%d_data(i) = TF_KSI*(df_c*f_a + f_c*df_a)
        
        table_C(KTEF)%data(i) = TF_WIJ*f_c
        table_C(KTEF)%d_data(i) = TF_WIJ*df_c
    enddo
    table_U(KTEF)%data(1:Np) = table_U(KTEF)%data(1:Np)/ENUNIT  ! convert to program units
    table_U(KTEF)%d_data(1:Np) = table_U(KTEF)%d_data(1:Np)/ENUNIT
    table_UA(KTEF)%data(1:Np) = table_UA(KTEF)%data(1:Np)/ENUNIT
    table_UA(KTEF)%d_data(1:Np) = table_UA(KTEF)%d_data(1:Np)/ENUNIT
    
    if(KE1 == KE2) then
        do i = 1, Np
            cos_th = -1.0d0 + 2.0d0*(i - 1)/(Np - 1)

            table_g(KTEF)%data(i) = 1.0d0 + (TF_C/TF_D)**2 - TF_C**2/(TF_D**2 + (cos_th - TF_H)**2)
            table_g(KTEF)%d_data(i) = 2.0d0*TF_C**2*(cos_th - TF_H)/(TF_D**2 + (cos_th - TF_H)**2)**2
        enddo  
    endif
    FullList = .true.
    
!    call print_table(table_U(KTEF))
!    call print_table(table_UA(KTEF))
!    call print_table(table_C(KTEF))
!    call print_table(table_g(KTEF))
    
    return
endsubroutine






