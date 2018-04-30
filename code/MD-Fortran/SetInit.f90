!     THIS SUBROUTINE SETS THE INITIAL PARAMETERS
!     MSE 6270, Leonid Zhigilei

subroutine SetInit()
    use GlobalVars
    implicit none
    integer :: i
    
    Npots = Ntype*(Ntype+1)/2      ! Number of types of particle pairs

!   Defining rigid atoms
    NRIGID = 0
    KRIGID(1:Natoms) = 0
    do i = 1, Natoms
        if((KHIST(I) == 3).and.(KBOUND == 1)) then
            KRIGID(I) = KBOUND
            NRIGID = NRIGID + 1
        endif
    enddo
    if(dim == 2) Q1D(3,1:Natoms) = 0.0d0

    call init_random_seed()
      
    return
endsubroutine

!initialize random seed based on the current time
subroutine init_random_seed()
    integer :: i, n, clock
    integer, allocatable :: seed(:)

    call random_seed(size = n)
    allocate(seed(n))

    call system_clock(count=clock)

    seed = clock + 37*(/(i - 1, i = 1, n)/)
    call random_seed(PUT = seed)
endsubroutine
