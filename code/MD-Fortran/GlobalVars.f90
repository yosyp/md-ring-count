!A module with global variables
module GlobalVars
    use GlobalParameters
    use unit_Table
    implicit none
    
    real*8, allocatable :: XD(:,:), FD(:,:)         !Coordinates and forces,
    real*8, allocatable :: Q1D(:,:)                 !Velocities*dt
    real*8, allocatable :: Q2D(:,:), Q3D(:,:)       !Higher time derivatives of the coordinates
    real*8, allocatable :: Q4D(:,:), Q5D(:,:)       !for Nordsieck integrator
    integer, allocatable :: KTYPE(:), KHIST(:)      !type and khist of atoms
    integer, allocatable :: KRIGID(:)               !rigid flags
    real*8, allocatable :: Ep(:), Ek(:), Ek_prev(:) !potential and kinetic energy on the current and previous time steps
    real*8, allocatable :: STEN(:,:,:)              !stress tensor (only xx, yy, zz components are set in the simulation
    character(len=:), allocatable :: OutDir         !directory for writing output snapshots
    integer, allocatable :: NNG(:), NNNG(:,:)       !number of neighbours, and neighbour list
    real*8, allocatable :: RhoDen(:)                !electron density, used for EAM
    
    real*8 :: CellSize(3), CellCenter(3)            !the size of the computational cell, and position of the cell center
    real*8 :: mass(KTMX)                            !mass of atoms
    real*8 :: time, T0, P0, dt                      !time, default temperature, default pressure, timeste
    integer :: Natoms, NRIGID, NSTEP                !number of atoms, number of ridged atoms, number of steps
    integer :: Boundary(3)                          !boundary condition (0 - free, 1 - periodic)
    integer :: KBOUND, dim                          !type of boundary (0 - free, 1 - rigid), number of dimensions
    integer :: NEWTAB, NEPRT, NWRITE, Ngather       !frequency of neighbour list renewal, printing information, writing snapshots, and gathering atoms
    integer :: NTYPE, NPOTS, KFLAG                  !number of atoms and atom pairs, simulation mode (1 - Quench, 2 - Vel, 3 - Heating, 5 - Const T)
    integer :: LFLAG, IPCON, KEYBS                  !pressure control flag, pressure control method, type of the interatomic potential
    
    logical :: FullList                             !full list flag
    real*8 :: Rskin                                 !skin layer thickness, used in neighbour search
    type(Table) :: table_U(KPMX)                    !tables for pair interatomic interaction U
    real*8 :: SW_sig(KPMX), SW_lam(KPMX), SW_gam(KPMX), SW_eps(KPMX) !parameters of SW potential
    type(Table) :: table_F(KTMX), table_fe(KTMX)    !F and f tables for EAM potential
    type(Table) :: table_UA(KPMX), table_C(KPMX), table_g(KPMX) !tables for Tersoff potential
    real*8 :: TF_N(KPMX), TF_Beta(KPMX)
    
contains
    subroutine allocate_system()                    !allocate memory for arrays used in the simulation according to Natoms
        if(allocated(XD)) deallocate(XD)
        if(allocated(FD)) deallocate(FD)
        if(allocated(Q1D)) deallocate(Q1D)
        if(allocated(Q2D)) deallocate(Q2D)
        if(allocated(Q3D)) deallocate(Q3D)
        if(allocated(Q4D)) deallocate(Q4D)
        if(allocated(Q5D)) deallocate(Q5D)
        if(allocated(KTYPE)) deallocate(KTYPE)
        if(allocated(KRIGID)) deallocate(KRIGID)
        if(allocated(KHIST)) deallocate(KHIST)
        if(allocated(Ep)) deallocate(Ep)
        if(allocated(Ek)) deallocate(Ek)
        if(allocated(Ek_prev)) deallocate(Ek_prev)
        if(allocated(STEN)) deallocate(STEN)
        if(allocated(NNG)) deallocate(NNG)
        if(allocated(NNNG)) deallocate(NNNG)
        if(allocated(RhoDen)) deallocate(RhoDen)

        allocate(XD(3,Natoms), FD(3,Natoms), Q1D(3,Natoms), Q2D(3,Natoms), Q3D(3,Natoms), Q4D(3,Natoms), Q5D(3,Natoms))
        allocate(KTYPE(Natoms), KRIGID(Natoms), KHIST(Natoms), Ep(Natoms), Ek(Natoms), Ek_prev(Natoms), STEN(3,3,Natoms))
        allocate(NNG(Natoms), NNNG(MAXNNB,Natoms), RhoDen(Natoms))

        Ek_prev(1:Natoms) = 0.0d0
        Q2D(1:3,1:Natoms) = 0.0d0
        Q3D(1:3,1:Natoms) = 0.0d0
        Q4D(1:3,1:Natoms) = 0.0d0
        Q5D(1:3,1:Natoms) = 0.0d0

        return
    endsubroutine
endmodule