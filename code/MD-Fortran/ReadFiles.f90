!   THIS SUBROUTINE READS MOST OF THE INPUT FILES
!   MSE 6270, Leonid Zhigilei

subroutine ReadFiles()
    use GlobalVars
    implicit none
    integer :: i
    
!   Read input data
    rewind 14
    read(14,*) NSTEP            !NSTEP (Number of steps)
    read(14,*) NEWTAB           !NEWTAB (Step of neighbour list renewal)
    read(14,*) NEPRT            !NEPRT (Step of printing output information)
    read(14,*) NWRITE           !NWRITE (Step of writing output information)
    read(14,*) Ngather          !Ngather (Step of gathering molecules to the comp. cell)
    read(14,*) KFLAG            !KFLAG (1-Quench, 2-Vel, 3-Heating, 5-Const T)
    read(14,*) LFLAG            !LFLAG (1-Berendsen Const Pressure)
    read(14,*) IPCON            !IPCON (P control: 1-3D P, 2-X,Y+Z, 3-X,Y,Z independent)
    read(14,*) KEYBS            !KEYBS (0-pair potential,1-Stillinger Weber)
    read(14,*) Boundary(1)      !LIDX (1-X periodic,0-free boundary conditions])
    read(14,*) Boundary(2)      !LIDY (1-Y periodic,0-free boundary conditions])
    read(14,*) Boundary(3)      !LIDZ (1-Z periodic,0-free boundary conditions])
    read(14,*) KBOUND           !KBOUND (Type of boundary 0-free,1-rigid,2-...)
    read(14,*) dim              !dim (3-3D simulation, 2-2D simulation)
    read(14,*) T0               !T0 (T distributed by VEL or controlled by Bere-T)
    read(14,*) P0               !P0 (P in GPa controlled by Bere-P)
    read(14,*) dt               !dt (Timestep of integration in pr.unit [psec])
    read(14,*) RSkin            !RSkin (Skin depth for neighbour list [A])

!   Read input coordinates/velocities file
    rewind 15
    read(15,*) Natoms, time
    call allocate_system()
    read(15,*) CellSize(1:3), CellCenter(1:3)
    read(15,*) (KTYPE(i), i = 1, Natoms)
    read(15,*) (KHIST(i), XD(1:3,i), i = 1, Natoms)
    read(15,*) (Q1D(1:3,i), i = 1, Natoms)
    Q1D(1:3,1:Natoms) = Q1D(1:3,1:Natoms)*dt !program uses v*dt for "velocity"
    
    return
endsubroutine










