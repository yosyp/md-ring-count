program MSE627_CG
!   This program builds one of the following crystalline structures:
!   (1) simple cubic
!   (2) body centered cubic
!   (3) body centered cubic with 2 types of atoms (CsCl)
!   (4) face centered cubic
!   (5) face centered cubic with 2 types of atoms (Ni3Al,Cu3Au)
!   (6) diamond strcture (Si)
!   (7) diamond strcture with 2 types of atoms (SiC, GaAs)
!   All the structures are described as simple cubic Bravais lattices
!   with a basis consisting of up to 8 atoms.
!   See Ashcroft & Mermin for reference.
!   A shift to have origin of the computational cell at (0,0,0).
    
    implicit none
!   Defining input parameters (you can change this part to generate a particular system)
    character(*), parameter :: outfile = "ArNew.data"   !Name of the output file
    integer, parameter :: nn(3) = (/5, 5, 5/)           !Number of unit cells in x,y,z
    real*8, parameter :: alat = 5.405d0                 !Lattice constant in Ang
    integer, parameter :: ltype = 4                     !Lattice type
    
    integer, parameter :: LPMX = 5000       !Maximum number of particles
    integer, parameter :: dim = 3
    integer, parameter :: lat_max = 8       !The maximum number of atoms in the basis
    character*5, parameter :: Crtype(7) = (/' SC  ',' BCC ','CsCl ',' FCC ','FCC31',' DL  ','GaAs '/)
    
    real*8 :: pos(dim,LPMX), vel(dim,LPMX)  !Coordinates, velocities
    integer :: KTYPE(LPMX), KHIST(LPMX)     !Type & History of each atom
    integer :: Natoms                       !Number of particles
    real*8 :: time                          !Starting time of the simulation
    real*8 :: CellSize(dim), CellCenter(dim)!Size of the computational cell & Position of the center of C.C.
    integer :: i, j, k, ic, atom_id, file_id
    
    integer :: ncell                        !Number of particles in the unit cell
    real*8 :: basis(dim, lat_max)           !Positions of atoms in unit cell
    integer:: ntype(lat_max)                !Number of types of atoms in the basis
    real*8 :: shift(dim)                    !Shift of the crystallite

    
!   Setting parameters by type
    select case(ltype)
        case (1)                            !Simple Cubic
            ncell = 1
            ntype(1) = 1
            basis(1:3,1) = 0.0d0
            shift(1:3) = 0.5d0
        case (2)                            !Body Centered Cubic
            ncell = 2
            ntype(1:2) = 1
            basis(1:3,1) = 0.0d0
            basis(1:3,2) = 0.5d0
            shift(1:3) = 0.25d0
        case (3)                            !CsCl - 2 atom bcc
            ncell = 2
            ntype(1:2) = (/1, 2/)
            basis(1:3,1) = 0.0d0
            basis(1:3,2) = 0.5d0
            shift(1:3) = 0.25d0
        case (4)                            !Face Centered Cubic
            ncell = 4
            ntype(1:4) = 1
            basis(1:3,1) = (/0.0d0, 0.0d0, 0.0d0/)
            basis(1:3,2) = (/0.0d0, 0.5d0, 0.5d0/)
            basis(1:3,3) = (/0.5d0, 0.0d0, 0.5d0/)
            basis(1:3,4) = (/0.5d0, 0.5d0, 0.0d0/)
            shift(1:3) = 0.25d0
        case (5)                            !Ni3Al-type FCC
            ncell = 4
            ntype(1:4) = (/2, 1, 1, 1/)
            basis(1:3,1) = (/0.0d0, 0.0d0, 0.0d0/)
            basis(1:3,2) = (/0.0d0, 0.5d0, 0.5d0/)
            basis(1:3,3) = (/0.5d0, 0.0d0, 0.5d0/)
            basis(1:3,4) = (/0.5d0, 0.5d0, 0.0d0/)
            shift(1:3) = 0.25d0
        case (6)                            !Diamond lattice
            ncell = 8
            ntype(1:8) = 1
            basis(1:3,1) = (/0.0d0, 0.0d0, 0.0d0/)
            basis(1:3,2) = (/0.0d0, 0.5d0, 0.5d0/)
            basis(1:3,3) = (/0.5d0, 0.0d0, 0.5d0/)
            basis(1:3,4) = (/0.5d0, 0.5d0, 0.0d0/)
            basis(1:3,5:8) = basis(1:3,1:4) + 0.25d0
            shift(1:3) = 0.125d0
        case (7)
            ncell = 8
            ntype(1:4) = 1
            ntype(5:8) = 2
            basis(1:3,1) = (/0.0d0, 0.0d0, 0.0d0/)
            basis(1:3,2) = (/0.0d0, 0.5d0, 0.5d0/)
            basis(1:3,3) = (/0.5d0, 0.0d0, 0.5d0/)
            basis(1:3,4) = (/0.5d0, 0.5d0, 0.0d0/)
            basis(1:3,5:8) = basis(1:3,1:4) + 0.25d0
            shift(1:3) = 0.125d0
        case default
            write(*,"('ltype = ', I3, ' is not implemented')") ltype
    endselect
    
    CellSize(1:dim) = nn(1:3)*alat          !Define the sizes of the computational cell
    CellCenter(1:dim) = 0.5d0*nn(1:3)*alat
    Natoms = product(nn(1:3))*ncell
    
    write(*,"('Building a crystal with nx = ',i3', ny = ',i3', nz = ',i3)") nn(1:dim)
    write(*,"('lattice constant = ',f8.5,' A, lattice type: ',A)") alat, Crtype(ltype)
    write(*,"('The system sizes are ',2(f8.5,' A, '),f8.5,' A')") CellSize(1:dim)
    write(*,"('The number of atoms is ', i5)") Natoms
    
    if(Natoms > LPMX) then
        write(*,"('Number of atoms ',i5,' is larger than the limit ',i5)") Natoms, LPMX
        stop
    endif

!   Generating the lattice
    atom_id = 1
    do i = 0, nn(1) - 1
        do j = 0, nn(2) - 1
            do k = 0, nn(3) - 1
                do ic = 1, ncell
                    ktype(atom_id) = ntype(ic)
                    pos(1:dim,atom_id) = basis(1:dim,ic) + shift(1:dim) + (/i,j,k/)
                    atom_id = atom_id + 1
                enddo
            enddo
        enddo
    enddo
    pos(1:dim,1:Natoms) = pos(1:dim,1:Natoms)*alat  !Scale by lattice constant after generating lattice
    KHIST(1:Natoms) = 0.0d0
    vel(1:dim,1:Natoms) = 0.0d0
    time = 0.0d0
    
!   Defining the rigid layers and "isotopes" for diffusion simulation (Homework #5)
!    do i = 1, Natoms
!        if(pos(3,i) > CellCenter(3)) KHIST(i) = 1
!        if(pos(3,i) < (alat + 0.01d0)) KHIST(i) = 3
!        if(pos(3,i) > (CellSize(3) - alat - 0.01d0)) KHIST(i) = 3
!    enddo

!   Write output coordinate file
    open (newunit = file_id, file = outfile)
    write(file_id,*) Natoms, time
    write(file_id,*) CellSize(1:3), CellCenter(1:3)
    write(file_id,*) (KTYPE(i),i = 1, Natoms)
    write(file_id,*) (KHIST(i), pos(1:3,i), i = 1, Natoms)
    write(file_id,*) (vel(1:3,i), i = 1, Natoms)
    close(file_id)
 
!   Write a file for system visualization
!    open (newunit = file_id, file = "snapshot.d")
!    do i = 1, Natoms
!       write(file_id,*) pos(1:3,i), KTYPE(i), KHIST(i)
!    enddo
!    close(file_id)
    
endprogram