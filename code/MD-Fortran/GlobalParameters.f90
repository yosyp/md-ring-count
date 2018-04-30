!A module with global parameters
module GlobalParameters
    implicit none
    
    integer, parameter :: KTMX = 3                  !Maximum number of particle types
    integer, parameter :: KPMX = KTMX*(KTMX+1)/2    !Maximum number of potential types
    integer, parameter :: MAXNNB = 250              !Max. number of neighbour pairs

    !---- Fundamental constants ----------------------------------------------
    real*8, parameter :: PI = 3.1415926535897932D+00 
    real*8, parameter :: EVTOJOU = 1.60219D-19      !J/eV
    real*8, parameter :: AMUTOKG = 1.6605402D-27    !kg/amu
    real*8, parameter :: BK = 8.617385D-05          !Boltzman constant, eV/K
    real*8, parameter :: XJOUTOEV = 1.0d0/EVTOJOU   !eV/J
    !-------------------------------------------------------------------------

    !---- Transfer to program units ------------------------------------------
    real*8, parameter :: ENUNIT = AMUTOKG*1.0d4*XJOUTOEV  !eV/pr.u.
    !-------------------------------------------------------------------------
    
!           KTYPE(J)    
!          | 1  2  3
!      K  ----------
!      T  1| 1  3  6
!      Y   |           IJINDEX(KTYPE(I),KTYPE(J)) = 1:Npots
!      P  2| 3  2  5
!      E   |
!     (I) 3| 6  5  4
    integer, parameter :: IJINDEX(KTMX,KTMX) = &    !indexes of the potential table
        reshape((/ 1, 3, 6, 3, 2, 5, 6, 5, 4 /), (/KTMX, KTMX/))
endmodule
