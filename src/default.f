C Subroutines to set default parameters i
C   which can be overwriten in the input file 
!
      SUBROUTINE STDFLT
!
      USE SPECIF
      USE PARAMS
      USE POTS
      USE STRUCTR
      USE MDSTP
!
      IMPLICIT none
!
!     PI=ACOS(-1.D0) 
      ITITLE = 'Calculation unknown'
      CTITLE = 'System unknown'
      IPOT   = 1           !Potential type, 1=REBO
      KVC    = 1500        !Number of steps in entire simulation
      MAXKB  = 100         !# of steps between data output
      NSTR = 1000          !# of structure writes
      PSEED  = 0.2         !Random seed number
      RLL    = 0.5         !Neighbor list
      TEM    = 297.0       !Temperature 
      IRSVDT = .FALSE.     !Variable time step
      DISMAX = 0.12768885  !Max displacement 
      TTIME  = 0.0d0       !Simulation start time
      DELTA  = 0.2d0       !Time step
!
      STPT = .FALSE.       !Read in potential parameters from pots.in
      LTORS = .TRUE.      !Use torsion
      ILJ   = .TRUE.       !Included Lenard-Jones interactions
!
      ENMIN = .FALSE.      !Do energy minimization
      MDRUN = .TRUE.       !Do MD run
      COLID = .FALSE.      !Set predicitor corrector for molecular collisions/depositios
!
! Write initial structure
       WSTRI =.FALSE.
!
! Pibond
      SYTP = 1 !1=CH,2=CHO,3=CHS
!
! Dihedral Terms
      NDIHED = 2 !Include dihedral terms
C     NDIHED =10 !Don't include dihedral terms
!
!
! Thermostats
      KFLAG =  -1           !Thermostat  Langevin
C      KFLAG =   1          !Thermostat  Berendsen
C      KFLAG =   2          !Thermostat  Zero velocity
C      KFLAG =   3          !Thermostat  Evan-Hoover
 !
C Prodictor Corrector Coefficients
!
      F02 = 1.0d0/6.0d0
      F12 = 5.0d0/6.0d0
      F32 = 1.0d0/3.0d0
!
! Neighbor list
!      NLMAX = 500000   !Max # of atoms in neighbor list
!
      PRNPR = .FALSE.  !Print pair information between atom 1 and 2
      PRNF  = .FALSE.            !Print out forces on each atom 
!
      PR1 = 1
      PR2 = 2
!
! Integrator 
      INTG = 1  ! 1- 3rd order PC,2 - Velocity Verlet
!
! Print out frames in xmol format
      PXMOL = .FALSE.
!
! Print out frames if cfg format for atomeye
      PCFG = .FALSE.
!      
! Heat 
      HT = .FALSE. !Heat structure
      TMAX = 600.0d0
      DT = 50 !K/ps
! Thermostate options
      RTHRM = .FALSE.
      NTHRM = 1
!     Initialize heat bath
      BXTHRM = .FALSE.
      TBOX(1,:) = (/ 0.0d16,0.0d16,0.0d16 /)
      TBOX(2,:) = (/ 0.0d16,0.0d16,0.0d16 /)
      RBOX(1,:) = (/ 0.0d16,0.0d16,0.0d16 /)
      RBOX(2,:) = (/ 0.0d16,0.0d16,0.0d16 /)
!
! Moving atoms velocity
      tpvel(:) = (/0.0d0,0.0d0,0.d0/)
      VCONV = 1d-5 !m/s to A/fs
!
      END SUBROUTINE STDFLT

