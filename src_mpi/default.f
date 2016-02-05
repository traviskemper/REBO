C Subroutines to set default parameters i
C   which can be overwriten in the input file 

      SUBROUTINE STDFLT
! 
      USE MPIvars
      USE STRUCTR
      USE SPECIF
      USE PARAMS
      USE ANALYSIS
      USE BEAM
      USE POTS
!
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
      IF ( mynode.EQ.0 ) THEN
C        PI=ACOS(-1.D0) 
!         Potential
          IPOT = 1
          ILJ  = .TRUE.
!         MD run
          CTITLE = 'System unknown'
          KVC = 1500           !Number of steps in entire simulation
          MAXKB = 100          !# of steps between data output
          NSTR = 1000          !# of xmol.d writes
          PSEED = 0.109824     !Random seed number
          RLL = 0.5            !Neighbor list
          TEM = 297.0          !Temperature 
          IRSVDT = .FALSE.     !Variable time step
          DISMAX = 0.12768885d0  !Max displacement 
          TTIME  = 0.0d0       !Simulation start time
          DELTA  = 0.2d0       !Time step
!
          STPT = .FALSE.       !Read in potential parameters from pots.in
          LTORS = .FALSE.      !Use torsion
!           LLJ   = .FALSE.      !Included Lenard-Jones interactions
          ENMIN = .FALSE.      !Do energy minimization
          MDRUN = .TRUE.       !Do MD run
          COLID = .FALSE.      !Set predicitor corrector for molecular collisions/depositios
!
C Pibond
          PRNLD = .FALSE.
!
C Dihedral Terms
C      NDIHED = 2 !Include dihedral terms
C     NDIHED =10 !Don't include dihedral terms

C Thermostats
          KFLAG =  -1          !Thermostat  Langevin
C      KFLAG =   1          !Thermostat  Berendsen
C      KFLAG =   2          !Thermostat  Zero velocity
C      KFLAG =   3          !Thermostat  Evan-Hoover
!
C Prodictor Corrector Coefficients
           F02 = 1.0d0/6.0d0
           F12 = 5.0d0/6.0d0
           F32 = 1.0d0/3.0d0
!
C Neighbor list
C      NLMAX = 500000   !Max # of atoms in neighbor list
          PRNPR = .FALSE.  !Print pair information between atom 1 and 2
          PRNF  = .FALSE.            !Print out forces on each atom !
          PR1 = 1
          PR2 = 2
! Moving atoms velocity
          tpvel(:) = (/0.0d0,0.0d0,0.d0/)
          VCONV = 1d-5 !m/s to A/fs
          MVTIP = .FALSE.
          ISTEPS = 1
          disp(:) =(/0.0d0,0.d0,0.d0/)
          PRNLD = .FALSE.
!
C Integrator 
          INTG = 1  ! 1- 3rd order PC,2 - Velocity Verlet
                
C Heat       
          HT = .FALSE. !Heat structure
          TMAX = 600.0d0
          DT = 50 !K/ps
! Beam 
      MKBM = .FALSE.
      CENTR = .FALSE.
      SHFT = .FALSE.
!     Set buff area between heat bath and deposited atoms
      THBUF = 5.0  !Angstroms
!     Set dimension of deposition
      DEPD=3  !Z
      DBUF = 7.0  !Hight of beam off substrate 
      MNDEPV = -0.001d0 ! minimum velocity to be considered a beam molecule

!
! Analysis
      DEPAN = .FALSE. !Moleculare analysis 
      MOLAN = .FALSE. !Moleculare analysis 
      REFSTR = .FALSE. !Use referance strucuture as point of comparison
      RECAN = .FALSE. !Moleculare analysis 
!     Set buffer region
      HBUF = 10.  !Max hight to be counted as substrate 
! Clear gas molecules
      CLRGAS = .FALSE.
! Atomic depth profile
      DEPPROF = .FALSE.
      DPDIV  = 5.0d0
      
!
      ENDIF
!
      END SUBROUTINE STDFLT

      
