C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                              C
C  THESE PROGRAMS PERFORM MOLECULAR DYNAMICS SIMULATIONS WITH  C
C  BOND-ORDER POTENTIALS FOR HYDROCARBON AND OXYGEN,           C
C  AND LENNARD-JONES WITH INPUTTED PARAMETERS.                 C
C                                                              C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C   Units: mass = AMU's, length = Angstroms, time = fs         C
C          energy = eV's                                       C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                              C
C   Details about development and properties of potential      C
C   is described at :                                          C
C   1. CH                                                      C
C      D.W. Brenner*, O.A. Shenderova,                         C
C      J.A. Harrison, S.J. Stewart, B. Ni, S.B. Sinnott,       C
C      Journal of Physics: Condensed Matter 14,                C
C      783-802 (2002).                                         C
C   2. CHO                                                     C
C      Boris Ni, Ki-Ho Lee and Susan B. Sinnott                C
C      J. Phys.: Condens. Matter 16 (2004) 7261-7275           C
C                                                              C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C     Programer:                                               C
C      Travis Kemper                                           C
C      traviskemper@ufl.ed                                     C
C    Advisor:                                                  C
C      Dr. Susan Sinnott                                       C
C      ssinn@mse.ufl.edu                                       C
C      Department of Materials Science and Engineering         C
C      University of FLorida                                   C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C

      INCLUDE 'module_TM.f'
      INCLUDE 'module_parameters.f'
      INCLUDE 'module_specif.f'
      INCLUDE 'module_structure.f'
      INCLUDE 'module_pot.f'
      INCLUDE 'module_tables.f'
      INCLUDE 'module_contm.f'
      INCLUDE 'module_md.f'

      PROGRAM rebomd 
! 
      USE TM
      USE PARAMS
      USE SPECIF 
      USE STRUCTR
      USE POTS
      USE TABS
      USE CONTM
      USE MDSTP
!
      IMPLICIT none
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Initialize simulation 
      CALL TIMEI      !Caclculate initial cpu time 
      CALL WRTHEAD    !Write head of output file
      CALL STDFLT     !Set simulations variables to default values
c open input/output files 
      CALL OPENIN
c initialize 
      call setin 
C read input data 
      CALL READIN
      call read_coord
      CALL OPENOUT
C setup potential parameters  
      call setpp 
c setup predictor-corrector coefficients
      call setpc 
C setup Langevin parameters  
      call setgle 
C initialize random number generator 
      call setran 
! 
      call setmd 

! Initialize cfg 
      IF(PCFG) CALL SETCFG
!     Write simulation specifications
      CALL WRITEOPT
c********************
c begin calculation *
c********************
      CALL WRITE_COORD 
      IF(PXMOL) CALL XMOL 
      IF(PCFG)  CALL CFG

      IF ( MDRUN ) THEN 
        WRITE(6,*) 'Starting molecular dynamics simulation'
        CALL MODEL 
        WRITE(6,*) 'Molecular dynamics simulation complete'
      ENDIF
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     End Calculation  
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      CALL TIMEF    !Caclculate finial cpu time 
!      CALL close
C
      STOP
      END  PROGRAM rebomd
C
C add included subroutines
C
      include 'subroutines.inc' 
C Analysis subroutines

