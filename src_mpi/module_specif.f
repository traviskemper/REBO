C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C     Version 1.0 01/21/2010 T. W. Kemper                      C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C  Module of simulations specifics  

      MODULE SPECIF
 
      INTEGER :: KVC,MAXKB,NDIHED,IPOT,NSTR
      REAL*8 :: PSEED,TEM,DISMAX,TTIME,DELTA,DELTSQ,TTCONV
     &,F02,F12,F32,F42,F52
C      LOGICAL :: IRSVDT,LTORS,LLJ,ENMIN,MDRUN,STPT,COLID
      CHARACTER(70) :: ITITLE
      LOGICAL :: LTORS,IRSVDT,ILJ,STPT,LJRUN,ENMIN,MDRUN,COLID
C Thermostat related variables 
      INTEGER :: KFLAG,NLR,NTHRM
      REAL*8  :: WD,TR,BET,GSIG,PI2,DNLA
     & ,TBOX(2,3),RBOX(2,3),RIG(2,3),THM(2,3)
      REAL*8, DIMENSION(:), ALLOCATABLE :: GL
      LOGICAL :: RTHRM,BXTHRM,CUSTHRM
!
C Random number genorator constants
      REAL*8 :: QBASE,QA1,QA2,QB1,QB2,RAN

C MD variables
      REAL*8 :: ENPR,XKEA,TIME
      INTEGER :: LCHK

C Pressing stuff 
      INTEGER :: ISTEPS
      REAL*8 ::  disp(3),tpvel(3),VCONV
      LOGICAL :: PRNLD,MVTIP

C Printing pair interactions
      LOGICAL :: PRNPR,PRNF
      INTEGER :: PR1,PR2
      REAL :: PREN,PRDIST

C Integrator
      INTEGER :: INTG

C Heat
      LOGICAL :: HT
      REAL*8 :: TMAX,DELTMP,DT
! Data outputs
      LOGICAL :: PXMOL,PCFG,WSTRI
! Add force 
      INTEGER :: FSTEP,SDIR,TATOMS,BATOMS
      REAL*8 :: ADFT,ADFB,PRADF,SAREA,PRBUF
      LOGICAL :: ADDFC,PTHRM
!     
      END MODULE SPECIF



