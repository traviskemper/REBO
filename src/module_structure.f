C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C     Version 1.0 01/21/2010 T. W. Kemper                      C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C  Module of structure information 

      MODULE STRUCTR
 
      INTEGER :: NP,NMA,NTA,SSS
      REAL*8 :: DCONV
      PARAMETER ( DCONV = 1.6605386d0 )  ! g A^3/AMU cm^3 
      REAL*8, DIMENSION(3) :: CUBE,CUBE2,CUBE_r
      REAL*8, DIMENSION(:,:), ALLOCATABLE ::  
     &   R0,R1,R2,R3,R4,R0L
C      REAL*8, DIMENSION(:), ALLOCATABLE ::  ITR
      CHARACTER(70) :: CTITLE
      INTEGER, DIMENSION(:), ALLOCATABLE :: NOA,MLIST,NLIST,ITR
     &,KTYPE
      REAL*8 :: VOL,RRL,RMN(2,3),VOLMXN,MXNDENS,CDENS,MTOT
! Cfg data
      REAL*8 :: cbox(3)

      END MODULE STRUCTR 



