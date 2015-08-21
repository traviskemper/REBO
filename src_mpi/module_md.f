C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C     Version 1.0 01/27/2010 T. W. Kemper                      C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C  Module of variablles for molecular dynamics steps 

      MODULE MDSTP 
 
      INTEGER :: LSTEP,KEND,KCHK
     &,KLIEND ,NNBI     
     &,KB       !Some count for load??
      REAL*8 :: TOTE,S3,T1,T2,TTOTE
     &,pvdw,RSS(3),RRS(3),dellj   !LJ 
      REAL*8, DIMENSION(:), ALLOCATABLE :: EATOM
      INTEGER , DIMENSION(:), ALLOCATABLE :: IVCT2B,JVCT2B  !Neighbor list
     & ,NABORSi(:),LISTi(:),RCLIST(:),IMOLi(:),IMOLo(:)
     & ,SIMOLi(:)
     &,IV,JV          !LJ

      END MODULE MDSTP



