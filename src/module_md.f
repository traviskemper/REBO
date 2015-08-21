!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
!     Version 1.0 01/27/2010 T. W. Kemper                      C
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
!  Module of variablles for molecular dynamics steps 
!
      MODULE MDSTP 
!
      INTEGER :: LSTEP,KEND,KCHK
     &,KLIEND ,NNBi             
      REAL*8 :: TOTE,S3,T1,T2,TEMPK,TTIME
     &,pvdw,RSS(3),RRS(3),dellj   !LJ 
      REAL*8, DIMENSION(:), ALLOCATABLE :: EATOM
     &,RCOR  !Neighbor list
     &,WW,DWW,EXX1,DEXX1
     &,IV,JV          !LJ
      INTEGER , DIMENSION(:), ALLOCATABLE :: NABORS,LIST,LCHECK
     & ,IVCT2B,JVCT2B,NABORSi(:),LISTi(:),RCLIST(:),IMOLo(:),IMOLi(:)
     & ,SIMOLi(:),SIMOLo(:)
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: RNP
     &,COR    !Neighbor list
     &,XHC
      INTEGER :: LCHK
!
      END MODULE MDSTP



