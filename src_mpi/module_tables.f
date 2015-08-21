C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C     Version 1.0 01/25/2010 T. W. Kemper                      C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C Look up tables module

      MODULE TABS

     
      INTEGER :: NTAB
      PARAMETER( NTAB =10000 )
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: 
     &  tabfc,tabdfc,atable,datable,rtable,drtable
     &,vlook,dlook
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: DDTAB
     &,XMM,XMMS,XM
     &,RMAXLJ,RSLJ,RSPL,RSPLS
     &,c2,c3

      END MODULE TABS
