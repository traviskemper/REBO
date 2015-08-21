!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Programer:
!     Travis Kemper
!     Department of Materials SCience and ENgineering
!     University of FLorida
!     traviskemper@ufl.edu
!
!     Version 1.0 4/19/11 T. W. Kemper
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Close files and deallocate arrays
!
      SUBROUTINE close 
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
      CLOSE(10)
      CLOSE(11)
      CLOSE(13)
      CLOSE(21)
      CLOSE(14)
      CLOSE(15)
      CLOSE(16)
      CLOSE(17)
      CLOSE(18)
      IF(PXMOL) CLOSE(1) 
!
      DEALLOCATE(AD
     & ,AXL,BD,BXL,CD,CXL
     & ,DD,DXL,ED,RB1,PID,RMAX,CHI
     & ,RLIST,XHC,DDTAB
     & ,tabfc,tabdfc,atable,datable,rtable,drtable
     & ,XM,XMM,XMMS,RSPL,RMAXLJ,RSLJ,RSPLS
     & ,C2,C3
     & ,vlook,dlook,iv,jv,RB2,SIG
     & ,EPS ,TAU,EPSS,SIGS,NOA,GL, EATOM
     & ,NABORS,LIST ,NABORSi,LISTi,RCLIST,IMOLo,IMOLi
     & ,SIMOLo,SIMOLi
     & ,IVCT2B,JVCT2B,LCHECK
     & ,COR,RCOR
     & ,WW,DWW,EXX1,DEXX1
     & ,KTYPE
     & ,R0,R1,R2,R3,R4,ITR
     & ,MLIST,NLIST,R0L,RNP
     & )
!
!
      RETURN
      END SUBROUTINE close
