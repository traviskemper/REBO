C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C     Version 1.0 01/22/2010 T. W. Kemper                      C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C  Module of potential parameters

      MODULE POTS

C Atom type numbers
      INTEGER :: icarb,ihyd,ioxygen,isulfur,iflor,isi,iger
     &,RTYPES,iarg,NTYPES,ikryp,ineon
      PARAMETER (  RTYPES = 5,NTYPES=6 )
      PARAMETER( icarb =1, ihyd = 2
     &,isulfur = 3
     &,ioxygen = 4
     &,iflor   = 5
     &,iarg    = 6
     &,ineon   = 7
     &,ikryp   = 8
     &,isi     = 11
     &,iger    = 12
     & )
      REAL*8 :: XMASS(100)
      CHARACTER(2), DIMENSION(100) :: ATYPE
! 
C Data statements from REBO CH
 
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: AD,AXL,BD,BXL
     & ,CD,CXL,DD,DXL,ED,RB1,RB2,PID,RMAX,CHI,RLIST,SIG,EPS
      REAL*8 :: 
     & XN1(2)
     &,SPGC(6,5),SPGH(6,3),XQ,ATT,XQM,PQ,SPGS(RTYPES,7)
     &,XTN2(4),XTN1(4),ADB(4),CDB(4),CDB2(4),DDB(4)
     &,DDB2(4),HDB(4)
     &,XXDB,XDB(RTYPES,RTYPES,RTYPES),REG(RTYPES,RTYPES,RTYPES)
     &,RHH,RCH,ROH
     &,RLL
      INTEGER :: IGC(25),IGH(25),KT(100),KT2(100)

C Spline 
      REAL*8 :: CLMN(3,10,10,10,64),CLM(2,10,10,16)
     &,TLMN(10,10,10,64),IN2(16,2),IN3(64,3),PIDT
     &,XH(2,10,10),XH1(2,10,10),XH2(2,10,10)
     &,XHO(RTYPES,RTYPES,10,10),XHO1(RTYPES,RTYPES,10,10)
     & ,XHO2(RTYPES,RTYPES,10,10),XHO12(RTYPES,RTYPES,10,10)
     &,CLMOX(RTYPES,RTYPES,10,10,16)
! Fluorine splines
     &,CLMF(RTYPES,10,10,10,64)
     &,XHF(RTYPES,10,10,10),XHF1(RTYPES,10,10,10)
     &,XHF2(RTYPES,10,10,10),XHF3(RTYPES,10,10,10)
!
! Sulfur Sline
     &,XHS(RTYPES,RTYPES,0:9,0:9,0:9),XHS1(RTYPES,RTYPES,0:9,0:9,0:9)
     &,XHS2(RTYPES,RTYPES,0:9,0:9,0:9),XHS3(RTYPES,RTYPES,0:9,0:9,0:9)
     &,XHS12(RTYPES,RTYPES,0:9,0:9,0:9),XHS23(RTYPES,RTYPES,0:9,0:9,0:9)
     &,XHS13(RTYPES,RTYPES,0:9,0:9,0:9)
     &,XHS123(RTYPES,RTYPES,0:9,0:9,0:9)
     &,CLMSX(64,0:9,0:9,0:9,RTYPES,RTYPES)


C Oxygen angular terms
      REAL*8 :: a_o0,a_o1,a_o2

C Sulfur angular terms
      REAL*8 :: a_s0,a_s1,a_s2

C Potential type
      INTEGER :: SYTP

      END MODULE POTS



