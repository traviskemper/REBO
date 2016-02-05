C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C     Version 1.0 06/22/2010 T. W. Kemper                      C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C  Module of analysis data
C
      MODULE ANALYSIS
C 
      USE PARAMS
C
      INTEGER :: RMX,STHS,MMX,NMX,ELMX,BTYPES
      PARAMETER ( MMX  = 10000 ) !Maximum number of molecule types 
      PARAMETER ( STHS = 100 ) !Number of atoms to be considered chain
      PARAMETER ( NMX  = 20 )  !Maximum number of neighbors 
      PARAMETER ( BTYPES = 10 ) 
      LOGICAL :: DEPAN,MOLAN,REFSTR,RECAN,RECSUB,SMOLAN
! Molecular analysis
      INTEGER :: MOLCNTi,MOLCNTo,SMCNTi
      INTEGER , ALLOCATABLE  :: MOLSTi(:),MPNTi(:),MOLSTo(:),MPNTo(:)
     & ,MIMOLi(:),SMPNTi(:),SMOLSTi(:),SMIMOLi(:)
      LOGICAL , DIMENSION(:), ALLOCATABLE :: ASUBi,ASUBo
! Mass spec
      INTEGER :: MNCNTi(4,MMX),MTCNTi,SMTCNTi,SMNCNTi(MMX)
      CHARACTER(30), DIMENSION(MMX):: MMFRMi,SMMFRMi
      REAL*8 :: MMASi(MMX),MMASo(MMX),SMMASi(MMX),GMVAVE(MMX)
! Reaction analysis
      INTEGER :: REC,NNBo
      INTEGER , ALLOCATABLE  ::NABORSo(:),LISTo(:),RECTi(:)
     & ,ELISTi(:),APLIST(:),ELISTo(:),ARLIST(:)
!
      INTEGER,ALLOCATABLE :: BID(:,:),BRB(:),FRB(:),EXB(:)
     &       ,BRBGS(:),FRBGS(:),EXBGS(:)
     &       ,BRBDP(:,:),FRBDP(:,:),EXBDP(:,:)
!     Substrate information
      CHARACTER(70) :: DTITLE 
      INTEGER :: SUBNA,NPS
      REAL*8 :: DBUF,HBUF,SUBMN,SUBMX,SUBTPEN,SUBTKEN,SUBTEMP,SUBCT
!
! Deposited analysis
      LOGICAL :: DEPPROF
      INTEGER :: VDEP,SDIV
      REAL*8 :: SMNo(2,3),SDEP,DPDIV
      INTEGER, ALLOCATABLE ::  GSA(:),SBA(:),DSBA(:,:) 
      REAL*8,ALLOCATABLE :: MXDEP(:)
! Clear gas
      REAL*8 :: GASD
      LOGICAL :: CLRGAS
! Crosslink analysis
      LOGICAL :: CRSAN
! Composition
      INTEGER, ALLOCATABLE :: ELCOMP(:,:)
      REAL*8, ALLOCATABLE :: GASKE(:)
!
      END MODULE  ANALYSIS
