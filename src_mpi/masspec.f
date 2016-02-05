C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C     Programer:
C     Travis Kemper
C     Department of Materials Science and Engineering
C     University of FLorida
C     traviskemper@ufl.edu
C
C     Version 1.0 01/12/11 T. W. Kemper
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Print out molecules according to mass
!
      SUBROUTINE masspec
!
      USE ANALYSIS
      USE POTS
      USE MDSTP
      USE MPIvars
      USE SPECIF
      USE BEAM 
!     
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
      INTEGER :: N,Io,Jo,Oo,MT,Mi,Mf,NMU,IM,NMB,NM,Ii,Ji,Oi
     &     ,ECNT,ELCNT(RTYPES),NGM,NSM,NDGM,NDSM,M
      REAL*8 :: MOLM,COM(3),MVEL(3),VAVE,VSQR,MOLVEL
      CHARACTER(30) :: PREL
      CHARACTER(4), DIMENSION(2) :: PREF
      LOGICAL :: DEPMOL
!
      IF ( mynode.EQ.0 )THEN
      MTCNTi = 0 
      MNCNTi(:,:) = 0
      GMVAVE(:) = 0.0d0
!
!     Loop over all molecules and output mass 
!     WRITE(*,*) MOLCNTi,' molecules where formed'
      DO 1001 NMU=1,MOLCNTi
        Ii=MPNTi(NMU)
        Ji=MPNTi(NMU+1)-1
        ECNT = 0 
!       Check if part of substrate 
        DEPMOL = .FALSE.
        DO Oi=Ii,Ji
          IM=MOLSTi(Oi)
          IF ( ASUBi(IM) ) GOTO 1001
          ECNT = ECNT + 1
          ELISTi(ECNT) = IM
          IF( IM.GT.SUBNA)  DEPMOL = .TRUE.
        ENDDO     
        CALL PRNEL(ECNT,PREL)
        CALL COMAS(ECNT,COM,MOLM,MVEL)
!       Make sure not beam  not sure if this works remove for now
!        IF( DEPMOL.AND.MVEL(DEPD).LT.MNDEPV ) GOTO 1001
        VSQR = 0.0d0
        DO M = 1,3
          VSQR = VSQR + MVEL(M)*MVEL(M)
        ENDDO
        MOLVEL = SQRT(VSQR)
!       Preform molecular mass analysis on Unbonded and gas phase molecules
        DO MT=1,MTCNTi  !Check molecular mass agianst recorded molecules
          Mi=INT(MOLM*1000)
          Mf=INT(MMASi(MT)*1000)
          IF ( Mi.EQ.Mf ) THEN
            IF (COM(DEPD).GT.SUBMN.AND.COM(DEPD).LT.SUBCT ) THEN
              MNCNTi(1,MT) = MNCNTi(1,MT) + 1
              IF(DEPMOL) MNCNTi(3,MT) = MNCNTi(3,MT) + 1
            ELSE
              MNCNTi(2,MT) = MNCNTi(2,MT) + 1
              IF(DEPMOL) MNCNTi(4,MT) = MNCNTi(4,MT) + 1
              GMVAVE(MT) = GMVAVE(MT) + MOLVEL
            ENDIF
            GOTO 1001
          ENDIF
        ENDDO
!       If not previously encountered molecule start a new list
        MTCNTi = MTCNTi  + 1
        IF ( MTCNTi.GE.MMX) THEN
          WRITE(*,*) 'More than ',MMX,' molc types found'
     &    ,'Increase MTPS'
          STOP
        ENDIF
        IF (COM(DEPD).GT.SUBMN.AND.COM(DEPD).LT.SUBCT ) THEN
          MNCNTi(1,MT) = MNCNTi(1,MT) + 1
          IF(DEPMOL) MNCNTi(3,MT) = MNCNTi(3,MT) + 1
        ELSE
          MNCNTi(2,MT) = MNCNTi(2,MT) + 1
          IF(DEPMOL) MNCNTi(4,MT) = MNCNTi(4,MT) + 1
          GMVAVE(MT) = GMVAVE(MT) + MOLVEL 
        ENDIF
        MMFRMi(MTCNTi) = PREL
        MMASi(MTCNTi) = MOLM
 1001 CONTINUE
!
      NMB = 0
      NMU = 0
      N = LSTEP / NSTR + 1
      DO MT=1,MTCNTi
        NSM =  MNCNTi(1,MT)      !Number of molc of this in sub
        NGM =  MNCNTi(2,MT)        !Number of molc of this in gas
        NDSM =  MNCNTi(3,MT)       !Number of molc of this in sub
        NDGM =  MNCNTi(4,MT)       !Number of molc of this in gas
        PREL = MMFRMi(MT)
        MOLM = MMASi(MT)
        VAVE = 0.0d0 
        IF (NGM.GT.0 ) VAVE = GMVAVE(MT)/FLOAT(NGM)
        WRITE(28,2701) N,TTIME,MT,PREL,MOLM,NGM,NSM,NDGM,NDSM,VAVE
      ENDDO
!
      ENDIF
      RETURN
 2701 FORMAT(I10,F14.2,I6,' ',A35,F10.6,' AMU',4I8,F12.6)
      END SUBROUTINE masspec
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Programer:
C     Travis Kemper
C     Department of Materials Science and Engineering
C     University of FLorida
C     traviskemper@ufl.edu
C
C     Version 1.0 02/24/11 T. W. Kemper
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Print out molecules which have formed on surface
!
      SUBROUTINE surfmspec
!
      USE ANALYSIS
      USE POTS
      USE MDSTP
      USE STRUCTR
      USE SPECIF
      USE MPIvars
!
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
      INTEGER ::  N,Io,Jo,Oo,MT,Mi,Mf,NMU,IM,NMB,NM,Ii,Ji,Oi
     & ,ECNT,ELCNT(RTYPES),PCNT,SCNT,I
      REAL*8 :: MOLM,COM(3),DELM,MVEL(3)
      CHARACTER(30) :: PREL
      CHARACTER(16), DIMENSION(:), ALLOCATABLE :: PRAN
      CHARACTER(4), DIMENSION(2) :: PREF
      CHARACTER(8) :: FLN
      LOGICAL :: PRIM
!
      IF ( mynode.EQ.0 )THEN
      SMTCNTi = 0 
!
!     Loop over all molecules and output mass 
  !    WRITE(*,*) MOLCNTi,' molecules where formed'
      DO 1001 NMU=1,SMCNTi
        Ii=SMPNTi(NMU)
        Ji=SMPNTi(NMU+1)-1
        ECNT = 0 
        IF( Ji-Ii .GT.STHS*10 ) THEN
          WRITE(*,*) 'Warning the # of molecules in suface molecule ',
     &           NMU,' is greater than STHS*10'
          STOP
        ENDIF
        DO Oi=Ii,Ji
          IM=SMOLSTi(Oi)
          ECNT = ECNT + 1
          ELISTi(ECNT) = IM
        ENDDO
        CALL PRNEL(ECNT,PREL)
        CALL COMAS(ECNT,COM,MOLM,MVEL)
!        WRITE(29,*) Ii,Ji,'Molecule ',NMU,PREL,MOLM
!       Preform molecular mass analysis on Unbonded and gas phase molecules
        DO MT=1,SMTCNTi  !Check molecular mass agianst recorded molecules
          Mi=INT(MOLM*1000)
          Mf=INT(SMMASi(MT)*1000)
          IF ( Mi.EQ.Mf ) THEN
            SMNCNTi(MT) = SMNCNTi(MT) + 1
            DO Oi=Ii,Ji
              IM=SMOLSTi(Oi)
              SMIMOLi(IM) = MT
            ENDDO     
            GOTO 1001
          ENDIF
        ENDDO
!       If not previously encountered molecule start a new list
        SMTCNTi = SMTCNTi  + 1
        IF ( SMTCNTi.GE.MMX) THEN
          WRITE(*,*) 'More than ',MMX,' molc types found'
     &    ,'Increase MTPS'
          STOP
        ENDIF
        SMNCNTi(SMTCNTi) = 1
        SMMFRMi(SMTCNTi) = PREL
        SMMASi(SMTCNTi) = MOLM
        DO Oi=Ii,Ji
          IM=SMOLSTi(Oi)
          SMIMOLi(IM) = SMTCNTi
        ENDDO 
 1001 CONTINUE
!
      NMB = 0
      NMU = 0
      N = LSTEP / NSTR + 1
      DO MT=1,SMTCNTi
   !     NMB =  MNCNT(1,MT)       !Number of molc of this type bonded 
        NM =  SMNCNTi(MT)       !Number of molc of this type unbonded 
   !     NM =  NMB + NMU + NMG(MT)
        PREL = SMMFRMi(MT)
        MOLM = SMMASi(MT)
        WRITE(29,2701) N,TTIME,MT,PREL,MOLM,NM
      ENDDO
!
      ENDIF
      RETURN
 2701 FORMAT(I10,F14.2,I6,' ',A35,F10.6,' AMU',4I6)
      END SUBROUTINE surfmspec

!
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Programer:
C     Travis Kemper
C     Department of Materials Science and Engineering
C     University of FLorida
C     traviskemper@ufl.edu
C
C     Version 1.0 01/12/11 T. W. Kemper 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Find crosslinks
!
      SUBROUTINE crslink
!
      USE MPIvars
      USE MDSTP
      USE ANALYSIS
      USE SPECIF
!
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
      INTEGER :: Ii,Ji,NMi,NMo,IM,Io,Jo,NAi,NAo,DA,CRSCNT,Oi
      LOGICAL :: FREF
 !
      IF ( mynode.EQ.0 )THEN
        CRSCNT = 0
        DO NMi=1,MOLCNTi
          Ii=MPNTi(NMi)
          IM=MOLSTi(Ii)
          IF ( ASUBi(IM) ) THEN
            Ji=MPNTi(NMi+1)-1
!           Find atom of chain that was also previously substrate
            FREF = .FALSE.
            DO Oi=Ii,Ji
                 IM=MOLSTi(Oi)
                 IF ( .NOT.ASUBo(IM) ) THEN
                   FREF = .TRUE.
                   EXIT
                 ENDIF
            ENDDO
            IF( FREF ) THEN
              NMo = IMOLo(IM)
              Io  = MPNTo(NMo)
              Jo = MPNTo(NMo+1)-1
              NAi = Ji - Ii + 1
              NAo = Jo - Io  + 1
              DA  = NAi - NAo
              IF( DA.GE.STHS) THEN
                CRSCNT = CRSCNT + 1
! begin debug
!              WRITE(*,*) NMi,NMo, NAi,NAo,DA
! end debug
              ENDIF 
            ENDIF 
          ENDIF 
        ENDDO
        WRITE(32,321) TTIME,LSTEP,CRSCNT
      ENDIF 
 321  FORMAT(F16.1,2I12)
      RETURN
      END SUBROUTINE crslink
