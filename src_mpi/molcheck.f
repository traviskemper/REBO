C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Programer:
C     Travis Kemper
C     Department of Materials SCience and ENgineering
C     University of FLorida
C     traviskemper@ufl.edu
C
C     Version 1.0 9/25/09 T. W. Kemper
C     Version 5.0 10/05/09 T. W. Kemper
C        - working independant version 
C     Version 6.2 10/06/09 T. W. Kemper
C        - copy of 5.0
C        - cut neighbor list genorator out to sepreate subroutine
C     Version 8 12/16/09 T. W. Kemper
C        - copy of 6.2
C        - fix cleanup (atoms being deleted)
C     Version 9 12/29/09 T. W. Kemper
C        - do not count H-H bonding 
C     Version 10 12/29/09 T. W. Kemper
C        - try new looping method 
C     Version 11 01/12/11 T. W. Kemper
C        - use rebo naborlist
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  Check if substrate atoms have formed any molecules
!
      SUBROUTINE molcheck
!
      USE ANALYSIS
      USE STRUCTR
      USE POTS
      USE MDSTP
      USE MPIvars
!
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
      INTEGER :: I,J,MI,MF,MJ,MCK,CK,A,B,C,N,II,NN
     & ,JBEGIN,JEND,CREF ,CNT,RPCNT
     & ,K,L,KI,KJ
      REAL*8 :: RIJSQ,DR,RCUTSQ
     & ,RI(3),RLIS,RSQ,RR(3)
      LOGICAL :: RPTCHK,VERBOSE,TERM,MLERR
      LOGICAL, ALLOCATABLE :: MOL(:)

      IF ( mynode.EQ.0 )THEN
!          Create neighbor list based on R0i
           K = 0
           DO 302 I=1,NP
                NABORSi(I)=K+1
                DO 299 L=1,3
                     RI(L)=R0i(L,I)
299             CONTINUE
                KI=KTYPEi(I)
                DO 301 J=1,NP
                     IF(I.EQ.J) GO TO 301
                     KJ=KTYPEi(J)
                     RLIS=RMIN(KI,KJ)
                     RSQ=0.0D0
                     DO 298 L=1,3
                          RR(L)=RI(L)-R0i(L,J)
                          RR(L)=RR(L) -
     &                          CUBE(L)*ANINT(RR(L)/CUBE(L))
                          RSQ=RSQ+RR(L)*RR(L)
                          IF(RSQ.GT.RLIS) GO TO 301
298                  CONTINUE
C
405                  CONTINUE
                     K=K+1
                     LISTi(K)=J
C
301             CONTINUE
302        CONTINUE
C
           NABORSi(NP+1)=K+1
!
      ALLOCATE( MOL(NP) ) 
!     Save previous list as o
      MOLCNTo = MOLCNTi
      MPNTo(:)   = MPNTi(:)
      MOLSTo(:)  = MOLSTi(:)
      ASUBo(:)   = ASUBi(:)
      IMOLo(:)   = IMOLi(:)
!     Zero out lists 
      MPNTi(:) = 0
      MOLSTi(:) = 0 
      IMOLi(:) = 0
C Set verbose output for debuging
      VERBOSE=.FALSE.
      MOL(:) = .FALSE.
      ASUBi(:) = .FALSE.
C
C Loop over all atoms and determine molecules
      MOLCNTi = 0
      CNT = 0
      DO I=1,NP
       RECTi(I) = 0 !initialize for react.f 
       IF ( .NOT.MOL(I).AND.KTYPEi(I).NE.ihyd ) THEN
        MOLCNTi=MOLCNTi+1
        CNT=CNT+1
        IF (CNT.GT.NP) THEN
         WRITE(13,*) 'In main molecule loop'
         WRITE(13,*) '# of atoms in molecules > max# of atoms',NP
         WRITE(13,*) 'When adding atom',I,CNT,R0i(:,I)
         WRITE(13,*)'Change dimension of MOLSTi in molcheck'
         WRITE(13,*)'but something else is probably wrong'
         OPEN(UNIT=87,file='mollist.err',status='unknown')
         DO N=1,MOLCNTi
             MI = MPNTi(N)
             MF = MPNTi(N+1)-1
             DO MJ=MI,MF
                II = MOLSTi(MJ)
                WRITE(87,*) N,II,MOL(II)
             ENDDO
         ENDDO 
         CLOSE(87)
         STOP
       ENDIF
       MOLSTi(CNT) = I
       MPNTi(MOLCNTi) = CNT
       MOL(I) =.TRUE.
       IMOLi(I) = MOLCNTi
C
C while termination is not present add non hydrogen neighbors to 
C  molecular/chain grouping 
       TERM=.TRUE.
       DO WHILE (TERM)
        MI = MPNTi(MOLCNTi)
        MF = CNT
        CREF = CNT
        DO MJ = MI,MF
          A = MOLSTi(MJ)
          JBEGIN=NABORSi(A)
          JEND=NABORSi(A+1)-1
          DO 1001 N = JBEGIN,JEND
             B = LISTi(N)
             IF ( MOL(B) ) GOTO 1001
             DO MCK = MI,CNT
              CK = MOLSTi(MCK)
              IF ( B.eq.CK.OR.KTYPEi(B).EQ.ihyd) GOTO 1001 
             ENDDO
             CNT = CNT + 1
             MOLSTi(CNT) = B
             MOL(B) =.TRUE.
             IMOLi(B) = MOLCNTi
 1001     CONTINUE
        ENDDO
        IF( CNT.EQ.CREF) TERM=.FALSE.
       ENDDO
C
C Add hydrogens to list
       MI = MPNTi(MOLCNTi)
       MF = CNT
       DO  MJ=MI,MF
         A = MOLSTi(MJ)
         JBEGIN=NABORSi(A)
         JEND=NABORSi(A+1)-1
         DO 2001 N = JBEGIN,JEND
            B=LISTi(N)
            IF( .NOT.MOL(B).AND. KTYPEi(B).EQ.ihyd ) THEN
              CNT = CNT + 1
              MOLSTi(CNT) = B
              MOL(B) =.TRUE.     
              IMOLi(B) = MOLCNTi
            ENDIF 
 2001    CONTINUE
        ENDDO
       ENDIF
      ENDDO
c$$$! begin debug
c$$$           OPEN(UNIT=88,file='mollist.dat',status='unknown')
c$$$           DO NN=1,MOLCNTi
c$$$             MI = MPNTi(NN)
c$$$             MF = MPNTi(NN+1)-1
c$$$             DO MJ=MI,MF
c$$$                II = MOLSTi(MJ)
c$$$                WRITE(88,*) NN,II,MJ,MOL(II)
c$$$             ENDDO
c$$$           ENDDO 
c$$$           OPEN(UNIT=89,file='atomlist.dat',status='unknown')
c$$$           DO I=1,NP
c$$$                WRITE(89,*) I,IMOLi(I),MOL(I)
c$$$           ENDDO 
c$$$           CLOSE(88)
c$$$           CLOSE(89)
! end debug
C
C Loop over uncounted atoms and find any H2's
      DO I = 1,NP
       IF ( .NOT.MOL(I) ) THEN
         MOLCNTi=MOLCNTi+1
         CNT=CNT+1
         IF (CNT.GT.NP) THEN
           WRITE(13,*)'When addin H2s'
           WRITE(13,*)'Change dimension of MOLSTi in molcheck.f'
           WRITE(13,*) ' Atom',I,CNT,R0i(:,I)
           WRITE(13,*)'but something else is probably wrong'
           OPEN(UNIT=87,file='mollist.err',status='unknown')
           DO NN=1,MOLCNTi
             MI = MPNTi(NN)
             MF = MPNTi(NN+1)-1
             DO MJ=MI,MF
                II = MOLSTi(MJ)
                WRITE(87,*) NN,II,MOL(II)
                IF( I.EQ.II) THEN
                  WRITE(13,*) 'Atom ',I,' was already in mol',NN
                ENDIF
             ENDDO
           ENDDO 
           CLOSE(87)
!          Check for repeats
           MPNTo(:) = 0 
           DO J =1,NP
             DO NN=1,MOLCNTi
               MI = MPNTi(NN)
               MF = MPNTi(NN+1)-1
               DO MJ=MI,MF
                 II = MOLSTi(MJ)
                 IF( J.EQ.II) THEN
                   IF ( MPNTo(J).NE.0 ) THEN
                   WRITE(13,*) 'Atom ',J,' was already in mol',MPNTo(J)
                   ELSE
                     MPNTo(J) = NN
                  ENDIF
                 ENDIF
               ENDDO
             ENDDO 
           ENDDO
           STOP
         ENDIF
         MOLSTi(CNT) = I
         MPNTi(MOLCNTi) = CNT
         MOL(I) =.TRUE.
         IMOLi(I) = MOLCNTi
         JBEGIN=NABORSi(I)
         JEND=NABORSi(I+1)-1
         DO N = JBEGIN,JEND
             B=LISTi(N)
             IF( KTYPEi(B).NE.ihyd) THEN
              WRITE(13,*) 'Molecule was missed intial loop'
              WRITE(13,*) 'Atom:',B,KTYPEi(B),R0i(:,B),MOL(B)
              WRITE(13,*) 'has atoms:',LISTi(JBEGIN:JEND)
              STOP
             ELSEIF( .NOT. MOL(B) ) THEN
              CNT = CNT + 1
              MOLSTi(CNT) = B
              MOL(B) =.TRUE.
              IMOLi(B) = MOLCNTi
             ENDIF
         ENDDO
       ENDIF        
      ENDDO
!      
C  Set pointer for end of last molecule as end of list 
      MPNTi(MOLCNTi+1) = CNT + 1
C
C Set molecules with greater than threshold number ofatoms as substrate 
      DO N=1,MOLCNTi
        MI = MPNTi(N)
        MF = MPNTi(N+1)-1
        C = MF - MI
        IF (C.GE.STHS) THEN
          DO MJ=MI,MF
            I = MOLSTi(MJ)
            ASUBi(I) = .TRUE.
!            IMOLi(I) = 0
          ENDDO
        ENDIF
      ENDDO
      DEALLOCATE( MOL ) 
      ENDIF
C
      END SUBROUTINE molcheck
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Determine molecules that have formed on the surface
!
      SUBROUTINE smolcheck
!
      USE ANALYSIS
      USE STRUCTR
      USE POTS
      USE MDSTP
      USE MPIvars
!
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
      INTEGER :: I,J,MI,MF,MJ,MCK,CK,A,B,C,N
     &   ,JBEGIN,JEND,CREF ,CNT,RPCNT
      REAL*8 :: RIJSQ,DR,RCUTSQ
      LOGICAL :: RPTCHK,VERBOSE,TERM,BTOBEAM
      LOGICAL, ALLOCATABLE :: MOL(:)
!
      IF ( mynode.EQ.0 )THEN
      ALLOCATE( MOL(NP) ) 
C Determine surface molecules
      MOL(:) = .TRUE.
      SMCNTi = 0
      CNT = 0
      DO I=1,NP
       IF( ASUBi(I).AND.MOL(I) ) THEN
         IF( .NOT.ASUBo(I).OR.IMOLo(I).NE.0 ) THEN
          SMCNTi = SMCNTi + 1
          CNT = CNT  + 1
! begin debug
       WRITE(88,*) 'Atom ',I,' initial atom in',SMCNTi
       WRITE(88,*) ASUBi(I),ASUBo(I),IMOLo(I)
! end debug
          SMOLSTi(CNT) = I 
          SMPNTi(SMCNTi) = CNT
          MOL(I) = .FALSE.
          SIMOLi(I) = -SMCNTi
          TERM=.TRUE.
          DO WHILE (TERM)
            MI = SMPNTi(SMCNTi)
            MF = CNT
            CREF = CNT
            DO MJ = MI,MF
              A = SMOLSTi(MJ)
              JBEGIN=NABORSi(A)
              JEND=NABORSi(A+1)-1
              DO 1001 N = JBEGIN,JEND
                 B = LISTi(N)
! begin debug
       WRITE(88,*) 'Possible: ',B,ASUBo(B),RECTi(B)
! end debug
                 IF(ASUBo(B).AND.RECTi(B).EQ.0 ) GOTO 1001
                 DO MCK = MI,CNT
                  CK = SMOLSTi(MCK)
                  IF ( B.eq.CK.OR.KTYPEi(B).EQ.ihyd) GOTO 1001
                 ENDDO
                 CNT = CNT + 1
                 SMOLSTi(CNT) = B
                 SIMOLi(B) = -SMCNTi
                 MOL(B) =.FALSE.
! begin debug
       WRITE(88,*) 'Atom ',B,' was added to list'
! end debug
 1001         CONTINUE
            ENDDO
            IF( CNT.EQ.CREF) TERM=.FALSE.
          ENDDO 
C Add hydrogens to list
          MI = SMPNTi(SMCNTi)
          MF = CNT
          DO  MJ=MI,MF
           A = SMOLSTi(MJ)
           JBEGIN=NABORSi(A)
           JEND=NABORSi(A+1)-1
           DO 2001 N = JBEGIN,JEND
             B=LISTi(N)
             IF( MOL(B).AND.KTYPEi(B).EQ.ihyd ) THEN
               CNT = CNT + 1
               SIMOLi(B) = -SMCNTi
               SMOLSTi(CNT) = B
               MOL(B) =.FALSE.     
             ENDIF 
 2001      CONTINUE
          ENDDO
         ENDIF 
        ENDIF 
      ENDDO !NP
C  Set pointer for end of last molecule as end of list 
      SMPNTi(SMCNTi+1) = CNT + 1
! begin debug
!      WRITE(13,*) SMCNTi,' surface species found',CNT
!      DO N=1,SMCNTi
!        Mi = SMPNTi(N)
!     WRITE(13,*) N,MI,SMOLSTi(Mi)
!      ENDDO
! end debug
      DEALLOCATE( MOL ) 
!
      ENDIF
      RETURN
      END SUBROUTINE smolcheck 
