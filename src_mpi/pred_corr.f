      subroutine setpc

      USE MPIvars
      USE SPECIF
      USE STRUCTR

      IMPLICIT none 

C
C NORDSIECK PREDICTOR-CORRECTOR COEFFICIENTS
C
      F02=1.0d0/6.0d0
      F12=5.0d0/6.0d0
      F32=1.0d0/3.0d0
c
      return 
      end 
c 
      SUBROUTINE CPRED


      USE SPECIF
      USE STRUCTR

      IMPLICIT none 

      INTEGER :: I,J,II

      IF((KFLAG.NE.6).AND.(KFLAG.NE.8)) THEN
           IF(NMA.NE.0) CALL PRED
           TTIME=TTIME+DELTA
           time = time + delta
      ENDIF
      return
      end 

      SUBROUTINE ccorr


      USE SPECIF
      USE STRUCTR

      IMPLICIT none 

      INTEGER :: I,J,II

      IF((KFLAG.NE.6).AND.(KFLAG.NE.8)) THEN
              IF(NMA.NE.0) CALL CORR
      endif
      return
      end 

      SUBROUTINE PRED

      USE dyn_array
      USE MPIvars
      USE SPECIF
      USE STRUCTR

      IMPLICIT none 

      INTEGER :: I,J,II

C
C THIRD ORDER PREDICTOR
C
C      WRITE(*,*) 'start predictor'
      DO 31 I=1,NPtot4_node
           DO 30 J=1,3
C                I = mlist(II)
              IF (itr_node(I).LT.2) THEN
                R0_node(J,I)=R0_node(J,I)+R1_node(J,I)+R2_node(J,I)+
     &          R3_node(J,I)
                R1_node(J,I)=R1_node(J,I) + 2.0d0*R2_node(J,I) +
     &          3.0d0*R3_node(J,I)
                R2_node(J,I)=R2_node(J,I) + 3.0d0*R3_node(J,I)
                R0_node(J,I)=R0_node(J,I)-CUBE(J)*
     &          ANINT(R0_node(J,I)/CUBE(J))
              ENDIF

C              WRITE(*,*) I,J
30         CONTINUE
31    CONTINUE
C      IF (LSTEP.EQ.22) THEN
C      DO i=1,NP_node
C        IF (mynode.EQ.0) THEN
C          WRITE(31,999) NA(i),R0_node(1,i),R0_node(2,i),
C     &    R0_node(3,i)
C        ENDIF
C        IF (mynode.EQ.1) THEN
C          WRITE(32,999) NA(i),R0_node(1,i),R0_node(2,i),
C     &    R0_node(3,i)
C        ENDIF
C        IF (mynode.EQ.2) THEN
C          WRITE(33,999) NA(i),R0_node(1,i),R0_node(2,i),
C     &    R0_node(3,i)
C        ENDIF
C        IF (mynode.EQ.3) THEN
C          WRITE(34,999) NA(i),R0_node(1,i),R0_node(2,i),
C     &    R0_node(3,i)
C        ENDIF
C      ENDDO
C      ENDIF
C999   FORMAT(I5,3(E20.11))
C      WRITE(*,*) 'finish predictor'

      RETURN
      END

C
C*****
C
      SUBROUTINE CORR

      USE dyn_array
      USE SPECIF
      USE STRUCTR
      USE PARAMS
      USE POTS
      USE MDSTP
      USE MPIvars

      IMPLICIT none

      INTEGER ::  I,J,II,MM,ICCC
C      REAL*8 :: DE,XXM,R,RI,RSQCRT,RSQ,XMT
      REAL*8 :: XXM,DE,RI
C      DIMENSION COM(3)
c
C  THIRD ORDER CORRECTOR
C
      DE=DELTSQ/ECONV
!      WRITE(*,*) mynode,'IN CORR',NPtot4_node
C      CALL write_error
      DO I=1,NPtot4_node !II=1,NMA
        DO J=1,3
C           I = mlist(II)
          IF (itr_node(I).LT.2) THEN
            XXM=XMASS(KTYPE_node(I))
            RI=R2_node(J,I) - DE*RNP(J,NA(i))/XXM
            R0_node(J,I)=R0_node(J,I) - RI*F02
            R1_node(J,I)=R1_node(J,I) - RI*F12
            R2_node(J,I)=R2_node(J,I) - RI
            R3_node(J,I)=R3_node(J,I) - RI*F32
C
            R0_node(J,I)=R0_node(J,I)-
     &      CUBE(J)*ANINT(R0_node(J,I)/CUBE(J))
          ENDIF
C
        ENDDO
      ENDDO
C      DO i=1,NPtot4_node
C        IF (mynode.EQ.0) THEN
C          WRITE(31,999) NA(i),R0_node(1,i),R0_node(2,i),
C     &    R0_node(3,i)
C        ENDIF
C        IF (mynode.EQ.1) THEN
C          WRITE(32,999) NA(i),R0_node(1,i),R0_node(2,i),
C     &    R0_node(3,i)
C        ENDIF
C        IF (mynode.EQ.2) THEN
C          WRITE(33,999) NA(i),R0_node(1,i),R0_node(2,i),
C     &    R0_node(3,i)
C        ENDIF
C        IF (mynode.EQ.3) THEN
C          WRITE(34,999) NA(i),R0_node(1,i),R0_node(2,i),
C     &    R0_node(3,i)
C        ENDIF
C      ENDDO
C999   FORMAT(I5,3(E20.11))
      RETURN
      END

      SUBROUTINE CHKNAB

      USE dyn_array
      USE MPIvars
      USE STRUCTR
      USE POTS 

      IMPLICIT none

      INTEGER :: LCHK,LCHKt,I,J,ierr
      REAL*8 :: R,RI,RSQ,RSQCRT

      INCLUDE 'mpif.h'
C
C UPDATE NEIGHBOR LIST?
C
!      WRITE(*,*) mynode,'CHKNAB starting'
      LCHK=0
      LCHKt=0
      RSQCRT=(0.5d0*RLL)**2
      DO 4 I=1,NP_node
        RSQ=0.0D0
        DO 3 J=1,3
           R=R0_node(J,I)-R0L(J,I)
           R = R - CUBE(J)*ANINT(R/CUBE(J))
           RSQ=RSQ+R*R
C
3       CONTINUE
        IF(RSQ.GT.RSQCRT) THEN
!          WRITE(*,*) mynode,NA(i)
          LCHKt=1
        ENDIF
4     CONTINUE
!      WRITE(*,*) mynode,'LCHKt=',LCHKt
      CALL MPI_REDUCE(LCHKt,LCHK,1,MPI_INTEGER,MPI_SUM,0,
     &MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(LCHK,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
C     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      WRITE(*,*) mynode,'CHKNAB MPI finished'
      IF (LCHK.GE.1)  GO TO 5
      LCHK=0
      GO TO 8
C

5     CONTINUE
      LCHK=1
C
C UPDATE NEIGHBOR LIST
C
C      DO 7 I=1,NP_node
C           DO 6 J=1,3
C                R0L(J,I)=R0_node(J,I)
C6          CONTINUE
C7     CONTINUE
C
C      KCHK=KCHK+1
C      NCHK(KCHK)=KB
C      LCHK=1

C
8     CONTINUE
C
C      IF(NP.EQ.NMA) RETURN
C      WRITE(*,*) mynode,'in CHKNAB - LCHK=',LCHK
      RETURN
      END

      SUBROUTINE MOVCTR
C
C***IMPORTANT********************************
C                                           *
C IF NO RIGID ATOMS, SUBTRACT COM VELOCITY  *
C                                           *
C REMOVE THIS CODE FOR MOLECULAR COLLISION  *
C                                           *
C********************************************
C
      USE dyn_array
      USE MPIvars
      USE POTS

      IMPLICIT none 
c
      INCLUDE 'mpif.h'
C
      INTEGER :: ICCC,MM,I,ierr
      REAL*8 :: COM(3),COMtmp(3),XMT,XMTemp

      WRITE(*,*) mynode,'ICCC=',ICCC
      ICCC=2
      IF(ICCC.EQ.1) THEN
      XMT=0.0D0
C
      DO 24 I=1,NP_node
           XMT=XMT+XMASS(KTYPE_node(I))
           XMTemp=XMT
24    CONTINUE
      XMT=0.0d0
      CALL MPI_REDUCE(XMTemp,XMT,1,MPI_REAL8,MPI_SUM,0,
     &MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(XMT,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
C
      COM(MM)=0.0d0
      DO 27 I=1,NP_node
        DO 25 MM=1,3
          COM(MM)=COM(MM)+R1_node(MM,I)*XMASS(KTYPE_node(I))
          COMtmp(MM)=COM(MM)
25      CONTINUE
27    CONTINUE
      CALL MPI_REDUCE(COMtmp(MM),COM(MM),3,MPI_REAL8,MPI_SUM,0,
     &MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(COM(MM),3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
c
      DO MM=1,3
        COM(MM)=COM(MM)/XMT
      ENDDO
c
      DO 26 I=1,NPtot4_node
        DO MM=1,3
          R1_node(MM,I)=R1_node(MM,I)-COM(MM)
        ENDDO
26    CONTINUE

      ENDIF
C
      RETURN
      END

