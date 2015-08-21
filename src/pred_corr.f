      subroutine setpc

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
      USE MDSTP

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

      USE SPECIF
      USE STRUCTR

      IMPLICIT none 

      INTEGER :: I,J,II
C
C THIRD ORDER PREDICTOR
C
      DO 31 J=1,3
           DO 30 II=1,NMA
                I = mlist(II)
                R0(I,J)=R0(I,J)+R1(I,J)+R2(I,J)+R3(I,J)
                R1(I,J)=R1(I,J) + 2.0d0*R2(I,J) + 3.0d0*R3(I,J)
                R2(I,J)=R2(I,J) + 3.0d0*R3(I,J)
                R0(I,J)=R0(I,J)-CUBE(J)*ANINT(R0(I,J)/CUBE(J))
   30      CONTINUE
   31 CONTINUE
      RETURN
      END
C
C*****
C
      SUBROUTINE CORR

      USE SPECIF
      USE STRUCTR
      USE PARAMS
      USE POTS
      USE MDSTP

      IMPLICIT none

      INTEGER ::  I,J,II,MM,ICCC
      REAL*8 :: COM(3),DE,XXM,R,RI,RSQCRT,RSQ,XMT

C  THIRD ORDER CORRECTOR
C
      DE=DELTSQ/ECONV
      DO 1 J=1,3
      DO 2 II=1,NMA
           I = mlist(II)
           XXM=XMASS(KTYPE(I))
           RI=R2(I,J) - DE*RNP(I,J)/XXM
           R0(I,J)=R0(I,J) - RI*F02
           R1(I,J)=R1(I,J) - RI*F12
           R2(I,J)=R2(I,J) - RI
           R3(I,J)=R3(I,J) - RI*F32
C
           R0(I,J)=R0(I,J)-CUBE(J)*ANINT(R0(I,J)/CUBE(J))


C
   2  CONTINUE
   1  CONTINUE
C
C UPDATE NEIGHBOR LIST?
C
      RSQCRT=(0.5d0*RLL)**2
      DO 4 I=1,NP
        RSQ=0.0D0
        DO 3 J=1,3
           R=R0(I,J)-R0L(I,J)
           R = R - CUBE(J)*ANINT(R/CUBE(J))
           RSQ=RSQ+R*R
C
3       CONTINUE
        IF(RSQ.GT.RSQCRT) GO TO 5
4     CONTINUE
      LCHK=2 !Do not update
      GO TO 8 
C

5     CONTINUE

C
C UPDATE NEIGHBOR LIST
C
      DO 7 J=1,3
           DO 6 I=1,NP
                R0L(I,J)=R0(I,J)
6          CONTINUE
7     CONTINUE
C
      KCHK=KCHK+1
C      NCHK(KCHK)=KB
      LCHK=1 !Update neighbor list in caguts
C
8     CONTINUE
C
      IF(NP.EQ.NMA) RETURN
c      RETURN
C
C***IMPORTANT********************************
C                                           *
C IF NO RIGID ATOMS, SUBTRACT COM VELOCITY  *
C                                           *
C REMOVE THIS CODE FOR MOLECULAR COLLISION  *
C                                           *
C********************************************
C
C
      ICCC=2
      IF(ICCC.EQ.1) THEN
      XMT=0.0D0
C
      DO 24 I=1,NP
           XMT=XMT+XMASS(KTYPE(I))
24    CONTINUE
C
      DO 27 MM=1,3
           COM(MM)=0.0d0
c
           DO 25 I=1,NP
               COM(MM)=COM(MM)+R1(I,MM)*XMASS(KTYPE(I))
25         CONTINUE
c
           COM(MM)=COM(MM)/XMT
c
           DO 26 I=1,NP
                R1(I,MM)=R1(I,MM)-COM(MM)
26         CONTINUE
27    CONTINUE
      ENDIF
C
      RETURN
      END

