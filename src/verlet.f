      SUBROUTINE VVa

      USE SPECIF
      USE STRUCTR
      USE PARAMS
      USE POTS
      USE MDSTP

      IMPLICIT none 

      INTEGER :: I,M,II
      REAL*8 :: RV,RA,XXM

           TTIME=TTIME+DELTA
           time = time + delta

C Step 1 of velocity verlet
C  Velocity half step
      DO M=1,3
        DO I=1,NP
           XXM = XMASS(KTYPE(I))
           R1(I,M) = R1(I,M) + R2(I,M)   !Half step v(r)
           R0(I,M) = R0(I,M) + R1(I,M)
           R0(I,M)=R0(I,M)-CUBE(M)*ANINT(R0(I,M)/CUBE(M))
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE VVa

      SUBROUTINE VVb

      USE SPECIF
      USE STRUCTR
      USE PARAMS
      USE POTS
      USE MDSTP

      IMPLICIT none 

      INTEGER :: I,M,II,MM,ICCC,J
      REAL*8 :: RA,XXM,COM(3),DE,R,RI,RSQCRT,RSQ,XMT


C Step 2 of velocity verlet
      DO M=1,3
        DO I=1,NP
           XXM = XMASS(KTYPE(I))
           RA = RNP(I,M)/XXM/ECONV*DELTSQ
           R1(I,M) = R1(I,M)  + RA
           R3(I,M) =(RA - R2(I,M) )/DELTA
           R2(I,M) = RA
        ENDDO
      ENDDO

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
C      IF(NP.EQ.NMA) RETURN
C***IMPORTANT********************************
C                                           *
C IF NO RIGID ATOMS, SUBTRACT COM VELOCITY  *
C                                           *
C REMOVE THIS CODE FOR MOLECULAR COLLISION  *
C                                           *
C********************************************
C
C
      ICCC=1
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



      RETURN

      RETURN
      END SUBROUTINE VVb


