      subroutine setgle

      USE SPECIF
      USE PARAMS
      USE STRUCTR

      IMPLICIT none

C Debye temperature converted to fs-1 X 2Pi
C
      WD=2230.0d0*(2.08365D-05)*2.0d0*PI

      nlr = 3*NTA/2

C not sure about this allocation may cause problems
C change nlr*2 to 3*NPMAX from AIREBO
C -travisk

c
C LANGEVIN PARAMETERS
C
      TR=TEM/EPSI/ECONV
      BET=WD*PI*ECONV/6.0d0/DELTA
      GSIG=SQRT(2.0d0*TR*ECONV*BET)
      PI2=PI*2.0d0
      DNLA=3.0d0*DELTA*DELTA*FLOAT(NTA)
      return
      end 

      SUBROUTINE thermos

      USE SPECIF

      IMPLICIT none
c
      IF(KFLAG.EQ.-1) CALL GLEQ
      IF(KFLAG.EQ.1) CALL BERE
      IF(KFLAG.EQ.2) CALL ZERO
      IF(KFLAG.EQ.3) CALL HOOV
      return
      end 


      SUBROUTINE GLEQ
C
C FRICTION AND RANDOM FORCE
C
      USE SPECIF
      USE STRUCTR 
      USE POTS
      USE MDSTP

      IMPLICIT none
C
      INTEGER :: I,II,J
      REAL*8 :: RR,PRE,BM,SM,RRZP,RANNUM
C
      nlr = 3*NTA/2
C
      DO 20 I=1,NLR
           RR = RANNUM(I)
           IF(RR.LT.1.0D-06) GO TO 20
           PRE=SQRT(-2.0d0*LOG(RR))
           GL(I)=PRE*COS(PI2*RANNUM(I))*GSIG
           GL(I+NLR)=PRE*SIN(PI2*RANNUM(I))*GSIG
C
   20 CONTINUE
C
      DO 30 II=1,NTA
           I=NLIST(II)
           BM=BET*XMASS(KTYPE(I))
           SM=SQRT(XMASS(KTYPE(I)))
           DO 29 J=1,3
                rrzp = RNP(I,J)
                RNP(I,J)=RNP(I,J)-BM*R1(I,J)-SM*GL(II+(J-1)*NTA)
29         CONTINUE
30    CONTINUE
      RETURN
      END

C
      SUBROUTINE BERE
C
C USE BERENDSEN SCHEME
C

      USE SPECIF
      USE STRUCTR
      USE POTS
      USE MDSTP

      IMPLICIT none

      INTEGER :: J,II,I
      REAL*8 :: XX,SC,SM

      XX=0.0d0
      DO 11 J=1,3
C
           DO 10 II=1,NTA
                I = NLIST(II)
                XX=XX+(R1(I,J)*R1(I,J))*XMASS(KTYPE(I))
   10      CONTINUE
C
   11 CONTINUE
c
      if(xx.lt.0.0d-7) then
            write(*,*) 'T=0, Reset Thermostat to other than 1'
            stop
      endif  
C
      IF(KFLAG.EQ.1) THEN
C
           SC=BET*(TR*DNLA/XX-1.0d0)
C
           DO 30 II=1,NTA
                I = NLIST(II)
                SM=XMASS(KTYPE(I))*SC
C
                DO 29 J=1,3
                     RNP(I,J)=RNP(I,J)+SM*R1(I,J)
29              CONTINUE
C
30         CONTINUE
C
      ELSE
C
           SC=SQRT(TR*6.0D0*DELTSQ*FLOAT(NTA)/XX)
           DO 32 II=1,NTA
                I = NLIST(II)
                DO 31 J=1,3
                     R1(I,J)=SC*R1(I,J)
31              CONTINUE
C
32         CONTINUE
C
      ENDIF
      RETURN
      END
C
      SUBROUTINE HOOV
C
C USE EVANS-HOOVER SCHEME
C
      USE SPECIF
      USE POTS
      USE CONTM
      USE STRUCTR
      USE MDSTP

      IMPLICIT none

      INTEGER :: J,I
      REAL*8 :: FF,DF,SC

C
C this used for all atoms
c
      FF=0.0D0
      DF=0.0D0
      DO 2 J=1,3
           DO 1  I=1,NP
                FF=FF+RNP(I,J)*R1(I,J)
                DF=DF+R1(I,J)*R1(I,J)*XMASS(KTYPE(I))
1          CONTINUE
2     CONTINUE
C
      SC=FF/DF
C
      DO 4 J=1,3
C
           DO 3 I=1,NP
                RNP(I,J)=RNP(I,J)-SC*R1(I,J)*XMASS(KTYPE(I))
3          CONTINUE
C
4     CONTINUE
      RETURN
      END
C
      SUBROUTINE ZERO
C
C ZERO VELOCITIES
C
      USE SPECIF
      USE STRUCTR
      USE MDSTP

      IMPLICIT none

      INTEGER :: I,J
C
C note - this works on all atoms
c
      DO 31 J=1,3
           DO 30 I=1,NP
                R1(I,J)=0.0d0
30         CONTINUE
31    CONTINUE
      RETURN
      END
