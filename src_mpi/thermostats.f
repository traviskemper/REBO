      subroutine setgle

      USE MPIvars
      USE SPECIF
      USE PARAMS
      USE STRUCTR

      IMPLICIT none
c
C Debye temperature converted to fs-1 X 2Pi
C
      WD=2230.0d0*(2.08365D-05)*2.0d0*PI

      nlr = 3*NTA/2


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

C      WRITE(*,*) mynode,'KFLAG=',KFLAG
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
      USE dyn_array
      USE MPIvars
      USE SPECIF
      USE STRUCTR 
      USE POTS
      USE MDSTP

      IMPLICIT none

      INTEGER :: I,II,J
      REAL*8 :: RR,PRE,BM,SM,RRZP,RANNUM
c
      NLR = 3*NTA/2

      DO 20 I=1,NLR
           RR = RANNUM(I)
           IF(RR.LT.1.0D-06) GO TO 20
           PRE=SQRT(-2.0d0*LOG(RR))
           GL(I)=PRE*COS(PI2*RANNUM(I))*GSIG
           GL(I+NLR)=PRE*SIN(PI2*RANNUM(I))*GSIG
   20 CONTINUE
C
C      DO i=1,NP_node
C        IF (mynode.EQ.0) THEN
C          WRITE(31,999) NA(i),RNP(1,NA(i)),RNP(2,NA(i)),
C     &    RNP(3,NA(i))
C        ENDIF
C        IF (mynode.EQ.1) THEN
C          WRITE(32,999) NA(i),RNP(1,NA(i)),RNP(2,NA(i)),
C     &    RNP(3,NA(i))
C        ENDIF
C        IF (mynode.EQ.2) THEN
C          WRITE(33,999) NA(i),RNP(1,NA(i)),RNP(2,NA(i)),
C     &    RNP(3,NA(i))
C        ENDIF
C        IF (mynode.EQ.3) THEN
C          WRITE(34,999) NA(i),RNP(1,NA(i)),RNP(2,NA(i)),
C     &    RNP(3,NA(i))
C        ENDIF
C      ENDDO
      DO 30 I=1,NPtot4_node !II=1,NTA
        IF (itr_node(I).EQ.1) THEN
C           II=NLIST(NA(I))
           BM=BET*XMASS(KTYPE_node(I))
           SM=SQRT(XMASS(KTYPE_node(I)))
           DO 29 J=1,3
                IF (disp(j).NE.0.0d0) GOTO 29
                rrzp = RNP(J,NA(i))
                RNP(J,NA(i))=RNP(J,NA(i))-BM*R1_node(J,I)-
     &          SM*GL(nlist(NA(i))+(J-1)*NTA)
29         CONTINUE
C        IF (I.GT.NP_node) THEN
C        IF (mynode.EQ.0) THEN
C          WRITE(33,999) NA(i),RNP(1,NA(i)),RNP(2,NA(i)),
C     &    RNP(3,NA(i))
C        ENDIF
C        IF (mynode.EQ.1) THEN
C          WRITE(34,999) NA(i),RNP(1,NA(i)),RNP(2,NA(i)),
C     &    RNP(3,NA(i))
C        ENDIF
C        ENDIF
        ENDIF
30    CONTINUE
C      DO i=1,NP_node
C        IF (mynode.EQ.0) THEN
C          WRITE(31,999) NA(i),RNP(1,NA(i)),RNP(2,NA(i)),
C     &    RNP(3,NA(i))
C        ENDIF
C        IF (mynode.EQ.1) THEN
C          WRITE(32,999) NA(i),RNP(1,NA(i)),RNP(2,NA(i)),
C     &    RNP(3,NA(i))
C        ENDIF
C        IF (mynode.EQ.2) THEN
C          WRITE(33,999) NA(i),RNP(1,NA(i)),RNP(2,NA(i)),
C     &    RNP(3,NA(i))
C        ENDIF
C        IF (mynode.EQ.3) THEN
C          WRITE(34,999) NA(i),RNP(1,NA(i)),RNP(2,NA(i)),
C     &    RNP(3,NA(i))
C        ENDIF
C      ENDDO
C999   FORMAT(I5,3(E20.11))

C      WRITE(*,*)mynode,'finish thermo'
      RETURN
      END

C
      SUBROUTINE BERE
C
C USE BERENDSEN SCHEME
C
      USE dyn_array
      USE MPIvars
      USE SPECIF
      USE STRUCTR
      USE POTS
      USE MDSTP

      IMPLICIT none

      INTEGER :: J,II,I
      REAL*8 :: XX,SC,SM

      INCLUDE 'mpif.h'
c
      XX=0.0d0
      DO 11 J=1,3
C
           DO 10 I=1,NPtot4_node !II=1,NTA
             IF (itr_node(I).EQ.1) THEN
C                I = NLIST(II)
                XX=XX+(R1_node(J,I)*R1_node(J,I))*XMASS(KTYPE_node(I))
             ENDIF
   10      CONTINUE
C
   11 CONTINUE

      if(xx.lt.0.0d-7) then
            write(*,*) 'T=0, Reset Thermostat to other than 1'
            CALL write_error
      endif  
C
      IF(KFLAG.EQ.1) THEN
C
           SC=BET*(TR*DNLA/XX-1.0d0)
C
           DO 30 I=1,NPtot4_node  !II=1,NTA
             IF (itr_node(I).EQ.1) THEN
C                I = NLIST(II)
                SM=XMASS(KTYPE_node(I))*SC
C
                DO 29 J=1,3
                     RNP(J,NA(i))=RNP(J,NA(i))+SM*R1_node(J,I)
29              CONTINUE
             ENDIF
30         CONTINUE
C
      ELSE
C
           SC=SQRT(TR*6.0D0*DELTSQ*FLOAT(NTA)/XX)
           DO 32 I=1,NPtot4_node !II=1,NTA
             IF (itr_node(I).EQ.1) THEN
C                I = NLIST(II)
C
                DO 31 J=1,3
                     R1_node(J,I)=SC*R1_node(J,I)
31              CONTINUE
             ENDIF
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
      USE dyn_array
      USE MPIvars
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
           DO 1  I=1,NPtot4_node
                FF=FF+RNP(J,NA(i))*R1_node(J,I)
                DF=DF+R1_node(J,I)*R1_node(J,I)*XMASS(KTYPE_node(I))
1          CONTINUE
2     CONTINUE
C
      SC=FF/DF
C
      DO 4 I=1,NPtot4_node
           DO 3 J=1,3
                RNP(J,NA(i))=RNP(J,NA(i))-SC*R1_node(J,I)*
     &          XMASS(KTYPE_node(I))
3          CONTINUE
4     CONTINUE
      RETURN
      END
C
      SUBROUTINE ZERO
C
C ZERO VELOCITIES
C
      USE dyn_array
      USE MPIvars
      USE SPECIF
      USE STRUCTR
      USE MDSTP

      IMPLICIT none

      INTEGER :: I,J
C
C note - this works on all atoms
c
      DO 31 J=1,3
           DO 30 I=1,NPtot4_node
                R1_node(J,I)=0.0d0
30         CONTINUE
31    CONTINUE
      RETURN
      END
