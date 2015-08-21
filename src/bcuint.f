C
      SUBROUTINE BCUINT(KL,KI,XX1,XX2,NH,NC,ANSY,ANSY1,ANSY2)
C
C Bicubic spline
C
      USE POTS

      IMPLICIT none
 
      INTEGER :: KL,KI,NH,NC,I,J,IC
      REAL*8 :: X,XX1,XX2,ANSY,ANSY1,ANSY2

c
      ANSY=0.0d0
      ANSY1=0.0d0
      ANSY2=0.0d0


C
      IF((KI.EQ.0).OR.(NH.EQ.0).OR.(NC.EQ.0)) THEN
           WRITE(6,*) KI,KL,NH,NC
           include 'close.inc'
           STOP
      ENDIF
C
      DO J=1,16
           X=CLM(KI,NH,NC,J)*
     &       (XX1**IN2(J,1))*(XX2**IN2(J,2))
           ANSY=ANSY+X
           ANSY1=ANSY1+X*IN2(J,1)/XX1
           ANSY2=ANSY2+X*IN2(J,2)/XX2
      ENDDO
C
      RETURN

      END

