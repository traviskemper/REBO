C
      SUBROUTINE TCUINT(I,JN,KL,KI,XX1,XX2,XX3,NH,NC,NF
     &   ,ANSY,ANSY1,ANSY2,ANSY3)
C
C Bicubic spline
C
      USE MPIvars
      USE POTS
     
      IMPLICIT REAL*8(A-H,O-Z)
c
      INCLUDE 'mpif.h'
c
      ANSY=0.0d0
      ANSY1=0.0d0
      ANSY2=0.0d0
      ANSY3=0.0d0
C
C      write(6,*) 'nf=',nf
      IF((KI.EQ.0).OR.(NH.EQ.0).OR.(NC.EQ.0)
     &   .or.(NF.EQ.0)) then
           WRITE(*,*) mynode,'tcunit',I,JN,KI,KL,NH,NC,NF
           include 'close.inc'
           CALL write_error
      ENDIF
C
      DO 32 J=1,64
        X=CLM(KI,NH,NC,NF,J)*
     &   (XX1**IN3(J,1))*(XX2**IN3(J,2))*(XX3**IN3(J,3))
           ANSY=ANSY+X
           ANSY1=ANSY1+X*IN3(J,1)/XX1
           ANSY2=ANSY2+X*IN3(J,2)/XX2
           ANSY3=ANSY3+X*IN3(J,3)/XX3
32    CONTINUE
C
      RETURN
      END

