C
      SUBROUTINE OXYBCUINT(KI,KJ,XX1,XX2,NTCH,NTO,ANSY,ANSY1,ANSY2)
C
C Bicubic spline 
C   NTCH- number of Carbons+ number of hydrogens, NTO-number of Oxygen
C

      USE POTS


      IMPLICIT REAL*8(A-H,O-Z)
c

c
      ANSY=0.0d0
      ANSY1=0.0d0
      ANSY2=0.0d0
C
      IF((KJ.EQ.0).OR.(NTCH.EQ.0).OR.(NTO.EQ.0)) THEN
          WRITE(*,*) "Warning, something is Zero in oxybcuint",
     +    KI,KJ,NTCH,NTO
C          call print_conf
C          include 'close.inc'
          STOP
      ENDIF
C
!dong, test 
!      write(*,*) CLMOX
!      stop
      DO 32 J=1,16
           X=CLMOX(KI,KJ,NTCH,NTO,J)*
     &       (XX1**IN2(J,1))*(XX2**IN2(J,2))
           ANSY=ANSY+X
           ANSY1=ANSY1+X*IN2(J,1)/XX1
           ANSY2=ANSY2+X*IN2(J,2)/XX2
32    CONTINUE
C
      RETURN
      END
C


