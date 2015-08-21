
      SUBROUTINE RADIC(KI,KJ,XNT1,XNT2,CONJUG,RAD,DRDL,DRDM,DRDN)

      USE POTS

      IMPLICIT none
 
      INTEGER :: KJ,KL,KI,L,M,N,KIKJ,J
      REAL*8 :: XNT1,XNT2,CONJUG,RAD,DRDL,DRDM,DRDN,X

C TRICUBIC SPLINE
C
C CONJUG=1: NONCONJUGATED
C       >1: CONJUGATED
C
C Account for posible round off errors
        IF ( XNT1.GT.4.999d0) THEN
          XNT1 = 4.999d0
        ELSE IF ( XNT1.LT.1.0d0) THEN
          XNT1 = 1.0d0
        ENDIF
        IF ( XNT2.GT.4.999d0) THEN
          XNT2 = 4.999d0
        ELSE IF ( XNT2.LT.1.0d0) THEN
          XNT2 = 1.0d0
        ENDIF

        IF (CONJUG.GT.9.999d0) THEN
           CONJUG = 9.999d0
        ELSE IF (CONJUG.LT.1.0d0) THEN
          CONJUG = 1.0d0
        ENDIF

      L=INT(XNT1+1.0D-12)
      M=INT(XNT2+1.0D-12)
      N=INT(CONJUG)
      IF(L*M*N.eq.0) THEN
           WRITE(6,*) "atom",KI,"bonded to atom",KJ
           WRITE(6,*) 'Zero appear',l,m,n
           include 'close.inc'
           STOP
      endif
      RAD=0.0d0
      DRDL=0.0d0
      DRDM=0.0d0
      DRDN=0.0d0
      KIKJ = KI + KJ - 1
      if ((ki.eq.iflor).and.(kj.ne.iflor)) kikj=kikj-2
      if ((ki.eq.iflor).and.(kj.eq.iflor)) kikj=kikj-4
      if ((kj.eq.iflor).and.(ki.ne.iflor)) kikj=kikj-2
C
      IF(L.GE.4) THEN
           L=4
           XNT1=4.0D0
      ENDIF
C
      IF(M.GE.4) THEN
           M=4
           XNT2=4.0D0
      ENDIF
C
      IF(N.GE.9) THEN
           N=9
           CONJUG=9.0D0
      ENDIF
C
      DO 32 J=1,64
           X=CLMN(KIKJ,L,M,N,J)*
     &       (XNT1**IN3(J,1))*(XNT2**IN3(J,2))*(CONJUG**IN3(J,3))
           RAD=RAD+X
           DRDL=DRDL+X*IN3(J,1)/XNT1
           DRDM=DRDM+X*IN3(J,2)/XNT2
           DRDN=DRDN+X*IN3(J,3)/CONJUG
32    CONTINUE
      RETURN
      END
C

