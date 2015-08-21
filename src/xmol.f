      SUBROUTINE XMOL

      USE POTS
      USE MDSTP
      USE STRUCTR

      IMPLICIT none

      INTEGER :: idum,N,I
      INTEGER :: NRA,NLA

      idum=3
      NRA=0
      NLA =0

      WRITE(1,100) NP,IDUM,NRA,NLA
      WRITE(1,1800) CTITLE
C
      DO 11 I=1,NP
           WRITE(1,350) KT2(KTYPE(I)),(R0(I,N),N=1,3)
11    CONTINUE
C
      call flush(1)
  100 FORMAT(4I6)
  300 FORMAT(3E20.11)
  350 FORMAT(I5,3E20.11)
  360 FORMAT(I5,3E20.11)
 1800 FORMAT(20A2)
      return
      end

