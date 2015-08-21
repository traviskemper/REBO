      subroutine overwr
      IMPLICIT REAL*8(A-H,O-Z)
      include 'common_files.inc'

      WRITE(50,1800) HEAD
      WRITE(50,100) NP,IDUM,NRA,NLA
      WRITE(50,300) TTIME,DELTA
      WRITE(50,300) (CUBE(N),N=1,3)
C
      DO 11 I=1,NP
           WRITE(50,350) I,KT2(KTYPE(I)),(R0(I,N),N=1,3)
11    CONTINUE
C
      DO 12 I=1,NP
           WRITE(50,360) I,((R1(I,N)/DELTA),N=1,3)
12    CONTINUE
C
      DO 13 I=1,NP
           WRITE(50,360) I,((R2(I,N)),N=1,3)
13    CONTINUE
C
      DO 14 I=1,NP
           WRITE(50,360) I,((R3(I,N)),N=1,3)
14    CONTINUE
C
      DO 15 I=1,NP
           WRITE(50,360) I,((R4(I,N)),N=1,3)
15    CONTINUE
      call flush(50)
      REWIND 50
      return
  100 FORMAT(4I6)
  200 FORMAT(4F12.6)
  300 FORMAT(3E20.11)
  350 FORMAT(2I5,3E20.11,I3)
  360 FORMAT(I5,3E20.11)
 1800 FORMAT(20A2)
      end

