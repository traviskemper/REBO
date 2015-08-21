C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Programer:
C     Travis Kemper
C     Department of Materials SCience and ENgineering
C     University of FLorida
C     traviskemper@ufl.edu
C
C     Version 1.0 06/01/11 T. W. Kemper
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Find number hydrogens and carbons connected to each atom
!
      SUBROUTINE cutfunc 
!
      USE STRUCTR
      USE POTS
      USE MDSTP
!
      IMPLICIT none
!
      INTEGER :: I,J,KJ,KI,JBEGIN,JEND,NBJ
      REAL*8 :: RC
!
      DO 500 I=1,NP
           KI =KTYPE(I)
           JBEGIN=NABORS(I)
           JEND=NABORS(I+1)-1
           XHC(I,1)=1
           XHC(I,2)=1
           IF(JBEGIN.GT.JEND) GO TO 500
           DO 490 NBJ=JBEGIN,JEND
                IF(LCHECK(NBJ).ne.1) GO TO 490
                J=LIST(NBJ)
                KJ =KTYPE(J)
                RC=RCOR(NBJ)
                XHC(I,KTYPE(J))=XHC(I,KTYPE(J))+WW(NBJ)
490        CONTINUE
C
500   CONTINUE
      RETURN
      END SUBROUTINE cutfunc 
