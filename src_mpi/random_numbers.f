      subroutine setran

      USE SPECIF

      IMPLICIT none

      REAL*8 RANNUM 

      CALL SETRN(PSEED)
      PSEED=RANNUM(1)

      return
      end 

      SUBROUTINE SETRN(Q)
      USE SPECIF

      IMPLICIT none

      REAL*8 :: Q,QC
c
C***  INITIALIZE WITH A CALL TO SETRN(0.D0-1.D0)
      QA1=2057713.0d0
      QA2=16676923.0d0
      QBASE=2.d0**24
      QC=DINT(QBASE*(QBASE*Q))
      QB1=DINT(QC/QBASE)
      QB2=QC-QB1*QBASE
      QB1=DMOD(QB1,QBASE)
      QB2=DINT(QB2/2.0d0) * 2.0d0 + 1.0d0
      RETURN
      END
C
C*************
C
      FUNCTION RANNUM(I)


      USE SPECIF

      IMPLICIT none

      REAL*8 :: I,QD2,QE2,QC2,RANNUM
C
C***  FROM CLAMPS AT NRCC - FROM KALOS
      QD2=QA2*QB2
      QE2=DINT(QD2/QBASE)
      QC2=QD2-QBASE*QE2
      QB1=DMOD(QE2+DMOD(QA1*QB2,QBASE)+DMOD(QA2*QB1,QBASE),QBASE)
      QB2=QC2
      RANNUM=QB1/QBASE
      RETURN
      END

