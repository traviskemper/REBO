C Calculates the center of mass
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  CNT - Number of atoms in ELIST
C  ELIST - list of atoms in molecule
C  ELCNT - Number of each element type
C  COM   - Center of mass
C  MOLM  - Molecular mass in AMU
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  Travis Kemper - MSE UF Dr. Sinnott

      SUBROUTINE COMAS(CNT,COM,MOLM,MVEL)

      USE ANALYSIS
      USE POTS
      USE STRUCTR
      USE SPECIF
      USE MPIvars
!
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
      INTEGER :: CNT,C,EN,AN,N,I,D
      REAL*8 :: EM,MOLM,MMA(3),COM(3),MVEL(3)

      IF ( mynode.EQ.0 )THEN
        IF ( CNT.GT.NP ) THEN
          WRITE(*,*) 'For list:'
          WRITE(*,*) ELISTi
          WRITE(*,*) 'Number of atoms passed to comas  is greater
     & than the number of atoms in system!',NP
         STOP
        ENDIF

      MOLM=0.0D0
      MMA(:)=(/0.0D0,0.0D0,0.0D0/)
      MVEL(:)= 0.0d0 

      DO C=1,CNT
        I=ELISTi(C)      !Passed list of atom numbers
        EN = KTYPEi(I)   !Element number internal
C        AN = ELTRN(EN)  !Internal atom # 
C        ELCNT(AN)= ELCNT(AN)+1
        EM  = XMASS(EN)
        MOLM = MOLM + EM
C        WRITE(*,*) AN,EN,EM,MOLM
        DO D=1,3   !average mass possition
          MMA(D) = MMA(D) + R0i(D,I)*EM
          MVEL(D) = MVEL(D) + R1i(D,I)/DELTA 
        ENDDO
      ENDDO
!
      DO D=1,3
        COM(D) = MMA(D)/MOLM
        MVEL(D) = MVEL(D)/FLOAT(CNT)
      ENDDO
!
      ENDIF
      RETURN
!
      END SUBROUTINE COMAS
