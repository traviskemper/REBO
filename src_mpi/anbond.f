C Bond analysis for PS  
!
      SUBROUTINE ANBOND
!
      USE MPIvars
      USE ANALYSIS
      USE STRUCTR
!
      IMPLICIT none
!
      INTEGER :: I,HCNT,ALAB(NPMAX)
!     
      IF ( mynode.EQ.0 ) THEN
      DO I=1,NP
       HCNT=NCTi(I) - CNNi(I)
       IF ( CNNi(I).eq.2 .and. HCNT.eq.2 ) THEN
         ALABi(I)=1
       ELSEIF (  CNNi(I).eq.3 .and. HCNT.eq.1 ) THEN
         ALABi(I)=2
       ELSEIF (  CNNi(I).eq.3 .and. HCNT.eq.0 ) THEN
         ALABi(I)=3
       ELSEIF (  CNNi(I).eq.2 .and. HCNT.eq.1 ) THEN
         ALABi(I)=4
       ENDIF
      ENDDO
      ENDIF
!
      ENDSUBROUTINE ANBOND 
