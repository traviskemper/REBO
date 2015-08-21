C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Programer:
C     Travis Kemper
C     Department of Materials SCience and ENgineering
C     University of FLorida
C     traviskemper@ufl.edu
C
C     Version 1.0 9/25/10 T. W. Kemper
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  Set thermostates
!
      SUBROUTINE setthrmbox
!
      USE STRUCTR
      USE SPECIF
!
      IMPLICIT none
      INTEGER :: M,I,D,N
!     Find rigid  amd thermostat box sizes 
      DO D=1,3
        M=1 !min
        RIG(M,D) = RMN(M,D) + RBOX(M,D)
        IF( RBOX(M,D).LT.0.1d0) RIG(M,D)=RMN(M,D)-1.0d0 !Account for round off
        THM(M,D) = RIG(M,D) + TBOX(M,D)
        IF( TBOX(M,D).LT.0.1d0) THM(M,D)=RMN(M,D)-1.0d0 !Account for round off
        M=2 !max
        RIG(M,D) = RMN(M,D) - RBOX(M,D)
        IF( RBOX(M,D).LT.0.1d0) RIG(M,D)=RMN(M,D)+1.0d0 !Account for round off
        THM(M,D) = RIG(M,D) - TBOX(M,D)
        IF( TBOX(M,D).LT.0.1d0) THM(M,D)=RMN(M,D)+1.0d0 !Account for round off
      ENDDO
!     
      DO I=1,NP
        ITR(I) = 2
        IF(R0(I,1).GT.RIG(1,1).AND.R0(I,1).LT.RIG(2,1)
     &    .AND.R0(I,2).GT.RIG(1,2).AND.R0(I,2).LT.RIG(2,2)
     &    .AND.R0(I,3).GT.RIG(1,3).AND.R0(I,3).LT.RIG(2,3)
     &       ) ITR(I) = 1
        IF(R0(I,1).GT.THM(1,1).AND.R0(I,1).LT.THM(2,1)
     &    .AND.R0(I,2).GT.THM(1,2).AND.R0(I,2).LT.THM(2,2)
     &    .AND.R0(I,3).GT.THM(1,3).AND.R0(I,3).LT.THM(2,3)
     &       ) ITR(I) = 0
      ENDDO
!
      RETURN
      END SUBROUTINE setthrmbox
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  Set thermostates
!
      SUBROUTINE setthrmcust
!
      USE STRUCTR
      USE SPECIF
!
      IMPLICIT none
      INTEGER :: M,I,D,N
!     Find rigid  amd thermostat box sizes 
      DO D=1,3
        M=1 !min
        RIG(M,D) = RMN(M,D) + RBOX(M,D)
        IF( RBOX(M,D).LT.0.1d0) RIG(M,D)=RMN(M,D)-1.0d0 !Account for round off
        THM(M,D) = RIG(M,D) + TBOX(M,D)
        IF( TBOX(M,D).LT.0.1d0) THM(M,D)=RMN(M,D)-1.0d0 !Account for round off
        M=2 !max
        RIG(M,D) = RMN(M,D) - RBOX(M,D)
        IF( RBOX(M,D).LT.0.1d0) RIG(M,D)=RMN(M,D)+1.0d0 !Account for round off
        THM(M,D) = RIG(M,D) - TBOX(M,D)
        IF( TBOX(M,D).LT.0.1d0) THM(M,D)=RMN(M,D)+1.0d0 !Account for round off
      ENDDO
!     
      DO I=1,NP
        ITR(I) = 2
        IF(R0(I,1).GT.RIG(1,1).AND.R0(I,1).LT.RIG(2,1)
     &    .AND.R0(I,2).GT.RIG(1,2).AND.R0(I,2).LT.RIG(2,2)
     &    .AND.R0(I,3).GT.RIG(1,3).AND.R0(I,3).LT.RIG(2,3)
     &       ) ITR(I) = 1
        IF(R0(I,1).GT.THM(1,1).AND.R0(I,1).LT.THM(2,1)
     &    .AND.R0(I,2).GT.THM(1,2).AND.R0(I,2).LT.THM(2,2)
     &    .AND.R0(I,3).GT.THM(1,3).AND.R0(I,3).LT.THM(2,3)
     &       ) ITR(I) = 0
      ENDDO
!
!     Set just bottom bit to rigid
      RIG(1,1) = RMN(1,3) + .5d0
      RIG(2,1) = RMN(1,3) + 4.8d0
      DO I =1,NP
        ITR(I) = 0
        IF(R0(I,3).LT.RIG(2,1) ) ITR(I) = 1
        IF(R0(I,3).LT.RIG(1,1) ) ITR(I) = 2
      ENDDO
!
      RETURN
      END SUBROUTINE setthrmcust
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
