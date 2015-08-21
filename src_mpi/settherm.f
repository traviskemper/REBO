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
      USE MPIvars
!
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
      INTEGER :: M,I,D,N
      IF ( mynode.EQ.0) THEN
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
        ITRi(I) = 2
        IF(R0i(1,I).GT.RIG(1,1).AND.R0i(1,I).LT.RIG(2,1)
     &    .AND.R0i(2,I).GT.RIG(1,2).AND.R0i(2,I).LT.RIG(2,2)
     &    .AND.R0i(3,I).GT.RIG(1,3).AND.R0i(3,I).LT.RIG(2,3)
     &       ) ITRi(I) = 1
        IF(R0i(1,I).GT.THM(1,1).AND.R0i(1,I).LT.THM(2,1)
     &    .AND.R0i(2,I).GT.THM(1,2).AND.R0i(2,I).LT.THM(2,2)
     &    .AND.R0i(3,I).GT.THM(1,3).AND.R0i(3,I).LT.THM(2,3)
     &       ) ITRi(I) = 0
      ENDDO
!
      ENDIF
      RETURN
      END SUBROUTINE setthrmbox
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  Set thermostates
!
      SUBROUTINE setthrmcust
!
      USE STRUCTR
      USE SPECIF
      USE MPIvars
!
!
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
      INTEGER :: M,I,D,N
!     Find rigid  amd thermostat box sizes 
      IF ( mynode.EQ.0) THEN
c$$$      DO D=1,3
c$$$        M=1 !min
c$$$        RIG(M,D) = RMN(M,D) + RBOX(M,D)
c$$$        IF( RBOX(M,D).LT.0.1d0) RIG(M,D)=RMN(M,D)-1.0d0 !Account for round off
c$$$        THM(M,D) = RIG(M,D) + TBOX(M,D)
c$$$        IF( TBOX(M,D).LT.0.1d0) THM(M,D)=RMN(M,D)-1.0d0 !Account for round off
c$$$        M=2 !max
c$$$        RIG(M,D) = RMN(M,D) - RBOX(M,D)
c$$$        IF( RBOX(M,D).LT.0.1d0) RIG(M,D)=RMN(M,D)+1.0d0 !Account for round off
c$$$        THM(M,D) = RIG(M,D) - TBOX(M,D)
c$$$        IF( TBOX(M,D).LT.0.1d0) THM(M,D)=RMN(M,D)+1.0d0 !Account for round off
c$$$      ENDDO
c$$$!     
c$$$      DO I=1,NP
c$$$        ITRi(I) = 2
c$$$        IF(R0i(1,I).GT.RIG(1,1).AND.R0i(1,I).LT.RIG(2,1)
c$$$     &    .AND.R0i(2,I).GT.RIG(1,2).AND.R0i(2,I).LT.RIG(2,2)
c$$$     &    .AND.R0i(3,I).GT.RIG(1,3).AND.R0i(3,I).LT.RIG(2,3)
c$$$     &       ) ITRi(I) = 1
c$$$        IF(R0i(1,I).GT.THM(1,1).AND.R0i(1,I).LT.THM(2,1)
c$$$     &    .AND.R0i(2,I).GT.THM(1,2).AND.R0i(2,I).LT.THM(2,2)
c$$$     &    .AND.R0i(3,I).GT.THM(1,3).AND.R0i(3,I).LT.THM(2,3)
c$$$     &       ) ITRi(I) = 0
c$$$      ENDDO
!
!     Set just bottom bit to rigid
      RIG(1,1) = RMN(1,3) + .5d0
      RIG(2,1) = RMN(1,3) + 4.8d0
      DO I =1,NP
        ITRi(I) = NTHRM
        IF(R0i(3,I).LT.RIG(2,1) ) ITRi(I) = 1
        IF(R0i(3,I).LT.RIG(1,1) ) ITRi(I) = 2
      ENDDO
!
      ENDIF
      RETURN
      END SUBROUTINE setthrmcust
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  Set thermostates
!
      SUBROUTINE prestherm
!
      USE STRUCTR
      USE SPECIF
      USE MPIvars
!
!
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
      INTEGER :: M,I,D,N
!     Find rigid  amd thermostat box sizes 
      IF ( mynode.EQ.0) THEN
!
!     Set just bottom bit to rigid
      THM(1,1) = RMN(1,SDIR) + PRBUF + 10.0d0
      THM(2,1) = RMN(2,SDIR) - PRBUF - 10.0d0
      RIG(1,1) = RMN(1,SDIR) + PRBUF
      RIG(2,1) = RMN(2,SDIR) - PRBUF
      DO I =1,NP
        ITRi(I) = 0 !NTHRM
        IF(R0i(3,I).LT.THM(1,1) ) ITRi(I) = 1
        IF(R0i(3,I).GT.THM(2,1) ) ITRi(I) = 1
        IF(R0i(3,I).LT.RIG(1,1) ) ITRi(I) = 2
        IF(R0i(3,I).GT.RIG(2,1) ) ITRi(I) = 3
      ENDDO
!
      ENDIF
      RETURN
      END SUBROUTINE prestherm
