C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Programer:
C     Travis Kemper
C     Department of Materials SCience and ENgineering
C     University of FLorida
C     traviskemper@ufl.edu
C
C     Version 1.0 9/25/09 T. W. Kemper
! shift the coordinates a small distance for parallel code
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      SUBROUTINE shift 
!
      USE STRUCTR
      USE MPIvars
!
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
      REAL*8  :: SHF(3)
      INTEGER :: J,N
!
      IF (mynode.EQ.0) THEN
      SHF(1) = 1.0d0
      SHF(2) = 1.0d0
      SHF(3) = 1.0d0
!
      RMN(1,:) = (/1.d16,1.0d16,1.0d16/)
      RMN(2,:) = (/-1.0d16,-1.0d16,-1.0d16/)
      DO J=1,NP
        DO N=1,3
          R0i(N,J) = R0i(N,J) + SHF(N) 
          R0i(N,J)=R0i(N,J)-CUBE(N)*ANINT(R0i(N,J)/CUBE(N))
!         Calculate r_min and r_max
          IF ( R0i(N,J) .GT. RMN(2,N) ) RMN(2,N)=R0i(N,J)
          IF ( R0i(N,J) .LT. RMN(1,N) ) RMN(1,N)=R0i(N,J)
        ENDDO
      ENDDO
!
      ENDIF
      END SUBROUTINE shift 
