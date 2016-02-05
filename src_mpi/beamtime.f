C  This program exports the time steps need to run a complete beam simulation
C   based on the distance and speed of the last atom in the coord.d file
C   syntax : beamtime < [Dep File coord].d
!
      SUBROUTINE beamtime
!
      USE beam
      USE STRUCTR
      USE SPECIF
      USE MPIvars
      USE dyn_array
!
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
      REAL :: V(2)
!
      IF (mynode.EQ.0) THEN
        V(1) = R0i(DEPD,NP) 
        V(2) = R1i(DEPD,NP) 
!
       STPS= INT( -V(1)/V(2)/DELTA )
!
      ENDIF
      RETURN
      END SUBROUTINE 
