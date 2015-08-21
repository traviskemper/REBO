c
c calculated load on simulated graphite layer
c
c*********************************************************************
      subroutine load

      USE MPIvars
      USE dyn_array
      USE SPECIF
      USE MDSTP
      
      IMPLICIT none

      INCLUDE 'mpif.h'


C     local variable
      INTEGER :: i,jj,ierr,navg
      REAL*8 floadtot(3)

      kb = kb + 1

C navg should be equal to or less than MAXKB.
      navg = 100

      IF (kb.GT.(MAXKB-navg)) THEN
C      WRITE(*,*) 'load.f',kb
      do 235 i = 1, NP_node     ! General use
           if (itr_node(i).eq.3) then
             DO jj=1,3
                fload_node(jj) = fload_node(jj) + rnp(jj,NA(i))
             ENDDO
           endif
235      continue

      ENDIF

      IF(kb.ge.MAXKB) THEN
        CALL MPI_REDUCE(fload_node,floadtot,3,MPI_REAL8,MPI_SUM,0,
     &       MPI_COMM_WORLD,ierr)
        IF (mynode.EQ.0) THEN
           write(51,100)LSTEP,1.602189*floadtot(1)/navg,
     &     1.602189*floadtot(2)/navg,1.602189*floadtot(3)/navg
        ENDIF 
        kb = 0
        fload_node = 0.0d0
      ENDIF
      return
100   format(I7,3f15.4)
      end
c

