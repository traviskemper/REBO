c
c move some atoms
c
c*********************************************************************
      subroutine move_tip

      USE dyn_array
      USE MPIvars
      USE STRUCTR
      USE SPECIF

      IMPLICIT none

      INTEGER :: i,j
C
C move tip down

      do i = 1,NPtot3_node
           if(itr_node(i).eq.3) then
                do j = 1, 3
                   r0_node(j,i) = r0_node(j,i) + disp(j)
                   r0_node(j,i) = r0_node(j,i) - 
     &             CUBE(j)*ANINT(R0_node(j,i)/CUBE(j))
                enddo
           endif
      enddo
c
c add velocity to carbon tip for first step
c      if(ttime.le.delta) then
c           do i=1,np
c                if((ktype(i).le.2).and.(itr(i).ne.2)) then
c                     r1(i,2) = -0.0005d0
c                endif
c           enddo
c      endif
      return
      end SUBROUTINE move_tip


c*********************************************************************
