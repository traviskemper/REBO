      SUBROUTINE MODEL

      USE dyn_array      
      USE MDSTP

      IMPLICIT none
c
c set forces to zero
C

      RNP = 0.0d0
      RNP_node=0.0d0
      RNP_lj=0.0d0


C set potential energy to zero
      tote = 0.0d0

C       sss = noa(1)+noa(2)+noa(3)+noa(4)+noa(5)
c       write(274,*) ipot,sss
C       if((ipot.eq.1).and.(sss.ne.0)) CALL CAGUTS
C       if((ipot.eq.2).and.((noa(1)+noa(2)+noa(3)).ne.0)) call tight_bind
C       if(ILJ.eq.1) call ljguts
c      call ljcont
      return
      end
