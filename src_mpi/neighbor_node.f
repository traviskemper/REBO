      SUBROUTINE nebor_node
C
C After predictor or corrector the atom's posiiton is changed
C When the atoms move accross the node region, that atom should be
C re-assign to the corresct node (only atoms in the neighbor region will move accross the node)
C And also when the atom move across the center region and neighbor region the atom also need
C to be re-assign to the correct region
C This work should do twise for pibond and Lennard Jones, because atoms are calculated seperatly
C This subroutine do this work for L-J part
C
C Keep NP_node_lj correct and transfer neighbor region atoms from neighbor node
C First Sort atoms in the first inner shell of center region
C Second check if the atom go out of the node
C Third make the list for atoms that go out of the  node
C Forth pass those atoms to the node they should be
C Fifth pass the neighbor region atoms from neighbor nodes
C
C The program should devide to three catalogy: nprocs_y=1, nprocs_x=1, and both nprocs_y nprocs_x GT 1
C

      USE MPIvars

      IMPLICIT none

      INCLUDE 'mpif.h'

C global variable

C local variable
      INTEGER L_node,R_node,U_node,D_node,RU_node,RD_node,LU_node,
     &LD_node,node_x,node_y
      INTEGER NL_node,NR_node,NU_node,ND_node,NRU_node,NRD_node,
     &NLU_node,NLD_node

CCC 1-D parallel
C Find the neighbor node
      IF (NNDP.EQ.1) THEN
        NL_node=mynode-1
        IF (NL_node.LT.0) THEN
          NL_node=NL_node+NDP(NPD_i)
        ENDIF
        nabor_node(1)=NL_node
        NR_node=mynode+1
        IF (NR_node.GT.nprocs-1) THEN
          NR_node=NR_node-NDP(NPD_i)
        ENDIF
        nabor_node(2)=NR_node
C        WRITE(*,*) mynode,'nabor_node',nabor_node(1),nabor_node(2)
CCC 2-D parallel
C Find the neighbor node
      ELSE
        node_x=MOD((mynode+1),NDP(NPD_i))
        IF (node_x.EQ.0) THEN
          node_x=NDP(NPD_i)
        ENDIF
        node_y=INT((mynode+1)/(NDP(NPD_i)+1.0d-12))+1
C        WRITE(*,*) mynode,node_x,node_y
        L_node=node_x-1
        IF (L_node.LE.0) THEN
          L_node=L_node+NDP(NPD_i)
        ENDIF
        NL_node=L_node+(node_y-1)*NDP(NPD_i)-1
        R_node=node_x+1
        IF (R_node.GT.NDP(NPD_i)) THEN
          R_node=R_node-NDP(NPD_i)
        ENDIF
        NR_node=R_node+(node_y-1)*NDP(NPD_i)-1
        D_node=node_y-1
        IF (D_node.LE.0) THEN
          D_node=D_node+NDP(NPD_j)
        ENDIF
        ND_node=node_x+(D_node-1)*NDP(NPD_i)-1
        U_node=node_y+1
        IF (U_node.GT.NDP(NPD_j)) THEN
          U_node=U_node-NDP(NPD_j)
        ENDIF
        NU_node=node_x+(U_node-1)*NDP(NPD_i)-1
        NLD_node=L_node+(D_node-1)*NDP(NPD_i)-1
        NLU_node=L_node+(U_node-1)*NDP(NPD_i)-1
        NRD_node=R_node+(D_node-1)*NDP(NPD_i)-1
        NRU_node=R_node+(U_node-1)*NDP(NPD_i)-1

        nabor_node(1)=NL_node
        nabor_node(2)=NR_node
        nabor_node(3)=NLD_node
        nabor_node(6)=NLU_node
        nabor_node(5)=NRD_node
        nabor_node(4)=NRU_node
        nabor_node(7)=ND_node
        nabor_node(8)=NU_node
C         WRITE(*,999) mynode,NL_node,NR_node,ND_node,NU_node,
C     &  NLD_node,NLU_node,NRD_node,NRU_node
      ENDIF
C      WRITE(*,*) 'finish nebor_node'
999   FORMAT(9(I5))
      RETURN
      END SUBROUTINE nebor_node
