      SUBROUTINE pass_node
C
C After predictor or corrector the atom's posiiton is changed
C When the atoms move accross the node region, that atom should be
C re-assign to the corresct node (only atoms in the neighbor region will move accross the node)
C And also when the atom move across the center region and neighbor region the atom also need
C to be re-assign to the correct region
C This work should do twise for pibond and Lennard Jones, because atoms are calculated seperatly
C This subroutine do this work for L-J part
C
C Keep NP_node correct and transfer neighbor region atoms from neighbor node
C First Sort atoms in the first inner shell of center region
C Second check if the atom go out of the node
C Third make the list for atoms that go out of the  node
C Forth pass those atoms to the node they should be
C Fifth pass the neighbor region atoms from neighbor nodes
C
C The program should devide to three catalogy: nprocs_y=1, nprocs_x=1, and both nprocs_y nprocs_x GT 1
C

      USE MPIvars
      USE dyn_array

      IMPLICIT none

      INCLUDE 'mpif.h'

C global variable

C local variable
      INTEGER istat(MPI_STATUS_SIZE,4),ireq(4)
      INTEGER iistat(MPI_STATUS_SIZE,12),iireq(12)
      INTEGER node_x,node_y,node,NN_data(6)
C N_data changed to global not sure if this correct
C -travisk

      INTEGER :: N_data,i,ierr

CCC 1-D parallel
C pass atoms go out of the center region to neighbor nodes
      IF (NNDP.EQ.1) THEN
C all nodes pass to it's left node
        DO i=1,2
          IF (MOD(mynode+i,2).EQ.0) THEN
            CALL MPI_ISEND(N_atom_passLL,3,MPI_INTEGER,nabor_node(1)
     &      ,1,MPI_COMM_WORLD,ireq(1),ierr)
            CALL MPI_ISEND(N_atom_passRR,3,MPI_INTEGER,nabor_node(2)
     &      ,2,MPI_COMM_WORLD,ireq(2),ierr)
          ELSE
            CALL MPI_IRECV(N_atom_recvRR,3,MPI_INTEGER,nabor_node(2)
     &      ,1,MPI_COMM_WORLD,ireq(3),ierr)
            CALL MPI_IRECV(N_atom_recvLL,3,MPI_INTEGER,nabor_node(1)
     &      ,2,MPI_COMM_WORLD,ireq(4),ierr)
          ENDIF
        ENDDO
        CALL MPI_WAITALL(4,ireq,istat,ierr)
C        WRITE(*,*) mynode,'N_atom_recvRR,LL',N_atom_recvRR(3),
C     &             N_atom_recvLL(3)
        N_data=15*N_atom_recvRR(3)
        DO i=1,2
          IF (MOD(mynode+i,2).EQ.0) THEN
            CALL MPI_ISEND(buf_atom_passLL,15*N_atom_passLL(3),
     &      MPI_REAL8,nabor_node(1),1,MPI_COMM_WORLD,ireq(1),ierr)
            CALL MPI_ISEND(buf_atom_passRR,15*N_atom_passRR(3),
     &      MPI_REAL8,nabor_node(2),2,MPI_COMM_WORLD,ireq(2),ierr)
          ELSE
            CALL MPI_IRECV(buf_atom_recv,15*N_atom_recvRR(3),
     &      MPI_REAL8,nabor_node(2),1,MPI_COMM_WORLD,ireq(3),ierr)
            CALL MPI_IRECV(buf_atom_recv(N_data+1),15*N_atom_recvLL(3),
     &      MPI_REAL8,nabor_node(1),2,MPI_COMM_WORLD,ireq(4),ierr)
          ENDIF
        ENDDO
        CALL MPI_WAITALL(4,ireq,istat,ierr)
CCC 2-D parallel
C pass atoms go out of the center region to neighbor nodes
      ELSE
C all nodes pass to its  left node
        DO i=1,2
          IF (MOD(mynode+i,2).EQ.0) THEN
            CALL MPI_ISEND(N_atom_passLL,3,MPI_INTEGER,nabor_node(1)
     &      ,1,MPI_COMM_WORLD,iireq(1),ierr)
            CALL MPI_ISEND(N_atom_passRR,3,MPI_INTEGER,nabor_node(2)
     &      ,2,MPI_COMM_WORLD,iireq(2),ierr)
            CALL MPI_ISEND(N_atom_passLD,3,MPI_INTEGER,nabor_node(3)
     &      ,3,MPI_COMM_WORLD,iireq(3),ierr)
            CALL MPI_ISEND(N_atom_passRU,3,MPI_INTEGER,nabor_node(4)
     &      ,4,MPI_COMM_WORLD,iireq(4),ierr)
            CALL MPI_ISEND(N_atom_passRD,3,MPI_INTEGER,nabor_node(5)
     &      ,5,MPI_COMM_WORLD,iireq(5),ierr)
            CALL MPI_ISEND(N_atom_passLU,3,MPI_INTEGER,nabor_node(6)
     &      ,6,MPI_COMM_WORLD,iireq(6),ierr)
          ELSE
            CALL MPI_IRECV(N_atom_recvRR,3,MPI_INTEGER,nabor_node(2)
     &      ,1,MPI_COMM_WORLD,iireq(7),ierr)
            CALL MPI_IRECV(N_atom_recvLL,3,MPI_INTEGER,nabor_node(1)
     &      ,2,MPI_COMM_WORLD,iireq(8),ierr)
            CALL MPI_IRECV(N_atom_recvRU,3,MPI_INTEGER,nabor_node(4)
     &      ,3,MPI_COMM_WORLD,iireq(9),ierr)
            CALL MPI_IRECV(N_atom_recvLD,3,MPI_INTEGER,nabor_node(3)
     &      ,4,MPI_COMM_WORLD,iireq(10),ierr)
            CALL MPI_IRECV(N_atom_recvLU,3,MPI_INTEGER,nabor_node(6)
     &      ,5,MPI_COMM_WORLD,iireq(11),ierr)
            CALL MPI_IRECV(N_atom_recvRD,3,MPI_INTEGER,nabor_node(5)
     &      ,6,MPI_COMM_WORLD,iireq(12),ierr)
          ENDIF
        ENDDO
        CALL MPI_WAITALL(12,iireq,iistat,ierr)
        NN_data(1)=15*N_atom_recvRR(3)
        NN_data(2)=NN_data(1)+15*N_atom_recvLL(3)
        NN_data(3)=NN_data(2)+15*N_atom_recvRU(3)
        NN_data(4)=NN_data(3)+15*N_atom_recvLD(3)
        NN_data(5)=NN_data(4)+15*N_atom_recvLU(3)
        NN_data(6)=NN_data(5)+15*N_atom_recvRD(3)
        DO i=1,2
          IF (MOD(mynode+i,2).EQ.0) THEN
            CALL MPI_ISEND(buf_atom_passLL,15*N_atom_passLL(3),
     &      MPI_REAL8,nabor_node(1),1,MPI_COMM_WORLD,iireq(1),ierr)
            CALL MPI_ISEND(buf_atom_passRR,15*N_atom_passRR(3),
     &      MPI_REAL8,nabor_node(2),2,MPI_COMM_WORLD,iireq(2),ierr)
            CALL MPI_ISEND(buf_atom_passLD,15*N_atom_passLD(3),
     &      MPI_REAL8,nabor_node(3),3,MPI_COMM_WORLD,iireq(2),ierr)
            CALL MPI_ISEND(buf_atom_passRU,15*N_atom_passRU(3),
     &      MPI_REAL8,nabor_node(4),4,MPI_COMM_WORLD,iireq(4),ierr)
            CALL MPI_ISEND(buf_atom_passRD,15*N_atom_passRD(3),
     &      MPI_REAL8,nabor_node(5),5,MPI_COMM_WORLD,iireq(5),ierr)
            CALL MPI_ISEND(buf_atom_passLU,15*N_atom_passLU(3),
     &      MPI_REAL8,nabor_node(6),6,MPI_COMM_WORLD,iireq(6),ierr)
          ELSE
            CALL MPI_IRECV(buf_atom_recv,
     &      15*N_atom_recvRR(3),MPI_REAL8,nabor_node(2),1,
     &      MPI_COMM_WORLD,iireq(7),ierr)
            CALL MPI_IRECV(buf_atom_recv(NN_data(1)+1),
     &      15*N_atom_recvLL(3),MPI_REAL8,nabor_node(1),2,
     &      MPI_COMM_WORLD,iireq(8),ierr)
            CALL MPI_IRECV(buf_atom_recv(NN_data(2)+1),
     &      15*N_atom_recvRU(3),MPI_REAL8,nabor_node(4),3,
     &      MPI_COMM_WORLD,iireq(9),ierr)
            CALL MPI_IRECV(buf_atom_recv(NN_data(3)+1),
     &      15*N_atom_recvLD(3),MPI_REAL8,nabor_node(3),4,
     &      MPI_COMM_WORLD,iireq(10),ierr)
            CALL MPI_IRECV(buf_atom_recv(NN_data(4)+1),
     &      15*N_atom_recvLU(3),MPI_REAL8,nabor_node(6),5,
     &      MPI_COMM_WORLD,iireq(11),ierr)
            CALL MPI_IRECV(buf_atom_recv(NN_data(5)+1),
     &      15*N_atom_recvRD(3),MPI_REAL8,nabor_node(5),6,
     &      MPI_COMM_WORLD,iireq(12),ierr)
          ENDIF
        ENDDO
        CALL MPI_WAITALL(12,iireq,iistat,ierr)
C        WRITE(*,*) mynode,'pass_node LL finish'
C transfer node number to node_x and node_y
        node_x=MOD((mynode+1),NDP(NPD_i))
        IF (node_x.EQ.0) THEN
          node_x=NDP(NPD_i)
        ENDIF
        node_y=INT((mynode+1)/(NDP(NPD_i)+1.0d-12))+1
C all nodes pass to its DOWN node
        DO i=1,2
          IF (MOD(node_y+i,2).EQ.0) THEN
            CALL MPI_ISEND(N_atom_passDD,3,MPI_INTEGER,nabor_node(7)
     &      ,1,MPI_COMM_WORLD,ireq(1),ierr)
            CALL MPI_ISEND(N_atom_passUU,3,MPI_INTEGER,nabor_node(8)
     &      ,2,MPI_COMM_WORLD,ireq(2),ierr)
          ELSE
            CALL MPI_IRECV(N_atom_recvUU,3,MPI_INTEGER,nabor_node(8)
     &      ,1,MPI_COMM_WORLD,ireq(3),ierr)
            CALL MPI_IRECV(N_atom_recvDD,3,MPI_INTEGER,nabor_node(7)
     &      ,2,MPI_COMM_WORLD,ireq(4),ierr)
          ENDIF
        ENDDO
        CALL MPI_WAITALL(4,ireq,istat,ierr)
        N_data=NN_data(6)+15*N_atom_recvUU(3)
        DO i=1,2
          IF (MOD(node_y+i,2).EQ.0) THEN
            CALL MPI_ISEND(buf_atom_passDD,15*N_atom_passDD(3),
     &      MPI_REAL8,nabor_node(7),1,MPI_COMM_WORLD,ireq(1),ierr)
            CALL MPI_ISEND(buf_atom_passUU,15*N_atom_passUU(3),
     &      MPI_REAL8,nabor_node(8),2,MPI_COMM_WORLD,ireq(2),ierr)
          ELSE
            CALL MPI_IRECV(buf_atom_recv(NN_data(6)+1),
     &      15*N_atom_recvUU(3),MPI_REAL8,nabor_node(8),1,
     &      MPI_COMM_WORLD,ireq(3),ierr)
            CALL MPI_IRECV(buf_atom_recv(N_data+1),
     &      15*N_atom_recvDD(3),MPI_REAL8,nabor_node(7),2,
     &      MPI_COMM_WORLD,ireq(4),ierr)
          ENDIF
        ENDDO
        CALL MPI_WAITALL(4,ireq,istat,ierr)
        N_data=N_data+15*N_atom_recvDD(3)
      ENDIF
C      WRITE(*,*) mynode,'pass_node finish'
      RETURN
      END
