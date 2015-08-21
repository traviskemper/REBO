       SUBROUTINE pass_force
       
C This subroutine pass the pi force
      USE MPIvars
      USE dyn_array

      IMPLICIT none

      INCLUDE 'mpif.h'

C local variable
      INTEGER :: i,ierr

      INTEGER istat(MPI_STATUS_SIZE,4),ireq(4)
      INTEGER iistat(MPI_STATUS_SIZE,12),iireq(12)
      INTEGER node_x,node_y,node,NN_force(6)

      N_force=0
      N_force_corner=0
      N_force_RL=0
      N_force_UD=0
CCC 1-D parallel
      IF (NNDP.EQ.1) THEN
C all nodes pass to it's left and right node
        DO i=1,2
          IF (MOD(mynode+i,2).EQ.0) THEN
            CALL MPI_ISEND(N_force_passLL(3),1,MPI_INTEGER,
     &      nabor_node(1),1,MPI_COMM_WORLD,ireq(1),ierr)
            CALL MPI_ISEND(N_force_passRR(3),1,MPI_INTEGER,
     &      nabor_node(2),2,MPI_COMM_WORLD,ireq(2),ierr)
          ELSE
            CALL MPI_IRECV(N_force_recvRR(3),1,MPI_INTEGER,
     &      nabor_node(2),1,MPI_COMM_WORLD,ireq(3),ierr)
            CALL MPI_IRECV(N_force_recvLL(3),1,MPI_INTEGER,
     &      nabor_node(1),2,MPI_COMM_WORLD,ireq(4),ierr)
          ENDIF
        ENDDO
        CALL MPI_WAITALL (4,ireq,istat,ierr)
        N_force=4*N_force_recvRR(3)
        DO i=1,2
          IF (MOD(mynode+i,2).EQ.0) THEN
            CALL MPI_ISEND(buf_force_LL,4*N_force_passLL(3),
     &      MPI_REAL8,nabor_node(1),1,MPI_COMM_WORLD,ireq(1),ierr)
            CALL MPI_ISEND(buf_force_RR,
     &      4*N_force_passRR(3),MPI_REAL8,nabor_node(2),2,
     *      MPI_COMM_WORLD,ireq(2),ierr)
          ELSE
            CALL MPI_IRECV(buf_force_recv,4*N_force_recvRR(3),
     &      MPI_REAL8,nabor_node(2),1,MPI_COMM_WORLD,ireq(3),ierr)
            CALL MPI_IRECV(buf_force_recv(N_force+1),
     &      4*N_force_recvLL(3),MPI_REAL8,nabor_node(1),2,
     &      MPI_COMM_WORLD,ireq(4),ierr)
          ENDIF
        ENDDO
        CALL MPI_WAITALL (4,ireq,istat,ierr)
        N_force=N_force+4*N_force_recvLL(3)
C        WRITE(*,*) 'pass_pi_force- NPtot4_node=',NPtot4_node,NPmax_node
CCC 2-D parallel
      ELSE
C all nodes pass to its  left node
        DO i=1,2
          IF (MOD(mynode+i,2).EQ.0) THEN
            CALL MPI_ISEND(N_force_passLL(3),1,MPI_INTEGER,nabor_node(1)
     &      ,1,MPI_COMM_WORLD,iireq(1),ierr)
            CALL MPI_ISEND(N_force_passRR(3),1,MPI_INTEGER,nabor_node(2)
     &      ,2,MPI_COMM_WORLD,iireq(2),ierr)
            CALL MPI_ISEND(N_force_passLD(3),1,MPI_INTEGER,nabor_node(3)
     &      ,3,MPI_COMM_WORLD,iireq(3),ierr)
            CALL MPI_ISEND(N_force_passRU(3),1,MPI_INTEGER,nabor_node(4)
     &      ,4,MPI_COMM_WORLD,iireq(4),ierr)
            CALL MPI_ISEND(N_force_passRD(3),1,MPI_INTEGER,nabor_node(5)
     &      ,5,MPI_COMM_WORLD,iireq(5),ierr)
            CALL MPI_ISEND(N_force_passLU(3),1,MPI_INTEGER,nabor_node(6)
     &      ,6,MPI_COMM_WORLD,iireq(6),ierr)
          ELSE
            CALL MPI_IRECV(N_force_recvRR(3),1,MPI_INTEGER,nabor_node(2)
     &      ,1,MPI_COMM_WORLD,iireq(7),ierr)
            CALL MPI_IRECV(N_force_recvLL(3),1,MPI_INTEGER,nabor_node(1)
     &      ,2,MPI_COMM_WORLD,iireq(8),ierr)
            CALL MPI_IRECV(N_force_recvRU(3),1,MPI_INTEGER,nabor_node(4)
     &      ,3,MPI_COMM_WORLD,iireq(9),ierr)
            CALL MPI_IRECV(N_force_recvLD(3),1,MPI_INTEGER,nabor_node(3)
     &      ,4,MPI_COMM_WORLD,iireq(10),ierr)
            CALL MPI_IRECV(N_force_recvLU(3),1,MPI_INTEGER,nabor_node(6)
     &      ,5,MPI_COMM_WORLD,iireq(11),ierr)
            CALL MPI_IRECV(N_force_recvRD(3),1,MPI_INTEGER,nabor_node(5)
     &      ,6,MPI_COMM_WORLD,iireq(12),ierr)
          ENDIF
        ENDDO
        CALL MPI_WAITALL(12,iireq,iistat,ierr)
        NN_force(1)=4*N_force_recvRR(3)
        NN_force(2)=NN_force(1)+4*N_force_recvLL(3)
        NN_force(3)=NN_force(2)+4*N_force_recvRU(3)
        NN_force(4)=NN_force(3)+4*N_force_recvLD(3)
        NN_force(5)=NN_force(4)+4*N_force_recvLU(3)
        NN_force(6)=NN_force(5)+4*N_force_recvRD(3)
        N_force_RL=4*N_force_recvRR(3)+4*N_force_recvLL(3)
        N_force_corner=4*N_force_recvRU(3)+4*N_force_recvLD(3)+
     &                 4*N_force_recvLU(3)+4*N_force_recvRD(3)
        DO i=1,2
          IF (MOD(mynode+i,2).EQ.0) THEN
            CALL MPI_ISEND(buf_force_LL,4*N_force_passLL(3),MPI_REAL8,
     &      nabor_node(1),1,MPI_COMM_WORLD,iireq(1),ierr)
            CALL MPI_ISEND(buf_force_RR,4*N_force_passRR(3),MPI_REAL8,
     &      nabor_node(2),2,MPI_COMM_WORLD,iireq(2),ierr)
            CALL MPI_ISEND(buf_force_LD,4*N_force_passLD(3),MPI_REAL8,
     &      nabor_node(3),3,MPI_COMM_WORLD,iireq(3),ierr)
            CALL MPI_ISEND(buf_force_RU,4*N_force_passRU(3),MPI_REAL8,
     &      nabor_node(4),4,MPI_COMM_WORLD,iireq(4),ierr)
            CALL MPI_ISEND(buf_force_RD,4*N_force_passRD(3),MPI_REAL8,
     &      nabor_node(5),5,MPI_COMM_WORLD,iireq(5),ierr)
            CALL MPI_ISEND(buf_force_LU,4*N_force_passLU(3),MPI_REAL8,
     &      nabor_node(6),6,MPI_COMM_WORLD,iireq(6),ierr)
          ELSE
            CALL MPI_IRECV(buf_force_recv,
     &      4*N_force_recvRR(3),MPI_REAL8,nabor_node(2),
     &      1,MPI_COMM_WORLD,iireq(7),ierr)
            CALL MPI_IRECV(buf_force_recv(NN_force(1)+1),
     &      4*N_force_recvLL(3),MPI_REAL8,nabor_node(1),
     &      2,MPI_COMM_WORLD,iireq(8),ierr)
            CALL MPI_IRECV(buf_force_recv(NN_force(2)+1),
     &      4*N_force_recvRU(3),MPI_REAL8,nabor_node(4),
     &      3,MPI_COMM_WORLD,iireq(9),ierr)
            CALL MPI_IRECV(buf_force_recv(NN_force(3)+1),
     &      4*N_force_recvLD(3),MPI_REAL8,nabor_node(3),
     &      4,MPI_COMM_WORLD,iireq(10),ierr)
            CALL MPI_IRECV(buf_force_recv(NN_force(4)+1),
     &      4*N_force_recvLU(3),MPI_REAL8,nabor_node(6),
     &      5,MPI_COMM_WORLD,iireq(11),ierr)
            CALL MPI_IRECV(buf_force_recv(NN_force(5)+1),
     &      4*N_force_recvRD(3),MPI_REAL8,nabor_node(5),
     &      6,MPI_COMM_WORLD,iireq(12),ierr)
          ENDIF
        ENDDO
        CALL MPI_WAITALL(12,iireq,iistat,ierr)
C transfer node number to node_x and node_y
        node_x=MOD((mynode+1),NDP(NPD_i))
        IF (node_x.EQ.0) THEN
          node_x=NDP(NPD_i)
        ENDIF
        node_y=INT((mynode+1)/(NDP(NPD_i)+1.0d-12))+1
C all nodes pass to its DOWN node
        DO i=1,2
          IF (MOD(node_y+i,2).EQ.0) THEN
            CALL MPI_ISEND(N_force_passDD(3),1,MPI_INTEGER,nabor_node(7)
     &      ,1,MPI_COMM_WORLD,ireq(1),ierr)
            CALL MPI_ISEND(N_force_passUU(3),1,MPI_INTEGER,nabor_node(8)
     &      ,2,MPI_COMM_WORLD,ireq(2),ierr)
          ELSE
            CALL MPI_IRECV(N_force_recvUU(3),1,MPI_INTEGER,nabor_node(8)
     &      ,1,MPI_COMM_WORLD,ireq(3),ierr)
            CALL MPI_IRECV(N_force_recvDD(3),1,MPI_INTEGER,nabor_node(7)
     &      ,2,MPI_COMM_WORLD,ireq(4),ierr)
          ENDIF
        ENDDO
        CALL MPI_WAITALL (4,ireq,istat,ierr)
        N_force=NN_force(6)+4*N_force_recvUU(3)
        N_force_UD=4*N_force_recvUU(3)+4*N_force_recvDD(3)
        DO i=1,2
          IF (MOD(node_y+i,2).EQ.0) THEN
            CALL MPI_ISEND(buf_force_DD,4*N_force_passDD(3),MPI_REAL8,
     &      nabor_node(7),1,MPI_COMM_WORLD,ireq(1),ierr)
            CALL MPI_ISEND(buf_force_UU,4*N_force_passUU(3),MPI_REAL8,
     &      nabor_node(8),2,MPI_COMM_WORLD,ireq(2),ierr)
          ELSE
            CALL MPI_IRECV(buf_force_recv(NN_force(6)+1),
     &      4*N_force_recvUU(3),MPI_REAL8,nabor_node(8),
     &      1,MPI_COMM_WORLD,ireq(3),ierr)
            CALL MPI_IRECV(buf_force_recv(N_force+1),
     &      4*N_force_recvDD(3),MPI_REAL8,nabor_node(7),
     &      2,MPI_COMM_WORLD,ireq(4),ierr)
          ENDIF
        ENDDO
        CALL MPI_WAITALL (4,ireq,istat,ierr)
        N_force=N_force+4*N_force_recvDD(3)
C        WRITE(*,*) 'pass_force UU finish'
      ENDIF
C      WRITE(*,*) mynode,'pass_pi_force - finish',N_force/4
C      WRITE(*,*) mynode,N_force_RL/4,N_force_corner/4,N_force_UD/4
      RETURN
      END
