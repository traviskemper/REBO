      SUBROUTINE XMOL
      
C This subroutine pass R0_node in every node to master node
      USE MPIvars
      USE dyn_array
      USE STRUCTR
      USE POTS

      IMPLICIT none

      INCLUDE 'mpif.h'


C local variable
      INTEGER  :: i,ii,NN,KK,N,K,ierr
     & ,NP_xmol,n_source
    
      INTEGER istat(MPI_STATUS_SIZE)
      REAL*8  r0_xmol(3)  

      NN=0
      DO i=1,NP_node
        NN=NN+1
        buf_xmol(NN)=float(NA(i))
        buf_xmol(NN+1)=float(ktype_node(i))
        buf_xmol(NN+2)=r0_node(1,i)
        buf_xmol(NN+3)=r0_node(2,i)
        buf_xmol(NN+4)=r0_node(3,i)
        NN=NN+4
      ENDDO

      IF (mynode.EQ.0) THEN
        WRITE(1,100) NP
        WRITE(1,1800) CTITLE
        DO 11 I=1,NP_node
           WRITE(1,350) NA(i),KT2(KTYPE_node(I)),
     &     (R0_node(N,I),N=1,3)
11      CONTINUE
        DO n_source=1,nprocs-1
          CALL MPI_RECV(NP_xmol,1,MPI_INTEGER,
     &    n_source,1,MPI_COMM_WORLD,istat,ierr)
          CALL MPI_RECV(buf_xmol_recv,5*NP_xmol,MPI_REAL8,
     &    n_source,2,MPI_COMM_WORLD,istat,ierr)
          NN=0
          DO ii=1,NP_xmol
            NN=NN+1
            K=INT(buf_xmol_recv(NN))
            KK=INT(buf_xmol_recv(NN+1))
            r0_xmol(1)=buf_xmol_recv(NN+2)
            r0_xmol(2)=buf_xmol_recv(NN+3)
            r0_xmol(3)=buf_xmol_recv(NN+4)
            NN=NN+4
            WRITE(1,350) K,KT2(KK),(r0_xmol(N),N=1,3)
          ENDDO
        ENDDO
      ELSE
        CALL MPI_SEND(NP_node,1,MPI_INTEGER,0
     &  ,1,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(buf_xmol,5*NP_node,MPI_REAL8,0
     &  ,2,MPI_COMM_WORLD,ierr)
      ENDIF

  100 FORMAT(I7)
  300 FORMAT(3E20.11)
  350 FORMAT(I7,I5,3E20.11)
  360 FORMAT(I7,3E20.11)
 1800 FORMAT(20A2)
      return
      end
