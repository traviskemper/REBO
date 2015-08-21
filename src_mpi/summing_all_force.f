       SUBROUTINE  sum_all_force
       
C sum of passed buffer force from neighbor nodes in each node
      USE dyn_array
      USE MPIvars

      IMPLICIT none

      INCLUDE 'mpif.h'


C local variable
      INTEGER    K,I,NN,MM

C      REAL*8 RNP_pi_buf(3,NP)

      RNP_buf=0.0d0
C      RNP_pi=0.0d0
C The center region force
C      DO I=1,NP_node
C        DO MM=1,3
C          RNP_pi(MM,NA(I))=RNP_node(MM,I)
C        ENDDO
C      ENDDO
C Unpack buffer region force
C      WRITE(*,*) mynode,'in sum_pi_force- N_force/4=',N_force/4
      IF (NNDP.EQ.1) THEN
        NN=0
        DO I=1,N_force/4
          NN=NN+1
          K=INT(buf_force_recv(NN))
          RNP_buf(1,K)=buf_force_recv(NN+1)
          RNP_buf(2,K)=buf_force_recv(NN+2)
          RNP_buf(3,K)=buf_force_recv(NN+3)
          NN=NN+3
        ENDDO
      ELSE
        NN=0
        DO I=1,N_force_RL/4
          NN=NN+1
          K=INT(buf_force_recv(NN))
          RNP_buf(1,K)=buf_force_recv(NN+1)
          RNP_buf(2,K)=buf_force_recv(NN+2)
          RNP_buf(3,K)=buf_force_recv(NN+3)
          NN=NN+3
        ENDDO
        DO I=1,N_force_corner/4
          NN=NN+1
          K=INT(buf_force_recv(NN))
          RNP_buf(1,K)=RNP_buf(1,K)+buf_force_recv(NN+1)
          RNP_buf(2,K)=RNP_buf(2,K)+buf_force_recv(NN+2)
          RNP_buf(3,K)=RNP_buf(3,K)+buf_force_recv(NN+3)
          NN=NN+3
        ENDDO
        DO I=1,N_force_UD/4
          NN=NN+1
          K=INT(buf_force_recv(NN))
          RNP_buf(1,K)=RNP_buf(1,K)+buf_force_recv(NN+1)
          RNP_buf(2,K)=RNP_buf(2,K)+buf_force_recv(NN+2)
          RNP_buf(3,K)=RNP_buf(3,K)+buf_force_recv(NN+3)
          NN=NN+3
        ENDDO
      ENDIF
C sum the pi force in center and from buffer

C      DO I=1,NP_node
C        K=NA(I)
C        DO MM=1,3
C          RNP_pi(MM,K)=RNP_node(MM,I)+RNP_buf(MM,K)
C        ENDDO
C      ENDDO
C       DO i=1,NP_node
C        IF (mynode.EQ.0) THEN
C          WRITE(31,999) NA(i),RNP_pi(1,NA(i)),RNP_pi(2,NA(i)),
C     &    RNP_pi(3,NA(i))
C        ENDIF
C        IF (mynode.EQ.1) THEN
C          WRITE(32,999) NA(i),RNP_pi(1,NA(i)),RNP_pi(2,NA(i)),
C     &    RNP_pi(3,NA(i))
C        ENDIF
C        IF (mynode.EQ.2) THEN
C          WRITE(33,999) NA(i),RNP_pi(1,NA(i)),RNP_pi(2,NA(i)),
C     &    RNP_pi(3,NA(i))
C        ENDIF
C        IF (mynode.EQ.3) THEN
C          WRITE(34,999) NA(i),RNP_pi(1,NA(i)),RNP_pi(2,NA(i)),
C     &    RNP_pi(3,NA(i))
C        ENDIF
C      ENDDO
C sum the force in caguts.f and pibond.f and ljguts.f
      DO I=1,NP_node
        K=NA(I)
        DO MM=1,3
          RNP(MM,K)=RNP(MM,K)+RNP_node(MM,I)+RNP_lj(MM,I)+RNP_buf(MM,K)
        ENDDO
      ENDDO
C       DO i=1,NP_node
C        IF (mynode.EQ.0) THEN
C          WRITE(31,999) NA(i),RNP(1,NA(i)),RNP(2,NA(i)),
C     &    RNP(3,NA(i))
C        ENDIF
C        IF (mynode.EQ.1) THEN
C          WRITE(32,999) NA(i),RNP(1,NA(i)),RNP(2,NA(i)),
C     &    RNP(3,NA(i))
C        ENDIF
C        IF (mynode.EQ.2) THEN
C          WRITE(33,999) NA(i),RNP(1,NA(i)),RNP(2,NA(i)),
C     &    RNP(3,NA(i))
C        ENDIF
C        IF (mynode.EQ.3) THEN
C          WRITE(34,999) NA(i),RNP(1,NA(i)),RNP(2,NA(i)),
C     &    RNP(3,NA(i))
C        ENDIF
C      ENDDO
C999   FORMAT(I5,3(E20.11))
      RETURN
      END

