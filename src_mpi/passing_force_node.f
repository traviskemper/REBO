       SUBROUTINE pass_force_node

C This subroutine find and pass pi force between nodes
      USE dyn_array
      USE MPIvars

      IMPLICIT none

      INCLUDE 'mpif.h'

C local variable
      INTEGER     K,I,NN

      CALL find_pi_force
      CALL pass_force
      CALL sum_all_force
      CALL find_pass_force
      CALL pass_force
C      WRITE(*,*) mynode,'pass_force_node',N_force/4
      NN=0
      DO I=1,N_force/4
        NN=NN+1
        K=INT(buf_force_recv(NN))
        RNP(1,K)=buf_force_recv(NN+1)
        RNP(2,K)=buf_force_recv(NN+2)
        RNP(3,K)=buf_force_recv(NN+3)
C        IF (mynode.EQ.0) THEN
C          WRITE(31,999) K,buf_force_recv(NN+1),buf_force_recv(NN+2),
C     &    buf_force_recv(NN+3)
C        ENDIF
C        IF (mynode.EQ.1) THEN
C          WRITE(32,999) K,buf_force_recv(NN+1),buf_force_recv(NN+2),
C     &    buf_force_recv(NN+3)
C        ENDIF
c        IF (mynode.EQ.2) THEN
c          WRITE(33,999) K,buf_force_recv(NN+1),buf_force_recv(NN+2),
c     &    buf_force_recv(NN+3)
c        ENDIF
c        IF (mynode.EQ.3) THEN
c          WRITE(34,999) K,buf_force_recv(NN+1),buf_force_recv(NN+2),
c     &    buf_force_recv(NN+3)
c        ENDIF
        NN=NN+3
      ENDDO
C      DO i=1,NP_node
C        IF (mynode.EQ.0) THEN
C          WRITE(33,999) NA(i),RNP(1,NA(i)),RNP(2,NA(i)),
C     &    RNP(3,NA(i))
C        ENDIF
C        IF (mynode.EQ.1) THEN
C          WRITE(34,999) NA(i),RNP(1,NA(i)),RNP(2,NA(i)),
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
999   FORMAT(I5,3(E20.11))

C      WRITE(*,*) 'pass_force_node - finish'
      RETURN
      END
