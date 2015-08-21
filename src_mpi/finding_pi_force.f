       SUBROUTINE find_pi_force
       
C This subroutine do the passing of pibond.f force
C First find the atoms that need to be passed
C These atoms are the atoms within 3 buffer layers (from NP_node+1 to NPtot3_node)
C Second make the list to pass
C The list will be made by the passing direction
C Third pass the forces
C Forth sum the forces

      USE dyn_array
      USE MPIvars   

      IMPLICIT none

      INCLUDE 'mpif.h'

C Local variables
      INTEGER :: NRR3,NLL3,NDD3,NUU3,NRD3,NRU3,NLD3,NLU3
     & ,ii,NN

C For 1-D parallel
      IF (NNDP.EQ.1) THEN
        NRR=NP_node+N_atom_recvRR(1)
        NLL=NRR+N_atom_recvLL(1)
        NRR3=NLL+(N_atom_recvRR(2)-N_atom_recvRR(1))
        NLL3=NRR3+(N_atom_recvLL(2)-N_atom_recvLL(1))
C Right buffer force to right node
        NN=0
        DO ii=NP_node+1,NRR
          NN=NN+1
C          IF (NN.GT.4*NBMAX1) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX1, increase Density in
C     &      parameters.inc fd_pi_force'
C            CALL write_error
C          ENDIF
          buf_force_RR(NN)=float(NA(ii))
          buf_force_RR(NN+1)=RNP_node(1,ii)
          buf_force_RR(NN+2)=RNP_node(2,ii)
          buf_force_RR(NN+3)=RNP_node(3,ii)
          NN=NN+3
        ENDDO
        DO ii=NLL+1,NRR3
          NN=NN+1
C          IF (NN.GT.4*NBMAX1) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX1, increase Density in
C     &      parameters.inc fd_pi_force'
C            CALL write_error
C          ENDIF
          buf_force_RR(NN)=float(NA(ii))
          buf_force_RR(NN+1)=RNP_node(1,ii)
          buf_force_RR(NN+2)=RNP_node(2,ii)
          buf_force_RR(NN+3)=RNP_node(3,ii)
          NN=NN+3
        ENDDO
        N_force_passRR(3)=NN/4
C        WRITE(*,*) mynode,'#right buffer forces=',NN
C Left buffer force to left node
        NN=0
        DO ii=NRR+1,NLL
          NN=NN+1
C          IF (NN.GT.4*NBMAX1) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX1, increase Density in
C     &      parameters.inc fd_pi_force'
C            CALL write_error
C          ENDIF
          buf_force_LL(NN)=float(NA(ii))
          buf_force_LL(NN+1)=RNP_node(1,ii)
          buf_force_LL(NN+2)=RNP_node(2,ii)
          buf_force_LL(NN+3)=RNP_node(3,ii)
          NN=NN+3
        ENDDO
        DO ii=NRR3+1,NPtot3_node
          NN=NN+1
C          IF (NN.GT.4*NBMAX1) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX1, increase Density in
C     &      parameters.inc fd_pi_force'
C            CALL write_error
C          ENDIF
          buf_force_LL(NN)=float(NA(ii))
          buf_force_LL(NN+1)=RNP_node(1,ii)
          buf_force_LL(NN+2)=RNP_node(2,ii)
          buf_force_LL(NN+3)=RNP_node(3,ii)
          NN=NN+3
        ENDDO
        N_force_passLL(3)=NN/4
C        WRITE(*,*) mynode,'#left buffer forces=',NN
C For 2-D paralle
      ELSE
        NRR=NP_node+N_atom_recvRR(1)
        NLL=NRR+N_atom_recvLL(1)
        NRU=NLL+N_atom_recvRU(1)
        NLD=NRU+N_atom_recvLD(1)
        NLU=NLD+N_atom_recvLU(1)
        NRD=NLU+N_atom_recvRD(1)
        NUU=NRD+N_atom_recvUU(1)
        NDD=NUU+N_atom_recvDD(1)
        NRR3=NDD+(N_atom_recvRR(2)-N_atom_recvRR(1))
        NLL3=NRR3+(N_atom_recvLL(2)-N_atom_recvLL(1))
        NRU3=NLL3+(N_atom_recvRU(2)-N_atom_recvRU(1))
        NLD3=NRU3+(N_atom_recvLD(2)-N_atom_recvLD(1))
        NLU3=NLD3+(N_atom_recvLU(2)-N_atom_recvLU(1))
        NRD3=NLU3+(N_atom_recvRD(2)-N_atom_recvRD(1))
        NUU3=NRD3+(N_atom_recvUU(2)-N_atom_recvUU(1))
        NDD3=NUU3+(N_atom_recvDD(2)-N_atom_recvDD(1))
C Right buffer force to right node
        NN=0
        DO ii=NP_node+1,NRR
          NN=NN+1
C          IF (NN.GT.4*NBMAX1) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX1, increase Density in
C     &      parameters.inc fd_pi_force'
C            CALL write_error
C          ENDIF
          buf_force_RR(NN)=float(NA(ii))
          buf_force_RR(NN+1)=RNP_node(1,ii)
          buf_force_RR(NN+2)=RNP_node(2,ii)
          buf_force_RR(NN+3)=RNP_node(3,ii)
          NN=NN+3
        ENDDO
        DO ii=NDD+1,NRR3
          NN=NN+1
C          IF (NN.GT.4*NBMAX1) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX1, increase Density in
C     &      parameters.inc fd_pi_force'
C            CALL write_error
C          ENDIF
          buf_force_RR(NN)=float(NA(ii))
          buf_force_RR(NN+1)=RNP_node(1,ii)
          buf_force_RR(NN+2)=RNP_node(2,ii)
          buf_force_RR(NN+3)=RNP_node(3,ii)
          NN=NN+3
        ENDDO
        N_force_passRR(3)=NN/4
C        WRITE(*,*) mynode,'find_pi_force- passRR',NN/4
C Left buffer force to left node
        NN=0
        DO ii=NRR+1,NLL
          NN=NN+1
C          IF (NN.GT.4*NBMAX1) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX1, increase Density in
C     &      parameters.inc fd_pi_force'
C            CALL write_error
C          ENDIF
          buf_force_LL(NN)=float(NA(ii))
          buf_force_LL(NN+1)=RNP_node(1,ii)
          buf_force_LL(NN+2)=RNP_node(2,ii)
          buf_force_LL(NN+3)=RNP_node(3,ii)
          NN=NN+3
        ENDDO
        DO ii=NRR3+1,NLL3
          NN=NN+1
C          IF (NN.GT.4*NBMAX1) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX1, increase Density in
C     &      parameters.inc fd_pi_force'
C            CALL write_error
C          ENDIF
          buf_force_LL(NN)=float(NA(ii))
          buf_force_LL(NN+1)=RNP_node(1,ii)
          buf_force_LL(NN+2)=RNP_node(2,ii)
          buf_force_LL(NN+3)=RNP_node(3,ii)
          NN=NN+3
        ENDDO
        N_force_passLL(3)=NN/4
C        WRITE(*,*) mynode,'find_pi_force- passLL',NN/4
C Right up buffer to right up node
        NN=0
        DO ii=NLL+1,NRU
          NN=NN+1
C          IF (NN.GT.4*NBMAX3) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX3, increase Density in
C     &      parameters.inc fd_pi_force'
C            CALL write_error
C          ENDIF
          buf_force_RU(NN)=float(NA(ii))
          buf_force_RU(NN+1)=RNP_node(1,ii)
          buf_force_RU(NN+2)=RNP_node(2,ii)
          buf_force_RU(NN+3)=RNP_node(3,ii)
          NN=NN+3
        ENDDO
        DO ii=NLL3+1,NRU3
          NN=NN+1
C          IF (NN.GT.4*NBMAX2) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX2, increase Density in
C     &      parameters.inc fd_pi_force'
C            CALL write_error
C          ENDIF
          buf_force_RU(NN)=float(NA(ii))
          buf_force_RU(NN+1)=RNP_node(1,ii)
          buf_force_RU(NN+2)=RNP_node(2,ii)
          buf_force_RU(NN+3)=RNP_node(3,ii)
          NN=NN+3
        ENDDO
        N_force_passRU(3)=NN/4
C        WRITE(*,*) mynode,'find_pi_force- passRU',NN/4
C Left down buffer force to left down node
        NN=0
        DO ii=NRU+1,NLD
          NN=NN+1
C          IF (NN.GT.4*NBMAX3) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX3, increase Density in
C     &      parameters.inc fd_pi_force'
C            CALL write_error
C         ENDIF
          buf_force_LD(NN)=float(NA(ii))
          buf_force_LD(NN+1)=RNP_node(1,ii)
          buf_force_LD(NN+2)=RNP_node(2,ii)
          buf_force_LD(NN+3)=RNP_node(3,ii)
          NN=NN+3
        ENDDO
        DO ii=NRU3+1,NLD3
          NN=NN+1
C          IF (NN.GT.4*NBMAX3) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX3, increase Density in
C    &      parameters.inc fd_pi_force'
C            CALL write_error
C          ENDIF
          buf_force_LD(NN)=float(NA(ii))
          buf_force_LD(NN+1)=RNP_node(1,ii)
          buf_force_LD(NN+2)=RNP_node(2,ii)
          buf_force_LD(NN+3)=RNP_node(3,ii)
          NN=NN+3
        ENDDO
        N_force_passLD(3)=NN/4
C        WRITE(*,*) mynode,'find_pi_force- passLD',NN/4
C Left up buffer to left up node
        NN=0
        DO ii=NLD+1,NLU
          NN=NN+1
C          IF (NN.GT.4*NBMAX3) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX3, increase Density in
C     &      parameters.inc fd_pi_force'
C            CALL write_error
C          ENDIF
          buf_force_LU(NN)=float(NA(ii))
          buf_force_LU(NN+1)=RNP_node(1,ii)
          buf_force_LU(NN+2)=RNP_node(2,ii)
          buf_force_LU(NN+3)=RNP_node(3,ii)
          NN=NN+3
        ENDDO
        DO ii=NLD3+1,NLU3
          NN=NN+1
C          IF (NN.GT.4*NBMAX3) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX3, increase Density in
C     &      parameters.inc fd_pi_force'
C            CALL write_error
C          ENDIF
          buf_force_LU(NN)=float(NA(ii))
          buf_force_LU(NN+1)=RNP_node(1,ii)
          buf_force_LU(NN+2)=RNP_node(2,ii)
          buf_force_LU(NN+3)=RNP_node(3,ii)
          NN=NN+3
        ENDDO
        N_force_passLU(3)=NN/4
C        WRITE(*,*) mynode,'find_pi_force- passLU',NN/4
C Right down buffer force to right down node
        NN=0
        DO ii=NLU+1,NRD
          NN=NN+1
C          IF (NN.GT.4*NBMAX3) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX3, increase Density in
C     &      parameters.inc fd_pi_force'
C            CALL write_error
C          ENDIF
          buf_force_RD(NN)=float(NA(ii))
          buf_force_RD(NN+1)=RNP_node(1,ii)
          buf_force_RD(NN+2)=RNP_node(2,ii)
          buf_force_RD(NN+3)=RNP_node(3,ii)
          NN=NN+3
        ENDDO
        DO ii=NLU3+1,NRD3
          NN=NN+1
C          IF (NN.GT.4*NBMAX3) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX3, increase Density in
C     &      parameters.inc fd_pi_force'
C            CALL write_error
C          ENDIF
          buf_force_RD(NN)=float(NA(ii))
          buf_force_RD(NN+1)=RNP_node(1,ii)
          buf_force_RD(NN+2)=RNP_node(2,ii)
          buf_force_RD(NN+3)=RNP_node(3,ii)
          NN=NN+3
        ENDDO
        N_force_passRD(3)=NN/4
C        WRITE(*,*) mynode,'find_pi_force- passRD',NN/4
C Up buffer force to up node
        NN=0
        DO ii=NRD+1,NUU
          NN=NN+1
C          IF (NN.GT.4*NBMAX2) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX2, increase Density in
C     &      parameters.inc fd_pi_force'
C            CALL write_error
C          ENDIF
          buf_force_UU(NN)=float(NA(ii))
          buf_force_UU(NN+1)=RNP_node(1,ii)
          buf_force_UU(NN+2)=RNP_node(2,ii)
          buf_force_UU(NN+3)=RNP_node(3,ii)
          NN=NN+3
        ENDDO
        DO ii=NRD3+1,NUU3
          NN=NN+1
C          IF (NN.GT.4*NBMAX2) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX2, increase Density in
C     &      parameters.inc fd_pi_force'
C            CALL write_error
C          ENDIF
          buf_force_UU(NN)=float(NA(ii))
          buf_force_UU(NN+1)=RNP_node(1,ii)
          buf_force_UU(NN+2)=RNP_node(2,ii)
          buf_force_UU(NN+3)=RNP_node(3,ii)
          NN=NN+3
        ENDDO
        N_force_passUU(3)=NN/4
C        WRITE(*,*) mynode,'find_pi_force- passUU',NN/4
C Down buffer force to down node
        NN=0
        DO ii=NUU+1,NDD
          NN=NN+1
C          IF (NN.GT.4*NBMAX2) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX2, increase Density in
C     &      parameters.inc fd_pi_force'
C            CALL write_error
C          ENDIF
          buf_force_DD(NN)=float(NA(ii))
          buf_force_DD(NN+1)=RNP_node(1,ii)
          buf_force_DD(NN+2)=RNP_node(2,ii)
          buf_force_DD(NN+3)=RNP_node(3,ii)
          NN=NN+3
        ENDDO
        DO ii=NUU3+1,NPtot3_node
          NN=NN+1
C          IF (NN.GT.4*NBMAX2) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX2, increase Density in
C     &      parameters.inc fd_pi_force'
C            CALL write_error
C          ENDIF
          buf_force_DD(NN)=float(NA(ii))
          buf_force_DD(NN+1)=RNP_node(1,ii)
          buf_force_DD(NN+2)=RNP_node(2,ii)
          buf_force_DD(NN+3)=RNP_node(3,ii)
          NN=NN+3
        ENDDO
        N_force_passDD(3)=NN/4
C        WRITE(*,*) mynode,'find_pi_force',NP_node,NDD,NDD3
      ENDIF
C      WRITE(*,*) mynode,'find_pi_force- finish'
      RETURN
      END
