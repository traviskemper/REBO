       SUBROUTINE find_pass_force
       
C This subroutine similar to finding_neighbor2.f
C find the atoms in the inner 4 layers and make the list to pass to neighbor nodes
      USE dyn_array
      USE MPIvars
 
      IMPLICIT none
      INCLUDE 'mpif.h'


C global variable

C local variable
      INTEGER :: i,ii,NN
     & ,NATOM

CCC 1-D parallel
      IF (NNDP.EQ.1) THEN
C Find the left side atoms that will be passed to left node
        NN=0
        DO i=1,N_atom_passLL(3)
          ii=(i-1)*15+1
          NN=NN+1
C          IF (NN.GT.4*NBMAX1) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX1, increase Density in
C     &      parameters.inc fd_pass_force'
C            CALL write_error
C          ENDIF
          NATOM=INT(buf_atom_passLL(ii))
          buf_force_LL(NN)=buf_atom_passLL(ii)
          buf_force_LL(NN+1)=RNP(1,NATOM)
          buf_force_LL(NN+2)=RNP(2,NATOM)
          buf_force_LL(NN+3)=RNP(3,NATOM)
          NN=NN+3
        ENDDO
        N_force_passLL(3)=NN/4
C        DO i=1,NLL
C        IF (mynode.EQ.0) THEN
C        WRITE(33,999) move_force_L(i),buf_force_L(1,i),buf_force_L(2,i),
C     &    buf_force_L(3,i)
C        ENDIF
C        IF (mynode.EQ.1) THEN
C        WRITE(34,999) move_force_L(i),buf_force_L(1,i),buf_force_L(2,i),
C     &    buf_force_L(3,i)
C        ENDIF
C        ENDDO
C Find the right side atoms that will be passed to right node
        NN=0
        DO i=1,N_atom_passRR(3)
          ii=(i-1)*15+1
          NN=NN+1
C          IF (NN.GT.4*NBMAX1) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX1, increase Density in
C     &      parameters.inc fd_pass_force'
C            CALL write_error
C          ENDIF
          NATOM=INT(buf_atom_passRR(ii))
          buf_force_RR(NN)=buf_atom_passRR(ii)
          buf_force_RR(NN+1)=RNP(1,NATOM)
          buf_force_RR(NN+2)=RNP(2,NATOM)
          buf_force_RR(NN+3)=RNP(3,NATOM)
          NN=NN+3
        ENDDO
        N_force_passRR(3)=NN/4
C        DO i=1,NRR
C        IF (mynode.EQ.0) THEN
C        WRITE(33,999) move_force_R(i),buf_force_R(1,i),buf_force_R(2,i),
C     &    buf_force_R(3,i)
C        ENDIF
C        IF (mynode.EQ.1) THEN
C        WRITE(34,999) move_force_R(i),buf_force_R(1,i),buf_force_R(2,i),
C     &    buf_force_R(3,i)
C        ENDIF
C        ENDDO
C999   FORMAT(I5,3(E20.11))        
C        WRITE(*,*) 'find_pass_force',mynode,NLL,NRR
CCC 2-D parallel
      ELSE
C Find the left side atoms that will be passed to left node
        NN=0
        DO i=1,N_atom_passLL(3)
          ii=(i-1)*15+1
          NN=NN+1
C          IF (NN.GT.4*NBMAX1) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX1, increase Density in
C     &      parameters.inc fd_pass_force'
C            CALL write_error
C          ENDIF
          NATOM=INT(buf_atom_passLL(ii))
          buf_force_LL(NN)=buf_atom_passLL(ii)
          buf_force_LL(NN+1)=RNP(1,NATOM)
          buf_force_LL(NN+2)=RNP(2,NATOM)
          buf_force_LL(NN+3)=RNP(3,NATOM)
          NN=NN+3
        ENDDO
        N_force_passLL(3)=INT(NN/4)
C Find the right side atoms that will be passed to right node
        NN=0
        DO i=1,N_atom_passRR(3)
          ii=(i-1)*15+1
          NN=NN+1
C          IF (NN.GT.4*NBMAX1) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX1, increase Density in
C     &      parameters.inc fd_pass_force'
C            CALL write_error
C          ENDIF
          NATOM=INT(buf_atom_passRR(ii))
          buf_force_RR(NN)=buf_atom_passRR(ii)
          buf_force_RR(NN+1)=RNP(1,NATOM)
          buf_force_RR(NN+2)=RNP(2,NATOM)
          buf_force_RR(NN+3)=RNP(3,NATOM)
          NN=NN+3
        ENDDO
        N_force_passRR(3)=INT(NN/4)
C Find the DOWN side atoms that will be passed to DOWN node
        NN=0
        DO i=1,N_atom_passDD(3)
          ii=(i-1)*15+1
          NN=NN+1
C          IF (NN.GT.4*NBMAX2) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX2, increase Density in
C     &      parameters.inc fd_pass_force'
C            CALL write_error
C          ENDIF
          NATOM=INT(buf_atom_passDD(ii))
          buf_force_DD(NN)=buf_atom_passDD(ii)
          buf_force_DD(NN+1)=RNP(1,NATOM)
          buf_force_DD(NN+2)=RNP(2,NATOM)
          buf_force_DD(NN+3)=RNP(3,NATOM)
          NN=NN+3
        ENDDO
        N_force_passDD(3)=INT(NN/4)
C Find the UP side atoms that will be passed to UP node
        NN=0
        DO i=1,N_atom_passUU(3)
          ii=(i-1)*15+1
          NN=NN+1
C          IF (NN.GT.4*NBMAX2) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX2, increase Density in
C     &      parameters.inc fd_pass_force'
C            CALL write_error
C          ENDIF
          NATOM=INT(buf_atom_passUU(ii))
          buf_force_UU(NN)=buf_atom_passUU(ii)
          buf_force_UU(NN+1)=RNP(1,NATOM)
          buf_force_UU(NN+2)=RNP(2,NATOM)
          buf_force_UU(NN+3)=RNP(3,NATOM)
          NN=NN+3
        ENDDO
        N_force_passUU(3)=INT(NN/4)
C Find the LEFT DOWN side atoms that will be passed to LEFT DOWN side
        NN=0
        DO i=1,N_atom_passLD(3)
          ii=(i-1)*15+1
          NN=NN+1
C          IF (NN.GT.4*NBMAX3) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX3, increase Density in
C     &      parameters.inc fd_pass_force'
C            CALL write_error
C          ENDIF
          NATOM=INT(buf_atom_passLD(ii))
          buf_force_LD(NN)=buf_atom_passLD(ii)
          buf_force_LD(NN+1)=RNP(1,NATOM)
          buf_force_LD(NN+2)=RNP(2,NATOM)
          buf_force_LD(NN+3)=RNP(3,NATOM)
          NN=NN+3
        ENDDO
        N_force_passLD(3)=INT(NN/4)
C Find the LEFT UP side atoms that will be passed to LEFT UP side
        NN=0
        DO i=1,N_atom_passLU(3)
          ii=(i-1)*15+1
          NN=NN+1
C          IF (NN.GT.4*NBMAX3) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX3, increase Density in
C     &      parameters.inc fd_pass_force'
C            CALL write_error
C          ENDIF
          NATOM=INT(buf_atom_passLU(ii))
          buf_force_LU(NN)=buf_atom_passLU(ii)
          buf_force_LU(NN+1)=RNP(1,NATOM)
          buf_force_LU(NN+2)=RNP(2,NATOM)
          buf_force_LU(NN+3)=RNP(3,NATOM)
          NN=NN+3
        ENDDO
        N_force_passLU(3)=INT(NN/4)
C Find the RIGHT DOWN side atoms that will be passed to RIGHT DOWN side
        NN=0
        DO i=1,N_atom_passRD(3)
          ii=(i-1)*15+1
          NN=NN+1
C          IF (NN.GT.4*NBMAX3) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX3, increase Density in
C     &      parameters.inc fd_pass_force'
C            CALL write_error
C          ENDIF
          NATOM=INT(buf_atom_passRD(ii))
          buf_force_RD(NN)=buf_atom_passRD(ii)
          buf_force_RD(NN+1)=RNP(1,NATOM)
          buf_force_RD(NN+2)=RNP(2,NATOM)
          buf_force_RD(NN+3)=RNP(3,NATOM)
          NN=NN+3
        ENDDO
        N_force_passRD(3)=INT(NN/4)
C Find the RIGHT UP side atoms that will be passed to RIGHT UP node
        NN=0
        DO i=1,N_atom_passRU(3)
          ii=(i-1)*15+1
          NN=NN+1
C          IF (NN.GT.4*NBMAX3) THEN
C            WRITE(*,*) NN,'NN is GT 4*NBMAX3, increase Density in
C     &      parameters.inc fd_pass_force'
C            CALL write_error
C          ENDIF
          NATOM=INT(buf_atom_passRU(ii))
          buf_force_RU(NN)=buf_atom_passRU(ii)
          buf_force_RU(NN+1)=RNP(1,NATOM)
          buf_force_RU(NN+2)=RNP(2,NATOM)
          buf_force_RU(NN+3)=RNP(3,NATOM)
          NN=NN+3
        ENDDO
        N_force_passRU(3)=INT(NN/4)

C        WRITE(*,*) 'find_pass_force',mynode,NLL,NRR,NDD,NUU,' LL'
C        WRITE(*,*) 'find_pass_force',mynode,NLD,NLU,NRD,NRU,' LD'
      ENDIF

      RETURN
      END


