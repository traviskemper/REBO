      SUBROUTINE pass_atom_node
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
C This subroutine pass the atom for pi bond

      USE dyn_array
      USE MPIvars
      USE STRUCTR

      IMPLICIT none

      INCLUDE 'mpif.h'


C global variable

C local variable
      REAL*8 ncell(3)
      INTEGER cellx,celly,cellz,Nadd
     &,i,j,NN,icell,N_data

C -travisk
       INTEGER :: NPDBTOT,ierr

C Find and check the center region atoms that need to be passed
      CALL find_center
C Pass the center region atoms that need to be passed
      CALL pass_node
C Unpack the received center region atoms on every node
      IF (NNDP.EQ.1) THEN
        N_data=15*(N_atom_recvRR(3)+N_atom_recvLL(3))
      ELSE
        N_data=15*(N_atom_recvRR(3)+N_atom_recvLL(3)+N_atom_recvDD(3)
     &            +N_atom_recvUU(3)+N_atom_recvRU(3)+N_atom_recvLU(3)
     &            +N_atom_recvRD(3)+N_atom_recvLD(3))
      ENDIF
C      WRITE(*,*) mynode,'pass_atom_node-N_data',N_data
        NN=0
      DO i=NP_node+1,NP_node+(N_data/15)
        NN=NN+1
        NA(i)=INT(buf_atom_recv(NN))
        ktype_node(i)=INT(buf_atom_recv(NN+1))
        itr_node(i)=INT(buf_atom_recv(NN+2))
        r0_node(1,i)=buf_atom_recv(NN+3)
        r0_node(2,i)=buf_atom_recv(NN+4)
        r0_node(3,i)=buf_atom_recv(NN+5)
        r1_node(1,i)=buf_atom_recv(NN+6)
        r1_node(2,i)=buf_atom_recv(NN+7)
        r1_node(3,i)=buf_atom_recv(NN+8)
        r2_node(1,i)=buf_atom_recv(NN+9)
        r2_node(2,i)=buf_atom_recv(NN+10)
        r2_node(3,i)=buf_atom_recv(NN+11)
        r3_node(1,i)=buf_atom_recv(NN+12)
        r3_node(2,i)=buf_atom_recv(NN+13)
        r3_node(3,i)=buf_atom_recv(NN+14)
        NN=NN+14
      ENDDO
      NP_node=NP_node+INT(N_data/15)
C      WRITE(*,*) 'pass_atom_node',mynode,NP_node
C Resort the center region atoms into cells
       DO icell=1,nncell(1)*nncell(2)*nncell(3)
         nchead_node(icell)=0
       ENDDO
       DO i=1,3
         ncell(i)=float(nncell(i))
       ENDDO
       DO i=1,NP_node
         icell=1+INT(((r0_node(1,i)/(CUBE(1)+1.0d-12))+0.5)*ncell(1))
     &          +INT(((r0_node(2,i)/(CUBE(2)+1.0d-12))+0.5)*ncell(2))
     &          *nncell(1)
     &          +INT(((r0_node(3,i)/(CUBE(3)+1.0d-12))+0.5)*ncell(3))
     &          *nncell(2)*nncell(1)
         nclist_node(i)=nchead_node(icell)
         nchead_node(icell)=i
       ENDDO
C Find the 1st and 2nd buffer region atoms
      CALL find_neighbor2
C Find the 3rd buffer region atoms
      CALL find_neighbor3 
C Find the 4th buffer region atoms
      CALL find_neighbor4
C Pass the all buffer region atoms
      CALL pass_node
C Unpack the buffer region atoms on every node
CCC 1-D parallel 
      IF (NNDP.EQ.1) THEN
CC Unpack 1st and 2nd buffer
C Unpack right 1st and 2nd buffer
        NN=0
        DO i=NP_node+1,NP_node+N_atom_recvRR(1)
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=NP_node+N_atom_recvRR(1)
C Unpack left 1st and 2nd buffer
        NN=15*N_atom_recvRR(3)
        DO i=Nadd+1,Nadd+N_atom_recvLL(1)
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+N_atom_recvLL(1)
        NPtot_node=Nadd
CC Unpack 3rd buffer
C Unpack 3rd right buffer
        NN=15*N_atom_recvRR(1)
        DO i=Nadd+1,Nadd+(N_atom_recvRR(2)-N_atom_recvRR(1))
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+(N_atom_recvRR(2)-N_atom_recvRR(1))
C Unpack 3rd left buffer
        NN=15*(N_atom_recvRR(3)+N_atom_recvLL(1))
        DO i=Nadd+1,Nadd+(N_atom_recvLL(2)-N_atom_recvLL(1))
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+(N_atom_recvLL(2)-N_atom_recvLL(1))
        NPtot3_node=Nadd
CC Unpack 4th buffer
C Unpack 4th right buffer
        NN=15*N_atom_recvRR(2)
        DO i=Nadd+1,Nadd+(N_atom_recvRR(3)-N_atom_recvRR(2))
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+(N_atom_recvRR(3)-N_atom_recvRR(2))
C Unpack 4th left buffer
        NN=15*(N_atom_recvRR(3)+N_atom_recvLL(2))
        DO i=Nadd+1,Nadd+(N_atom_recvLL(3)-N_atom_recvLL(2))
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+(N_atom_recvLL(3)-N_atom_recvLL(2))
        NPtot4_node=Nadd
C        WRITE(*,*) mynode,NP_node,NPtot_node,NPtot3_node,NPtot4_node
C 2-D parallel
      ELSE
CC Unpack 1st and 2nd buffer
C Unpack right 1st and 2nd buffer
        NN=0
        DO i=NP_node+1,NP_node+N_atom_recvRR(1)
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=NP_node+N_atom_recvRR(1)
C Unpack left 1st and 2nd buffer
        NN=15*N_atom_recvRR(3)
        DO i=Nadd+1,Nadd+N_atom_recvLL(1)
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+N_atom_recvLL(1)
C Unpack right-up 1st and 2nd buffer
        NN=15*(N_atom_recvRR(3)+N_atom_recvLL(3))
        DO i=Nadd+1,Nadd+N_atom_recvRU(1)
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+N_atom_recvRU(1)
C Unpack left-down 1st and 2nd buffer
        NN=15*(N_atom_recvRR(3)+N_atom_recvLL(3)+N_atom_recvRU(3))
        DO i=Nadd+1,Nadd+N_atom_recvLD(1)
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+N_atom_recvLD(1)
C Unpack left-up 1st and 2nd buffer
        NN=15*(N_atom_recvRR(3)+N_atom_recvLL(3)+N_atom_recvRU(3)+
     &         N_atom_recvLD(3))
        DO i=Nadd+1,Nadd+N_atom_recvLU(1)
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+N_atom_recvLU(1)
C Unpack right-down 1st and 2nd buffer
        NN=15*(N_atom_recvRR(3)+N_atom_recvLL(3)+N_atom_recvRU(3)+
     &         N_atom_recvLD(3)+N_atom_recvLU(3))
        DO i=Nadd+1,Nadd+N_atom_recvRD(1)
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+N_atom_recvRD(1)
C Unpack up 1st and 2nd buffer
        NN=15*(N_atom_recvRR(3)+N_atom_recvLL(3)+N_atom_recvRU(3)+
     &         N_atom_recvLD(3)+N_atom_recvLU(3)+N_atom_recvRD(3))
        DO i=Nadd+1,Nadd+N_atom_recvUU(1)
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+N_atom_recvUU(1)
C Unpack down 1st and 2nd buffer
        NN=15*(N_atom_recvRR(3)+N_atom_recvLL(3)+N_atom_recvRU(3)+
     &         N_atom_recvLD(3)+N_atom_recvLU(3)+N_atom_recvRD(3)+
     &         N_atom_recvUU(3))
        DO i=Nadd+1,Nadd+N_atom_recvDD(1)
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+N_atom_recvDD(1)
        NPtot_node=Nadd
CC Unpack 3rd buffer
C Unpack right 3rd buffer
        NN=15*N_atom_recvRR(1)
        DO i=Nadd+1,Nadd+(N_atom_recvRR(2)-N_atom_recvRR(1))
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+(N_atom_recvRR(2)-N_atom_recvRR(1))
C Unpack left 3rd buffer
        NN=15*(N_atom_recvRR(3)+N_atom_recvLL(1))
        DO i=Nadd+1,Nadd+(N_atom_recvLL(2)-N_atom_recvLL(1))
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+(N_atom_recvLL(2)-N_atom_recvLL(1))
C Unpack right-up 3rd buffer
        NN=15*(N_atom_recvRR(3)+N_atom_recvLL(3)+N_atom_recvRU(1))
        DO i=Nadd+1,Nadd+(N_atom_recvRU(2)-N_atom_recvRU(1))
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+(N_atom_recvRU(2)-N_atom_recvRU(1))
C Unpack left-down 3rd buffer
        NN=15*(N_atom_recvRR(3)+N_atom_recvLL(3)+N_atom_recvRU(3)+
     &         N_atom_recvLD(1))
        DO i=Nadd+1,Nadd+(N_atom_recvLD(2)-N_atom_recvLD(1))
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+(N_atom_recvLD(2)-N_atom_recvLD(1))
C Unpack left-up 3rd buffer
        NN=15*(N_atom_recvRR(3)+N_atom_recvLL(3)+N_atom_recvRU(3)+
     &         N_atom_recvLD(3)+N_atom_recvLU(1))
        DO i=Nadd+1,Nadd+(N_atom_recvLU(2)-N_atom_recvLU(1))
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+(N_atom_recvLU(2)-N_atom_recvLU(1))
C Unpack right-down 3rd buffer
        NN=15*(N_atom_recvRR(3)+N_atom_recvLL(3)+N_atom_recvRU(3)+
     &         N_atom_recvLD(3)+N_atom_recvLU(3)+N_atom_recvRD(1))
        DO i=Nadd+1,Nadd+(N_atom_recvRD(2)-N_atom_recvRD(1))
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+(N_atom_recvRD(2)-N_atom_recvRD(1))
C Unpack up 3rd buffer
        NN=15*(N_atom_recvRR(3)+N_atom_recvLL(3)+N_atom_recvRU(3)+
     &         N_atom_recvLD(3)+N_atom_recvLU(3)+N_atom_recvRD(3)+
     &         N_atom_recvUU(1))
        DO i=Nadd+1,Nadd+(N_atom_recvUU(2)-N_atom_recvUU(1))
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+(N_atom_recvUU(2)-N_atom_recvUU(1))
C Unpack down 3rd buffer
        NN=15*(N_atom_recvRR(3)+N_atom_recvLL(3)+N_atom_recvRU(3)+
     &         N_atom_recvLD(3)+N_atom_recvLU(3)+N_atom_recvRD(3)+
     &         N_atom_recvUU(3)+N_atom_recvDD(1))
        DO i=Nadd+1,Nadd+(N_atom_recvDD(2)-N_atom_recvDD(1))
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+(N_atom_recvDD(2)-N_atom_recvDD(1))
        NPtot3_node=Nadd
CC Unpack 4th buffer
C Unpack right 4th buffer
        NN=15*N_atom_recvRR(2)
        DO i=Nadd+1,Nadd+(N_atom_recvRR(3)-N_atom_recvRR(2))
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+(N_atom_recvRR(3)-N_atom_recvRR(2))
C Unpack left 4th buffer
        NN=15*(N_atom_recvRR(3)+N_atom_recvLL(2))
        DO i=Nadd+1,Nadd+(N_atom_recvLL(3)-N_atom_recvLL(2))
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+(N_atom_recvLL(3)-N_atom_recvLL(2))
C Unpack right-up 4th buffer
        NN=15*(N_atom_recvRR(3)+N_atom_recvLL(3)+N_atom_recvRU(2))
        DO i=Nadd+1,Nadd+(N_atom_recvRU(3)-N_atom_recvRU(2))
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+(N_atom_recvRU(3)-N_atom_recvRU(2))
C Unpack left-down 4th buffer
        NN=15*(N_atom_recvRR(3)+N_atom_recvLL(3)+N_atom_recvRU(3)+
     &         N_atom_recvLD(2))
        DO i=Nadd+1,Nadd+(N_atom_recvLD(3)-N_atom_recvLD(2))
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+(N_atom_recvLD(3)-N_atom_recvLD(2))
C Unpack left-up 4th buffer
        NN=15*(N_atom_recvRR(3)+N_atom_recvLL(3)+N_atom_recvRU(3)+
     &         N_atom_recvLD(3)+N_atom_recvLU(2))
        DO i=Nadd+1,Nadd+(N_atom_recvLU(3)-N_atom_recvLU(2))
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+(N_atom_recvLU(3)-N_atom_recvLU(2))
C Unpack right-down 4th buffer
        NN=15*(N_atom_recvRR(3)+N_atom_recvLL(3)+N_atom_recvRU(3)+
     &         N_atom_recvLD(3)+N_atom_recvLU(3)+N_atom_recvRD(2))
        DO i=Nadd+1,Nadd+(N_atom_recvRD(3)-N_atom_recvRD(2))
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+(N_atom_recvRD(3)-N_atom_recvRD(2))
C Unpack up 4th buffer
        NN=15*(N_atom_recvRR(3)+N_atom_recvLL(3)+N_atom_recvRU(3)+
     &         N_atom_recvLD(3)+N_atom_recvLU(3)+N_atom_recvRD(3)+
     &         N_atom_recvUU(2))
        DO i=Nadd+1,Nadd+(N_atom_recvUU(3)-N_atom_recvUU(2))
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+(N_atom_recvUU(3)-N_atom_recvUU(2))
C Unpack down 4th buffer
        NN=15*(N_atom_recvRR(3)+N_atom_recvLL(3)+N_atom_recvRU(3)+
     &         N_atom_recvLD(3)+N_atom_recvLU(3)+N_atom_recvRD(3)+
     &         N_atom_recvUU(3)+N_atom_recvDD(2))
        DO i=Nadd+1,Nadd+(N_atom_recvDD(3)-N_atom_recvDD(2))
          NN=NN+1
          NA(i)=INT(buf_atom_recv(NN))
          ktype_node(i)=INT(buf_atom_recv(NN+1))
          itr_node(i)=INT(buf_atom_recv(NN+2))
          r0_node(1,i)=buf_atom_recv(NN+3)
          r0_node(2,i)=buf_atom_recv(NN+4)
          r0_node(3,i)=buf_atom_recv(NN+5)
          r1_node(1,i)=buf_atom_recv(NN+6)
          r1_node(2,i)=buf_atom_recv(NN+7)
          r1_node(3,i)=buf_atom_recv(NN+8)
          r2_node(1,i)=buf_atom_recv(NN+9)
          r2_node(2,i)=buf_atom_recv(NN+10)
          r2_node(3,i)=buf_atom_recv(NN+11)
          r3_node(1,i)=buf_atom_recv(NN+12)
          r3_node(2,i)=buf_atom_recv(NN+13)
          r3_node(3,i)=buf_atom_recv(NN+14)
          NN=NN+14
        ENDDO
        Nadd=Nadd+(N_atom_recvDD(3)-N_atom_recvDD(2))
        NPtot4_node=Nadd
      ENDIF
      IF (NPMAX.LT.NPtot4_node) THEN
        WRITE(*,*) 'NPMAX is too small',NPMAX,NPtot4_node
        WRITE(*,*) 'increase NPMAX in allocate_arrays.f'
        CALL write_error
      ENDIF
C     WRITE(*,*) mynode,'pass_atom_node',NP_node,NPtot_node,
C    & NPtot3_node,NPtot4_node
C UPDATE NEIGHBOR LIST
C
C -travisk
c$$$  
c$$$      NPDBTOT = 0
c$$$      
c$$$      DO 7 I=1,NP_node
c$$$           DO 6 J=1,3
c$$$                R0L(J,I)=R0_node(J,I)
c$$$C -travisk
c$$$           NPDBTOT = NPDBTOT + 1
c$$$6          CONTINUE
c$$$7     CONTINUE
C      IF (mynode.EQ.0 ) WRITE(*,*) 'POOP1' 
C      CALL MPI_REDUCE(NPDBTOT,NP_node,1,MPI_INTEGER,MPI_SUM,0,
C     &   MPI_COMM_WORLD,ierr)
C      IF (mynode.EQ.0 ) THEN
C        WRITE(*,*) mynode,'total atoms',NPDBTOT
C      ENDIF
C      IF ( NPDBTOT.NE.NP )THEN
C        WRITE(*,*) 'Atoms have gone missing'
C        CALL write_error
C      ENDIF

      RETURN

      END SUBROUTINE pass_atom_node
