      SUBROUTINE find_center
C
C After predictor or corrector the atom's posiiton is changed
C When the atoms move accross the node region, that atom should be
C re-assign to the corresct node (only atoms in the neighbor region will move accross the node)
C And also when the atom move across the center region and neighbor region the atom also need
C to be re-assign to the correct region
C This work should do twise for pibond and Lennard Jones, because atoms are calculated seperatly
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
      USE STRUCTR

      IMPLICIT none

      INCLUDE 'mpif.h'

C global variable

C local variable
      INTEGER  cell,move_mark(NPMAX)
     &,i,j,k,II,NN,MM
      REAL*8 ncell(3),pmax_node_temp(3),pmin_node_temp(3),
     &       pmin_node_range(3),pmax_node_range(3)


      move_mark=0
C      OPEN(35,file='ii.dat',status='unknown')

CCC Check and make list for atoms going out of center region
CCC 1-D parallel
      IF (NNDP.EQ.1) THEN
        N_atom_passLL=0
        N_atom_passRR=0
C Sort the first inner shell atoms in the center region
        pmin_node_temp(NPD_i)=pmin_node(NPD_i)-1.0d-12
        pmin_node_temp(NPD_i)=pmin_node_temp(NPD_i)-CUBE(NPD_i)*
     &  ANINT(pmin_node_temp(NPD_i)/CUBE(NPD_i))
        pmin_node_range(NPD_i)=pmin_node(NPD_i)-2.0d0
        pmin_node_range(NPD_i)=pmin_node_range(NPD_i)-CUBE(NPD_i)*
     &  ANINT(pmin_node_range(NPD_i)/CUBE(NPD_i))
        pmax_node_temp(NPD_i)=pmax_node(NPD_i)+1.0d-12
        pmax_node_temp(NPD_i)=pmax_node_temp(NPD_i)-CUBE(NPD_i)*
     &  ANINT(pmax_node_temp(NPD_i)/CUBE(NPD_i))
        pmax_node_range(NPD_i)=pmax_node(NPD_i)+2.0d0
        pmax_node_range(NPD_i)=pmax_node_range(NPD_i)-CUBE(NPD_i)*
     &  ANINT(pmax_node_range(NPD_i)/CUBE(NPD_i))

C        WRITE(*,*) 'NPD_i=',NPD_i
C        WRITE(*,999)mynode,'pmin',pmin_node_temp(1),
C     &  pmin_node_temp(2),pmin_node_temp(3)
C        WRITE(*,999)mynode,'pmax',pmax_node_temp(1),
C     &  pmax_node_temp(2),pmax_node_temp(3)
C        WRITE(*,999)mynode,'pminR',pmin_node_range(1),
C     &  pmin_node_range(2),pmin_node_range(3)
C        WRITE(*,999)mynode,'pmaxR',pmax_node_range(1),
C     &  pmax_node_range(2),pmax_node_range(3)
        NLL=0
        NRR=0
        DO i=node_cell_min(1),node_cell_max(1)
          DO j=node_cell_min(2),node_cell_max(2)
            DO k=node_cell_min(3),node_cell_max(3)
              IF (NPD_i.EQ.1) THEN
                IF ((i.GT.(node_cell_min(1)+1)).AND.
     &              (i.LT.(node_cell_max(1)-1))) GOTO 100
              ELSEIF (NPD_i.EQ.2) THEN
                IF ((j.GT.(node_cell_min(2)+1)).AND.
     &              (j.LT.(node_cell_max(2)-1))) GOTO 101
              ELSE
                IF ((k.GT.(node_cell_min(3)+1)).AND.
     &              (k.LT.(node_cell_max(3)-1))) GOTO 102
              ENDIF
              cell=i+(j-1)*nncell(1)+(k-1)*nncell(2)
     &        *nncell(1)
              ii=nchead_node(cell)
              DO while (ii.NE.0)
                IF (ii.GT.NP_node) GOTO 110
C check if atom i is in the same cell
C If not in the same cell, check if the atom move out of the node then make the move atom list
                IF ((R0_node(NPD_i,ii).GE.pmin_node(NPD_i)).AND.
     &              (R0_node(NPD_i,ii).LT.pmax_node(NPD_i))) GOTO 110
                  IF ((R0_node(NPD_i,ii).LE.pmin_node_temp(NPD_i)).AND.
     &            (R0_node(NPD_i,ii).GT.pmin_node_range(NPD_i))) THEN
C these atoms are moving to left node
                    move_mark(ii)=1
                    NLL=NLL+1
                    buf_atom_passLL(NLL)=float(NA(ii))
                    buf_atom_passLL(NLL+1)=float(ktype_node(ii))
                    buf_atom_passLL(NLL+2)=float(itr_node(ii))
                    buf_atom_passLL(NLL+3)=r0_node(1,ii)
                    buf_atom_passLL(NLL+4)=r0_node(2,ii)
                    buf_atom_passLL(NLL+5)=r0_node(3,ii)
                    buf_atom_passLL(NLL+6)=r1_node(1,ii)
                    buf_atom_passLL(NLL+7)=r1_node(2,ii)
                    buf_atom_passLL(NLL+8)=r1_node(3,ii)
                    buf_atom_passLL(NLL+9)=r2_node(1,ii)
                    buf_atom_passLL(NLL+10)=r2_node(2,ii)
                    buf_atom_passLL(NLL+11)=r2_node(3,ii)
                    buf_atom_passLL(NLL+12)=r3_node(1,ii)
                    buf_atom_passLL(NLL+13)=r3_node(2,ii)
                    buf_atom_passLL(NLL+14)=r3_node(3,ii)
                    NLL=NLL+14
                  ELSEIF ((R0_node(NPD_i,ii).GE.pmax_node_temp(NPD_i))
     &            .AND.(R0_node(NPD_i,ii).LT.pmax_node_range(NPD_i)))
     &            THEN
C these atoms are moving to right node
                    move_mark(ii)=1
                    NRR=NRR+1
                    buf_atom_passRR(NRR)=float(NA(ii))
                    buf_atom_passRR(NRR+1)=float(ktype_node(ii))
                    buf_atom_passRR(NRR+2)=float(itr_node(ii))
                    buf_atom_passRR(NRR+3)=r0_node(1,ii)
                    buf_atom_passRR(NRR+4)=r0_node(2,ii)
                    buf_atom_passRR(NRR+5)=r0_node(3,ii)
                    buf_atom_passRR(NRR+6)=r1_node(1,ii)
                    buf_atom_passRR(NRR+7)=r1_node(2,ii)
                    buf_atom_passRR(NRR+8)=r1_node(3,ii)
                    buf_atom_passRR(NRR+9)=r2_node(1,ii)
                    buf_atom_passRR(NRR+10)=r2_node(2,ii)
                    buf_atom_passRR(NRR+11)=r2_node(3,ii)
                    buf_atom_passRR(NRR+12)=r3_node(1,ii)
                    buf_atom_passRR(NRR+13)=r3_node(2,ii)
                    buf_atom_passRR(NRR+14)=r3_node(3,ii)
                    NRR=NRR+14
                  ELSE
                    WRITE(*,*) mynode,' should NOT have this atom out
     &              of range',ii
                    CALL write_error 
                  ENDIF
110             ii=nclist_node(ii)
              ENDDO
102           CONTINUE
            ENDDO
101         CONTINUE
          ENDDO
100       CONTINUE
        ENDDO
        N_atom_passLL(3)=NLL/15
        N_atom_passRR(3)=NRR/15
C        WRITE(*,*) mynode,'fd_center',' LL=',NLL,' RR=',NRR
C resort the correct atoms
        NN=0
        DO i=1,NP_node
           IF (move_mark(i).EQ.0) THEN
             NN=NN+1
             IF (NN.NE.i) THEN
               NA(NN)=NA(i)
               KTYPE_node(NN)=KTYPE_node(i)
               itr_node(NN)=itr_node(i)
               DO MM=1,3
                 r0_node(MM,NN)=r0_node(MM,i)
                 r1_node(MM,NN)=r1_node(MM,i)
                 r2_node(MM,NN)=r2_node(MM,i)
                 r3_node(MM,NN)=r3_node(MM,i)
               ENDDO
             ENDIF
           ENDIF
         ENDDO
         NP_node=NN
C         WRITE(*,*) mynode,'NN=',NN,'164'
C         WRITE(*,*) 'fd_center',mynode,NP_node
CCC 2-D parallel
      ELSE
        N_atom_passLL=0
        N_atom_passRR=0
        N_atom_passLD=0
        N_atom_passRU=0
        N_atom_passRD=0
        N_atom_passLU=0
        N_atom_passDD=0
        N_atom_passUU=0
C Sort the first inner shell atoms in the center region
        pmin_node_temp(NPD_i)=pmin_node(NPD_i)-1.0d-12
        pmin_node_temp(NPD_i)=pmin_node_temp(NPD_i)-CUBE(NPD_i)*
     &  ANINT(pmin_node_temp(NPD_i)/CUBE(NPD_i))
        pmin_node_range(NPD_i)=pmin_node(NPD_i)-2.0d0
        pmin_node_range(NPD_i)=pmin_node_range(NPD_i)-CUBE(NPD_i)*
     &  ANINT(pmin_node_range(NPD_i)/CUBE(NPD_i))
        pmax_node_temp(NPD_i)=pmax_node(NPD_i)+1.0d-12
        pmax_node_temp(NPD_i)=pmax_node_temp(NPD_i)-CUBE(NPD_i)*
     &  ANINT(pmax_node_temp(NPD_i)/CUBE(NPD_i))
        pmax_node_range(NPD_i)=pmax_node(NPD_i)+2.0d0
        pmax_node_range(NPD_i)=pmax_node_range(NPD_i)-CUBE(NPD_i)*
     &  ANINT(pmax_node_range(NPD_i)/CUBE(NPD_i))

        pmin_node_temp(NPD_j)=pmin_node(NPD_j)-1.0d-12
        pmin_node_temp(NPD_j)=pmin_node_temp(NPD_j)-CUBE(NPD_j)*
     &  ANINT(pmin_node_temp(NPD_j)/CUBE(NPD_j))
        pmin_node_range(NPD_j)=pmin_node(NPD_j)-2.0d0
        pmin_node_range(NPD_j)=pmin_node_range(NPD_j)-CUBE(NPD_j)*
     &  ANINT(pmin_node_range(NPD_j)/CUBE(NPD_j))
        pmax_node_temp(NPD_j)=pmax_node(NPD_j)+1.0d-12
        pmax_node_temp(NPD_j)=pmax_node_temp(NPD_j)-CUBE(NPD_j)*
     &  ANINT(pmax_node_temp(NPD_j)/CUBE(NPD_j))
        pmax_node_range(NPD_j)=pmax_node(NPD_j)+2.0d0
        pmax_node_range(NPD_j)=pmax_node_range(NPD_j)-CUBE(NPD_j)*
     &  ANINT(pmax_node_range(NPD_j)/CUBE(NPD_j))
C        WRITE(*,*)mynode,'pmini',pmin_node_temp(NPD_i),
C     &  pmin_node_range(NPD_i)
C        WRITE(*,*)mynode,'pmaxi',pmax_node_temp(NPD_i),
C     &  pmax_node_range(NPD_i)
C        WRITE(*,*)mynode,'pminj',pmin_node_temp(NPD_j),
C     &  pmin_node_range(NPD_j)
C        WRITE(*,*)mynode,'pmaxj',pmax_node_temp(NPD_j),
C     &  pmax_node_range(NPD_j)
        NRR=0
        NLL=0
        NDD=0
        NUU=0
        NRD=0
        NRU=0
        NLD=0
        NLU=0
        DO i=node_cell_min(1),node_cell_max(1)
          DO j=node_cell_min(2),node_cell_max(2)
            DO k=node_cell_min(3),node_cell_max(3)
              IF ((NPD_i.EQ.1).AND.(NPD_j.EQ.2)) THEN
                IF ((i.GT.(node_cell_min(NPD_i)+1)).AND.
     &              (i.LT.(node_cell_max(NPD_i)-1)).AND.
     &              (j.GT.(node_cell_min(NPD_j)+1)).AND.
     &              (j.LT.(node_cell_max(NPD_j)-1))) GOTO 300
              ELSEIF ((NPD_i.EQ.1).AND.(NPD_j.EQ.3)) THEN
                IF ((i.GT.(node_cell_min(NPD_i)+1)).AND.
     &              (i.LT.(node_cell_max(NPD_i)-1)).AND.
     &              (k.GT.(node_cell_min(NPD_j)+1)).AND.
     &              (k.LT.(node_cell_max(NPD_j)-1))) GOTO 301
              ELSE
                IF ((j.GT.(node_cell_min(NPD_i)+1)).AND.
     &              (j.LT.(node_cell_max(NPD_i)-1)).AND.
     &              (k.GT.(node_cell_min(NPD_j)+1)).AND.
     &              (k.LT.(node_cell_max(NPD_j)-1))) GOTO 301
              ENDIF
              cell=i+(j-1)*nncell(1)+(k-1)*nncell(2)
     &        *nncell(1)
              ii=nchead_node(cell)

C           IF (ii.GT.0) WRITE(35,*) mynode,'ii=', ii
C -travisk Debug
C              IF (mynode.EQ.0) THEN
C                IF ( ii.GT.NPMAX) WRITE(*,*)'Error','ii=',ii
C              ENDIF

              DO while (ii.NE.0)
                IF (ii.GT.NP_node) GOTO 330
C check if atom i is in the same cell
C If not in the same cell, check if the atom move out of the node then make the move atom list
                IF ((R0_node(NPD_i,ii).GE.pmin_node(NPD_i)).AND.
     &              (R0_node(NPD_i,ii).LT.pmax_node(NPD_i)).AND.
     &              (R0_node(NPD_j,ii).GE.pmin_node(NPD_j)).AND.
     &              (R0_node(NPD_j,ii).LT.pmax_node(NPD_j))) GOTO 330
                  IF ((R0_node(NPD_i,ii).LT.pmin_node_temp(NPD_i)).AND.
     &            (R0_node(NPD_i,ii).GT.pmin_node_range(NPD_i))) THEN
                    IF ((R0_node(NPD_j,ii).GE.pmin_node(NPD_j))
     &              .AND.(R0_node(NPD_j,ii).LT.pmax_node(NPD_j)))
     &              THEN
C these atoms are moving to left node
                      move_mark(ii)=1
                      NLL=NLL+1
                      buf_atom_passLL(NLL)=float(NA(ii))
                      buf_atom_passLL(NLL+1)=float(ktype_node(ii))
                      buf_atom_passLL(NLL+2)=float(itr_node(ii))
                      buf_atom_passLL(NLL+3)=r0_node(1,ii)
                      buf_atom_passLL(NLL+4)=r0_node(2,ii)
                      buf_atom_passLL(NLL+5)=r0_node(3,ii)
                      buf_atom_passLL(NLL+6)=r1_node(1,ii)
                      buf_atom_passLL(NLL+7)=r1_node(2,ii)
                      buf_atom_passLL(NLL+8)=r1_node(3,ii)
                      buf_atom_passLL(NLL+9)=r2_node(1,ii)
                      buf_atom_passLL(NLL+10)=r2_node(2,ii)
                      buf_atom_passLL(NLL+11)=r2_node(3,ii)
                      buf_atom_passLL(NLL+12)=r3_node(1,ii)
                      buf_atom_passLL(NLL+13)=r3_node(2,ii)
                      buf_atom_passLL(NLL+14)=r3_node(3,ii)
                      NLL=NLL+14
C                    WRITE(*,*) mynode,'fd_center_L',move_atom_L(NLL),
C     &              move_ktype_L(NLL),icell,cellx,celly,cellz
                    ELSEIF((R0_node(NPD_j,ii).LE.pmin_node(NPD_j)))
     &              THEN
C these atoms are moving left down node
                      move_mark(ii)=1
                      NLD=NLD+1
                      buf_atom_passLD(NLD)=float(NA(ii))
                      buf_atom_passLD(NLD+1)=float(ktype_node(ii))
                      buf_atom_passLD(NLD+2)=float(itr_node(ii))
                      buf_atom_passLD(NLD+3)=r0_node(1,ii)
                      buf_atom_passLD(NLD+4)=r0_node(2,ii)
                      buf_atom_passLD(NLD+5)=r0_node(3,ii)
                      buf_atom_passLD(NLD+6)=r1_node(1,ii)
                      buf_atom_passLD(NLD+7)=r1_node(2,ii)
                      buf_atom_passLD(NLD+8)=r1_node(3,ii)
                      buf_atom_passLD(NLD+9)=r2_node(1,ii)
                      buf_atom_passLD(NLD+10)=r2_node(2,ii)
                      buf_atom_passLD(NLD+11)=r2_node(3,ii)
                      buf_atom_passLD(NLD+12)=r3_node(1,ii)
                      buf_atom_passLD(NLD+13)=r3_node(2,ii)
                      buf_atom_passLD(NLD+14)=r3_node(3,ii)
                      NLD=NLD+14
                    ELSEIF((R0_node(NPD_j,ii).GE.pmax_node(NPD_j)))
     &              THEN
C these atoms are moving left up node
                      move_mark(ii)=1
                      NLU=NLU+1
                      buf_atom_passLU(NLU)=float(NA(ii))
                      buf_atom_passLU(NLU+1)=float(ktype_node(ii))
                      buf_atom_passLU(NLU+2)=float(itr_node(ii))
                      buf_atom_passLU(NLU+3)=r0_node(1,ii)
                      buf_atom_passLU(NLU+4)=r0_node(2,ii)
                      buf_atom_passLU(NLU+5)=r0_node(3,ii)
                      buf_atom_passLU(NLU+6)=r1_node(1,ii)
                      buf_atom_passLU(NLU+7)=r1_node(2,ii)
                      buf_atom_passLU(NLU+8)=r1_node(3,ii)
                      buf_atom_passLU(NLU+9)=r2_node(1,ii)
                      buf_atom_passLU(NLU+10)=r2_node(2,ii)
                      buf_atom_passLU(NLU+11)=r2_node(3,ii)
                      buf_atom_passLU(NLU+12)=r3_node(1,ii)
                      buf_atom_passLU(NLU+13)=r3_node(2,ii)
                      buf_atom_passLU(NLU+14)=r3_node(3,ii)
                      NLU=NLU+14
                    ELSE
                      WRITE(*,*) mynode,' should NOT have this atom:
     &                LL,LD,LU',i
                    ENDIF
                  ELSEIF ((R0_node(NPD_i,ii).GE.pmax_node_temp(NPD_i))
     &            .AND.(R0_node(NPD_i,ii).LT.pmax_node_range(NPD_i)))
     &            THEN
                    IF ((R0_node(NPD_j,ii).GE.pmin_node(NPD_j))
     &              .AND.(R0_node(NPD_j,ii).LT.pmax_node(NPD_j)))
     &              THEN
C these atoms are moving to right node
                      move_mark(ii)=1
                      NRR=NRR+1
                      buf_atom_passRR(NRR)=float(NA(ii))
                      buf_atom_passRR(NRR+1)=float(ktype_node(ii))
                      buf_atom_passRR(NRR+2)=float(itr_node(ii))
                      buf_atom_passRR(NRR+3)=r0_node(1,ii)
                      buf_atom_passRR(NRR+4)=r0_node(2,ii)
                      buf_atom_passRR(NRR+5)=r0_node(3,ii)
                      buf_atom_passRR(NRR+6)=r1_node(1,ii)
                      buf_atom_passRR(NRR+7)=r1_node(2,ii)
                      buf_atom_passRR(NRR+8)=r1_node(3,ii)
                      buf_atom_passRR(NRR+9)=r2_node(1,ii)
                      buf_atom_passRR(NRR+10)=r2_node(2,ii)
                      buf_atom_passRR(NRR+11)=r2_node(3,ii)
                      buf_atom_passRR(NRR+12)=r3_node(1,ii)
                      buf_atom_passRR(NRR+13)=r3_node(2,ii)
                      buf_atom_passRR(NRR+14)=r3_node(3,ii)
                      NRR=NRR+14
                    ELSEIF((R0_node(NPD_j,ii).LE.pmin_node(NPD_j)))
     &              THEN
C these atoms are moving to right down node
                      move_mark(ii)=1
                      NRD=NRD+1
                      buf_atom_passRD(NRD)=float(NA(ii))
                      buf_atom_passRD(NRD+1)=float(ktype_node(ii))
                      buf_atom_passRD(NRD+2)=float(itr_node(ii))
                      buf_atom_passRD(NRD+3)=r0_node(1,ii)
                      buf_atom_passRD(NRD+4)=r0_node(2,ii)
                      buf_atom_passRD(NRD+5)=r0_node(3,ii)
                      buf_atom_passRD(NRD+6)=r1_node(1,ii)
                      buf_atom_passRD(NRD+7)=r1_node(2,ii)
                      buf_atom_passRD(NRD+8)=r1_node(3,ii)
                      buf_atom_passRD(NRD+9)=r2_node(1,ii)
                      buf_atom_passRD(NRD+10)=r2_node(2,ii)
                      buf_atom_passRD(NRD+11)=r2_node(3,ii)
                      buf_atom_passRD(NRD+12)=r3_node(1,ii)
                      buf_atom_passRD(NRD+13)=r3_node(2,ii)
                      buf_atom_passRD(NRD+14)=r3_node(3,ii)
                      NRD=NRD+14
                    ELSEIF((R0_node(NPD_j,ii).GE.pmax_node(NPD_j)))
     &              THEN
C these atoms are moving to right up node
                      move_mark(ii)=1
                      NRU=NRU+1
                      buf_atom_passRU(NRU)=float(NA(ii))
                      buf_atom_passRU(NRU+1)=float(ktype_node(ii))
                      buf_atom_passRU(NRU+2)=float(itr_node(ii))
                      buf_atom_passRU(NRU+3)=r0_node(1,ii)
                      buf_atom_passRU(NRU+4)=r0_node(2,ii)
                      buf_atom_passRU(NRU+5)=r0_node(3,ii)
                      buf_atom_passRU(NRU+6)=r1_node(1,ii)
                      buf_atom_passRU(NRU+7)=r1_node(2,ii)
                      buf_atom_passRU(NRU+8)=r1_node(3,ii)
                      buf_atom_passRU(NRU+9)=r2_node(1,ii)
                      buf_atom_passRU(NRU+10)=r2_node(2,ii)
                      buf_atom_passRU(NRU+11)=r2_node(3,ii)
                      buf_atom_passRU(NRU+12)=r3_node(1,ii)
                      buf_atom_passRU(NRU+13)=r3_node(2,ii)
                      buf_atom_passRU(NRU+14)=r3_node(3,ii)
                      NRU=NRU+14
                    ELSE
                      WRITE(*,*) mynode,' should NOT have this atom:
     &                RR,RD,RU',i
                    ENDIF
                  ELSEIF ((R0_node(NPD_j,ii).LT.pmin_node_temp(NPD_j))
     &            .AND.(R0_node(NPD_j,ii).GT.pmin_node_range(NPD_j)))
     &            THEN
                    IF ((R0_node(NPD_i,ii).GE.pmin_node(NPD_i))
     &              .AND.(R0_node(NPD_i,ii).LT.pmax_node(NPD_i)))
     &              THEN
C these atoms are moving to down node
                      move_mark(ii)=1
                      NDD=NDD+1
                      buf_atom_passDD(NDD)=float(NA(ii))
                      buf_atom_passDD(NDD+1)=float(ktype_node(ii))
                      buf_atom_passDD(NDD+2)=float(itr_node(ii))
                      buf_atom_passDD(NDD+3)=r0_node(1,ii)
                      buf_atom_passDD(NDD+4)=r0_node(2,ii)
                      buf_atom_passDD(NDD+5)=r0_node(3,ii)
                      buf_atom_passDD(NDD+6)=r1_node(1,ii)
                      buf_atom_passDD(NDD+7)=r1_node(2,ii)
                      buf_atom_passDD(NDD+8)=r1_node(3,ii)
                      buf_atom_passDD(NDD+9)=r2_node(1,ii)
                      buf_atom_passDD(NDD+10)=r2_node(2,ii)
                      buf_atom_passDD(NDD+11)=r2_node(3,ii)
                      buf_atom_passDD(NDD+12)=r3_node(1,ii)
                      buf_atom_passDD(NDD+13)=r3_node(2,ii)
                      buf_atom_passDD(NDD+14)=r3_node(3,ii)
                      NDD=NDD+14
                    ELSE
                      WRITE(*,*) mynode,' should NOT have this atom:
     &                DD',NA(ii)
                WRITE(*,*)mynode,NA(ii),r0_node(1,ii),r0_node(2,ii)
     &          ,r0_node(3,ii),'DD'
                    ENDIF
                  ELSEIF ((R0_node(NPD_j,ii).GE.pmax_node_temp(NPD_j))
     &            .AND.(R0_node(NPD_j,ii).LT.pmax_node_range(NPD_j)))
     &            THEN
                    IF ((R0_node(NPD_i,ii).GE.pmin_node(NPD_i))
     &              .AND.(R0_node(NPD_i,ii).LT.pmax_node(NPD_i)))
     &              THEN
C these atoms are moving to up node
                      move_mark(ii)=1
                      NUU=NUU+1
                      buf_atom_passUU(NUU)=float(NA(ii))
                      buf_atom_passUU(NUU+1)=float(ktype_node(ii))
                      buf_atom_passUU(NUU+2)=float(itr_node(ii))
                      buf_atom_passUU(NUU+3)=r0_node(1,ii)
                      buf_atom_passUU(NUU+4)=r0_node(2,ii)
                      buf_atom_passUU(NUU+5)=r0_node(3,ii)
                      buf_atom_passUU(NUU+6)=r1_node(1,ii)
                      buf_atom_passUU(NUU+7)=r1_node(2,ii)
                      buf_atom_passUU(NUU+8)=r1_node(3,ii)
                      buf_atom_passUU(NUU+9)=r2_node(1,ii)
                      buf_atom_passUU(NUU+10)=r2_node(2,ii)
                      buf_atom_passUU(NUU+11)=r2_node(3,ii)
                      buf_atom_passUU(NUU+12)=r3_node(1,ii)
                      buf_atom_passUU(NUU+13)=r3_node(2,ii)
                      buf_atom_passUU(NUU+14)=r3_node(3,ii)
                      NUU=NUU+14
                    ELSE
                      WRITE(*,*) mynode,' should NOT have this atom:
     &                UU',NA(ii)
                WRITE(*,*)mynode,NA(ii),r0_node(1,ii),r0_node(2,ii),
     &          r0_node(3,ii),'UU'
                    ENDIF
                  ELSE
                    WRITE(*,*) mynode,' should NOT have this atom:
     &              out of range',NA(ii)
                WRITE(*,*)mynode,NA(ii),r0_node(1,ii),r0_node(2,ii),
     &          r0_node(3,ii),'UU'
                    CALL write_error
                  ENDIF
330             ii=nclist_node(ii)
              ENDDO
301           CONTINUE
            ENDDO
300       CONTINUE
          ENDDO
        ENDDO
        N_atom_passLL(3)=NLL/15
        N_atom_passRR(3)=NRR/15
        N_atom_passLD(3)=NLD/15
        N_atom_passRU(3)=NRU/15
        N_atom_passRD(3)=NRD/15
        N_atom_passLU(3)=NLU/15
        N_atom_passDD(3)=NDD/15
        N_atom_passUU(3)=NUU/15
C        WRITE(*,*) mynode,'fd_center',NRR,NLL,NDD,NUU,' RR'
C        WRITE(*,*) mynode,'fd_center',NRU,NRD,NLU,NLD,' RU'
C resort the correct atoms
        NN=0
C        WRITE(*,*) mynode,'NP_node=',NP_node,'491'
        DO i=1,NP_node
C           IF (I.GT.0) WRITE(35,*) mynode,'i=',i
           IF (move_mark(i).EQ.0) THEN
             NN=NN+1
             IF (NN.NE.i) THEN
               NA(NN)=NA(i)
               KTYPE_node(NN)=KTYPE_node(i)
               itr_node(NN)=itr_node(i)
               DO MM=1,3
                 r0_node(MM,NN)=r0_node(MM,i)
                 r1_node(MM,NN)=r1_node(MM,i)
                 r2_node(MM,NN)=r2_node(MM,i)
                 r3_node(MM,NN)=r3_node(MM,i)
               ENDDO
             ENDIF
           ENDIF
         ENDDO
C        WRITE(*,*) mynode,'NP_node=',NP_node,'509',NN
         NP_node=NN
C        WRITE(*,*) mynode,'NP_node=',NP_node,'511',NN
      ENDIF
999   FORMAT(I4,A8,3F20.11)
      RETURN
      END SUBROUTINE find_center
