      SUBROUTINE find_neighbor2
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
C This subroutine find the neighbor region atoms for Pi
      USE MPIvars
      USE dyn_array
      USE STRUCTR

      IMPLICIT none
      INCLUDE 'mpif.h'

C global variable

C local variable
      INTEGER  icell,cell
      INTEGER  node_cell_min_temp(3),node_cell_max_temp(3)
     &,i,j,k,II,NN
      REAL*8 R_i,R_j

CCC pass_atom
CCC 1-D parallel
      IF (NNDP.EQ.1) THEN
        NLL=0
        NRR=0
        node_cell_min_temp=node_cell_min
        node_cell_max_temp=node_cell_max
        node_cell_max_temp(NPD_i)=node_cell_min(NPD_i)+2
C Find the left side atoms that will be passed to left node
        DO i=node_cell_min_temp(1),node_cell_max_temp(1)
          DO j=node_cell_min_temp(2),node_cell_max_temp(2)
            DO k=node_cell_min_temp(3),node_cell_max_temp(3)
              cell=i+(j-1)*nncell(1)+(k-1)*nncell(2)
     &        *nncell(1)
              ii=nchead_node(cell)
              DO while (ii.NE.0)
                R_i=ABS(R0_node(NPD_i,ii)-pmin_node(NPD_i))
                IF (R_i.LE.2*RLIMAX) THEN
C these atoms are moving to left node
                  NLL=NLL+1
C                  IF (NLL.GT.15*NBMAX1) THEN
C                    WRITE(*,*) NLL,' NLL is GT 15*NBMAX1,
C     &              increase Density in parameters.inc fd_nabor2'
C                    STOP
C                  ENDIF
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
                ENDIF
                ii=nclist_node(ii)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        node_cell_min_temp=node_cell_min
        node_cell_max_temp=node_cell_max
        node_cell_min_temp(NPD_i)=node_cell_max(NPD_i)-2
C Find the right side atoms that will be passed to right node
        DO i=node_cell_min_temp(1),node_cell_max_temp(1)
          DO j=node_cell_min_temp(2),node_cell_max_temp(2)
            DO k=node_cell_min_temp(3),node_cell_max_temp(3)
              cell=i+(j-1)*nncell(1)+(k-1)*nncell(2)
     &        *nncell(1)
              ii=nchead_node(cell)
              DO while (ii.NE.0)
                R_i=ABS(R0_node(NPD_i,ii)-pmax_node(NPD_i))
                IF (R_i.LE.2*RLIMAX) THEN
C these atoms are moving to RIGHT node
                  NRR=NRR+1
C                  IF (NRR.GT.15*NBMAX1) THEN
C                    WRITE(*,*) NRR,' NRR is GT 15*NBMAX1,
C     &              increase Density in parameters.inc fd_nabor2'
C                    STOP
C                  ENDIF
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
                ENDIF
                ii=nclist_node(ii)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        N_atom_passLL(1)=NLL/15
        N_atom_passRR(1)=NRR/15
c        WRITE(*,*) 'find_nabor2,NLL,NRR,',mynode,NLL/15,NRR/15,NN
CCC 2-D parallel
      ELSE
        NLL=0
        NRR=0
        NDD=0
        NUU=0
        NLD=0
        NLU=0
        NRD=0
        NRU=0
        node_cell_min_temp=node_cell_min
        node_cell_max_temp=node_cell_max
        node_cell_max_temp(NPD_i)=node_cell_min(NPD_i)+2
C Find the left side atoms that will be passed to left node
        DO i=node_cell_min_temp(1),node_cell_max_temp(1)
          DO j=node_cell_min_temp(2),node_cell_max_temp(2)
            DO k=node_cell_min_temp(3),node_cell_max_temp(3)
              cell=i+(j-1)*nncell(1)+(k-1)*nncell(2)
     &        *nncell(1)
              ii=nchead_node(cell)
              DO while (ii.NE.0)
                R_i=ABS(R0_node(NPD_i,ii)-pmin_node(NPD_i))
                IF (R_i.LE.2*RLIMAX) THEN
C these atoms are moving to left node
                  NLL=NLL+1
C                  IF (NLL.GT.15*NBMAX1) THEN
C                    WRITE(*,*) NLL,' NLL is GT 15*NBMAX1,
C     &              increase Density in parameters.inc fd_nabor2'
C                    STOP
C                  ENDIF
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
                ENDIF
                ii=nclist_node(ii)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        node_cell_min_temp=node_cell_min
        node_cell_max_temp=node_cell_max
        node_cell_min_temp(NPD_i)=node_cell_max(NPD_i)-2
C Find the right side atoms that will be passed to right node
        DO i=node_cell_min_temp(1),node_cell_max_temp(1)
          DO j=node_cell_min_temp(2),node_cell_max_temp(2)
            DO k=node_cell_min_temp(3),node_cell_max_temp(3)
              cell=i+(j-1)*nncell(1)+(k-1)*nncell(2)
     &        *nncell(1)
              ii=nchead_node(cell)
              DO while (ii.NE.0)
                R_i=ABS(R0_node(NPD_i,ii)-pmax_node(NPD_i))
                IF (R_i.LE.2*RLIMAX) THEN
C these atoms are moving to RIGHT node
                  NRR=NRR+1
C                  IF (NRR.GT.15*NBMAX1) THEN
C                    WRITE(*,*) NRR,' NRR is GT 15*NBMAX1,
C     &              increase Density in parameters.inc fd_nabor2'
C                    STOP
C                  ENDIF
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
                ENDIF
                ii=nclist_node(ii)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        node_cell_min_temp=node_cell_min
        node_cell_max_temp=node_cell_max
        node_cell_max_temp(NPD_j)=node_cell_min(NPD_j)+2
C Find the DOWN side atoms that will be passed to DOWN node
        DO i=node_cell_min_temp(1),node_cell_max_temp(1)
          DO j=node_cell_min_temp(2),node_cell_max_temp(2)
            DO k=node_cell_min_temp(3),node_cell_max_temp(3)
              cell=i+(j-1)*nncell(1)+(k-1)*nncell(2)
     &        *nncell(1)
              ii=nchead_node(cell)
              DO while (ii.NE.0)
                R_j=ABS(R0_node(NPD_j,ii)-pmin_node(NPD_j))
                IF (R_j.LE.2*RLIMAX) THEN
C these atoms are moving to DOWN node
                  NDD=NDD+1
C                  IF (NDD.GT.15*NBMAX2) THEN
C                    WRITE(*,*) NDD,' NDD is GT 15*NBMAX2,
C     &              increase Density in parameters.inc fd_nabor2'
C                    STOP
C                  ENDIF
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
                ENDIF
                ii=nclist_node(ii)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        node_cell_min_temp=node_cell_min
        node_cell_max_temp=node_cell_max
        node_cell_min_temp(NPD_j)=node_cell_max(NPD_j)-2
C Find the UP side atoms that will be passed to UP node
        DO i=node_cell_min_temp(1),node_cell_max_temp(1)
          DO j=node_cell_min_temp(2),node_cell_max_temp(2)
            DO k=node_cell_min_temp(3),node_cell_max_temp(3)
              cell=i+(j-1)*nncell(1)+(k-1)*nncell(2)
     &        *nncell(1)
              ii=nchead_node(cell)
              DO while (ii.NE.0)
                R_j=ABS(R0_node(NPD_j,ii)-pmax_node(NPD_j))
                IF (R_j.LE.2*RLIMAX) THEN
C these atoms are moving to UP node
                  NUU=NUU+1
C                  IF (NUU.GT.15*NBMAX2) THEN
C                    WRITE(*,*) NUU,' NUU is GT 15*NBMAX2,
C     &              increase Density in parameters.inc fd_nabor2'
C                    STOP
C                  ENDIF
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
                ENDIF
                ii=nclist_node(ii)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        node_cell_min_temp=node_cell_min
        node_cell_max_temp=node_cell_max
        node_cell_max_temp(NPD_i)=node_cell_min(NPD_i)+2
        node_cell_max_temp(NPD_j)=node_cell_min(NPD_j)+2
C Find the LEFT DOWN side atoms that will be passed to LEFT DOWN side
        Do i=node_cell_min(1),node_cell_max_temp(1)
          DO j=node_cell_min(2),node_cell_max_temp(2)
            DO k=node_cell_min(3),node_cell_max_temp(3)
              cell=i+(j-1)*nncell(1)+(k-1)*nncell(2)
     &        *nncell(1)
              ii=nchead_node(cell)
              DO while (ii.NE.0)
                R_i=ABS(R0_node(NPD_i,ii)-pmin_node(NPD_i))
                R_j=ABS(R0_node(NPD_j,ii)-pmin_node(NPD_j))
                IF ((R_i.LE.2*RLIMAX).AND.(R_j.LE.2*RLIMAX)) THEN
C these atoms are moving to LEFT DOWN node
                  NLD=NLD+1
C                  IF (NLD.GT.15*NBMAX3) THEN
C                    WRITE(*,*) NLD,'NLD is GT 15*NBMAX3,
C     &              increase Density in parameters.inc fd_nabor2'
C                    STOP
C                  ENDIF
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
                ENDIF
                ii=nclist_node(ii)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        node_cell_min_temp=node_cell_min
        node_cell_max_temp=node_cell_max
        node_cell_max_temp(NPD_i)=node_cell_min(NPD_i)+2
        node_cell_min_temp(NPD_j)=node_cell_max(NPD_j)-2
C Find the LEFT UP side atoms that will be passed to LEFT UP side
        DO i=node_cell_min_temp(1),node_cell_max_temp(1)
          DO j=node_cell_min_temp(2),node_cell_max_temp(2)
            DO k=node_cell_min_temp(3),node_cell_max_temp(3)
              cell=i+(j-1)*nncell(1)+(k-1)*nncell(2)
     &        *nncell(1)
              ii=nchead_node(cell)
              DO while (ii.NE.0)
                R_i=ABS(R0_node(NPD_i,ii)-pmin_node(NPD_i))
                R_j=ABS(R0_node(NPD_j,ii)-pmax_node(NPD_j))
                IF ((R_i.LE.2*RLIMAX).AND.(R_j.LE.2*RLIMAX)) THEN
C these atoms are moving to LEFT UP node
                  NLU=NLU+1
C                  IF (NLU.GT.15*NBMAX3) THEN
C                    WRITE(*,*) NLU,'NLU is GT 15*NBMAX3,
C     &              increase Density in parameters.inc fd_nabor2'
C                    STOP
C                  ENDIF
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
                ENDIF
                ii=nclist_node(ii)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        node_cell_min_temp=node_cell_min
        node_cell_max_temp=node_cell_max
        node_cell_min_temp(NPD_i)=node_cell_max(NPD_i)-2
        node_cell_max_temp(NPD_j)=node_cell_min(NPD_j)+2
C Find the RIGHT DOWN side atoms that will be passed to RIGHT DOWN side
        DO i=node_cell_min_temp(1),node_cell_max_temp(1)
          DO j=node_cell_min_temp(2),node_cell_max_temp(2)
            DO k=node_cell_min_temp(3),node_cell_max_temp(3)
              cell=i+(j-1)*nncell(1)+(k-1)*nncell(2)
     &        *nncell(1)
              ii=nchead_node(cell)
              DO while (ii.NE.0)
                R_i=ABS(R0_node(NPD_i,ii)-pmax_node(NPD_i))
                R_j=ABS(R0_node(NPD_j,ii)-pmin_node(NPD_j))
                IF ((R_i.LE.2*RLIMAX).AND.(R_j.LE.2*RLIMAX)) THEN
C these atoms are moving to RIGHT DOWN node
                  NRD=NRD+1
C                  IF (NRD.GT.15*NBMAX3) THEN
C                    WRITE(*,*) NRD,'NRD is GT 15*NBMAX3,
C     &              increase Density in parameters.inc fd_nabor2'
C                    STOP
C                  ENDIF
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
                ENDIF
                ii=nclist_node(ii)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        node_cell_min_temp=node_cell_min
        node_cell_max_temp=node_cell_max
        node_cell_min_temp(NPD_i)=node_cell_max(NPD_i)-2
        node_cell_min_temp(NPD_j)=node_cell_max(NPD_j)-2
C Find the RIGHT UP side atoms that will be passed to RIGHT UP node
        DO i=node_cell_min_temp(1),node_cell_max_temp(1)
          DO j=node_cell_min_temp(2),node_cell_max_temp(2)
            DO k=node_cell_min_temp(3),node_cell_max_temp(3)
              cell=i+(j-1)*nncell(1)+(k-1)*nncell(2)
     &        *nncell(1)
              ii=nchead_node(cell)
              DO while (ii.NE.0)
                R_i=ABS(R0_node(NPD_i,ii)-pmax_node(NPD_i))
                R_j=ABS(R0_node(NPD_j,ii)-pmax_node(NPD_j))
                IF ((R_i.LE.2*RLIMAX).AND.(R_j.LE.2*RLIMAX)) THEN
C these atoms are moving to RIGHT UP node
                  NRU=NRU+1
C                  IF (NRU.GT.15*NBMAX3) THEN
C                    WRITE(*,*) NRU,'NRU is GT 15*NBMAX3,
C     &              increase Density in parameters.inc fd_nabor2'
C                    STOP
C                  ENDIF
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
                ENDIF
                ii=nclist_node(ii)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        N_atom_passLL(1)=NLL/15
        N_atom_passRR(1)=NRR/15
        N_atom_passLD(1)=NLD/15
        N_atom_passRU(1)=NRU/15
        N_atom_passRD(1)=NRD/15
        N_atom_passLU(1)=NLU/15
        N_atom_passDD(1)=NDD/15
        N_atom_passUU(1)=NUU/15
C        WRITE(*,*)'find_nabor2',mynode,NLL/15,NRR/15,NDD/15,NUU/15,'LL'
C        WRITE(*,*)'find_nabor2',mynode,NLD/15,NLU/15,NRD/15,NRU/15,'LD'
      ENDIF

      RETURN
      END


