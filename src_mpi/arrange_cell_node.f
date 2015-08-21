      SUBROUTINE arrange_cell
C
C Right now all the nodes know all the informations of all atoms
C Here atoms will be rearranged in 2-D by spacial relationship
C This arrangement will be applied only for pi bond calculation
C First sort all the atoms into cells
C Second assign the cells to nodes
C third assign the atoms in the same cell find in second step to nodes
C Forth find the atoms in the first and second neighbor cells
C Fifth add those atoms found in forth step to the final of the atoms found in the third step
C
      USE MPIvars
      USE POTS
      USE PARAMS
      USE TABS
      USE STRUCTR 

      IMPLICIT none

      INCLUDE 'mpif.h'

C local variables
      INTEGER MNDP,NNNDP(2),node_i,node_j
     &,NN,I,J,KI,KJ
      REAL*8 RLIcell(3)
C MPI variables

C Check 1-D or 2-D parallel read from input.d
      NN=0
      DO i=1,3
        IF (NDP(i).NE.1) THEN
          NN=NN+1
        ENDIF
      ENDDO
      IF (NN.EQ.1) THEN
        NNDP=1
      ENDIF
      IF (NN.EQ.2) THEN
        NNDP=2
      ENDIF
      IF (NN.EQ.3) THEN
        WRITE(*,*) 'This version only has 1-D and 2-D parallel'
        CALL write_error
      ENDIF
      MNDP=NDP(1)*NDP(2)*NDP(3)
      IF (MNDP.NE.nprocs) THEN
        WRITE(*,*) '# of CPU set in input.d is not the same as #s
     &  set in the script'
        CALL write_error
      ENDIF
      IF (MNDP.EQ.1) THEN
        WRITE(*,*) 'It is a parallel code, 
     &  mutiple-CPU usage is necessary'
        CALL write_error
      ENDIF
      IF (NNDP.EQ.1) THEN
        IF (MOD(nprocs,2).NE.0) THEN
          WRITE(*,*) '1-D parallel, # of CPUs should be devide by 2'
          CALL write_error
        ENDIF
      ENDIF
      IF (NNDP.EQ.2) THEN
        IF (MOD(nprocs,4).NE.0) THEN
          WRITE(*,*) '2-D parallel, # of CPUs should be devide by 4'
          CALL write_error
        ENDIF
      ENDIF
C Find the maximun cut-off to make the cell
       RLIMAX=0.0
       DO KI=1,3
         DO KJ=1,3
C           WRITE(*,*) RLIST(KI,KJ)
           IF (RLIMAX.LT.SQRT(RLIST(KI,KJ))) THEN
             RLIMAX=SQRT(RLIST(KI,KJ))
           ENDIF
         ENDDO
       ENDDO
C       WRITE(*,*) 'arrange_cell-RLIMAX=',RLIMAX
C Find the maximum cut-off for buffer region
       RLJMAX=0.0
       DO KI=1,KTMAX
         DO KJ=1,KTMAX
           IF (RLJMAX.LT.SQRT(RSLJ(KI,KJ))) THEN
             RLJMAX=SQRT(RSLJ(KI,KJ))
           ENDIF
         ENDDO
       ENDDO
C Determine how to devide the system
      DO i=1,3
        pmin_node(i)=-CUBE(i)/2
        pmax_node(i)=CUBE(i)/2+1.0d-12
        node_cell_min(i)=1
        node_cell_max(i)=INT(CUBE(i)/RLIMAX)
        nncell(i)=INT(CUBE(i)/RLIMAX)
C        WRITE(*,*) i,CUBE(i),RLIMAX
C        WRITE(*,*) nncell(i)
        IF (nncell(i).LT.3) THEN
          WRITE(*,*) 'System size is too small',nncell(i)
          WRITE(*,*) 'Minimum size=',3*RLIMAX
          CALL write_error
        ENDIF
        nncell_lj(i)=INT(CUBE(i)/RLJMAX)
        IF (nncell_lj(i).LT.3) THEN
          WRITE(*,*) 'LJ System size is too small',nncell_lj(i)
          WRITE(*,*) 'Minimum size=',3*RLJMAX
          CALL write_error
        ENDIF
      ENDDO
C      WRITE(*,*) 'nncell',nncell(1),nncell(2),nncell(3)
      IF (NNDP.EQ.1) THEN
        DO i=1,3
          IF (NDP(i).NE.1) THEN
            IF ((CUBE(i)/NDP(i)).LE.RLJMAX) THEN
                WRITE(*,*) 'use too many CPUs, center region smaller
     &          than buffer region'
              CALL write_error
            ENDIF
            NPD_i=i
            RLIcell(i)=CUBE(i)/nncell(i)
            pmin_node(i)=(-CUBE(i)/2)+(mynode*CUBE(i)/NDP(i))
            pmax_node(i)=(-CUBE(i)/2)+((mynode+1)*CUBE(i)/NDP(i))
            IF (mynode.EQ.nprocs) THEN
              pmax_node(i)=CUBE(i)/2+1.0d-12
            ENDIF
            node_cell_min(i)=INT((pmin_node(i)-1.0d-12+CUBE(i)/2)/
     &                       RLIcell(i))+1
            node_cell_max(i)=INT((pmax_node(i)-1.0d-12+CUBE(i)/2)/
     &                       RLIcell(i))+1
            IF (node_cell_max(i).GT.nncell(i)) THEN
              node_cell_max(i)=nncell(i)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
C      WRITE(*,*) 'NPD_i=',NPD_i
      IF (NNDP.EQ.2) THEN
        DO i=1,2
          DO j=i+1,3
            IF ((NDP(i).NE.1).AND.(NDP(j).NE.1)) THEN
              IF ((CUBE(i)/NDP(i)).LE.RLJMAX) THEN
                WRITE(*,*) 'use too many CPUs, center region smaller
     &          than buffer region'
                CALL write_error
              ENDIF
              IF ((CUBE(j)/NDP(j)).LE.RLJMAX) THEN
                WRITE(*,*) 'use too many CPUs, center region smaller
     &          than buffer region'
                CALL write_error
              ENDIF
              NPD_i=i
              NPD_j=j
              NPD_k=INT(6-i-j)
              RLIcell(i)=CUBE(i)/nncell(i)
              RLIcell(j)=CUBE(j)/nncell(j)
              node_i=MOD((mynode+1),NDP(i))
              IF (node_i.EQ.0) THEN
                node_i=NDP(i)
              ENDIF
              node_j=INT((mynode+1)/(NDP(i)+1.0d-12))+1
              pmin_node(i)=(-CUBE(i)/2)+((node_i-1)*CUBE(i)/NDP(i))
              pmax_node(i)=(-CUBE(i)/2)+((node_i)*CUBE(i)/NDP(i))
              pmin_node(j)=(-CUBE(j)/2)+((node_j-1)*CUBE(j)/NDP(j))
              pmax_node(j)=(-CUBE(j)/2)+((node_j)*CUBE(j)/NDP(j))
              IF (node_i.EQ.NDP(i)) THEN
                pmax_node(i)=CUBE(i)/2+1.0d-12
              ENDIF
              IF (node_j.EQ.NDP(j)) THEN
                pmax_node(j)=CUBE(j)/2+1.0d-12
              ENDIF
              node_cell_min(i)=INT((pmin_node(i)-1.0d-12+CUBE(i)/2)/
     &                         RLIcell(i))+1
              node_cell_max(i)=INT((pmax_node(i)-1.0d-12+CUBE(i)/2)/
     &                         RLIcell(i))+1
              IF (node_cell_max(i).GT.nncell(i)) THEN
                node_cell_max(i)=nncell(i)
              ENDIF
              node_cell_min(j)=INT((pmin_node(j)-1.0d-12+CUBE(j)/2)/
     &                         RLIcell(j))+1
              node_cell_max(j)=INT((pmax_node(j)-1.0d-12+CUBE(j)/2)/
     &                         RLIcell(j))+1
              IF (node_cell_max(j).GT.nncell(j)) THEN
                node_cell_max(j)=nncell(j)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDIF

c$$$      WRITE(*,*) mynode,'devide in x',pmin_node(1),pmax_node(1)
c$$$      WRITE(*,*) mynode,'devide in y',pmin_node(2),pmax_node(2)
c$$$      WRITE(*,*) mynode,'devide in z',pmin_node(3),pmax_node(3)
c$$$      WRITE(*,*) mynode,'node_cell x',node_cell_min(1),node_cell_max(1)
c$$$      WRITE(*,*) mynode,'node_cell y',node_cell_min(2),node_cell_max(2)
c$$$      WRITE(*,*) mynode,'node_cell z',node_cell_min(3),node_cell_max(3)
c$$$      WRITE(*,*) 'finish devide system by processors'


      RETURN
      END
