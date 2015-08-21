       SUBROUTINE cell_range
       
C This subroutine find the node_cell_min(NPD_i) and node_cell_max(NPD_i)
C on each node

      USE MPIvars
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'mpif.h'
      include 'common_files.inc'
      
C Local variable
      INTEGER icell_min,icell_max,cell_i(3)
      REAL*8 ncell(3)
C Find node_cell_min AND node_cell_max
       icell_min=nncell(1)*nncell(2)*nncell(3)
       icell_max=1
       DO i=1,3
         ncell(i)=float(nncell(i))
       ENDDO
       DO i=1,NP_node
         icell=1+INT(((r0_node(1,i)/(CUBE(1)+1.0d-12))+0.5)*ncell(1))
     &          +INT(((r0_node(2,i)/(CUBE(2)+1.0d-12))+0.5)*ncell(2))
     &          *nncell(1)
     &          +INT(((r0_node(3,i)/(CUBE(3)+1.0d-12))+0.5)*ncell(3))
     &          *nncell(2)*nncell(1)
         IF (icell.GT.icell_max) THEN
           icell_max=icell
         ENDIF
         IF (icell.LT.icell_min) THEN
           icell_min=icell
         ENDIF
       ENDDO
       IF (NNDP.EQ.1) THEN
C Find node_cell_max(NPD_i)
         cell_i(1)=MOD(icell_max,nncell(1))
         IF (cell_i(1).EQ.0) THEN
           cell_i(1)=nncell(1)
         ENDIF
         cell_i(3)=INT(icell_max/(nncell(1)*nncell(2)+1.0d-12))+1
         cell_i(2)=INT((icell_max-(cell_i(3)-1)*nncell(1)*nncell(2))/
     &   (nncell(1)+1.0d-12))+1
         node_cell_max(NPD_i)=cell_i(NPD_i)
         cell_i(1)=MOD(icell_min,nncell(1))
         IF (cell_i(1).EQ.0) THEN
           cell_i(1)=nncell(1)
         ENDIF
         cell_i(3)=INT(icell_min/(nncell(1)*nncell(2)+1.0d-12))+1
         cell_i(2)=INT((icell_min-(cell_i(3)-1)*nncell(1)*nncell(2))/
     &   (nncell(1)+1.0d-12))+1
         node_cell_min(NPD_i)=cell_i(NPD_i)
       ELSE
         cell_i(1)=MOD(icell_max,nncell(1))
         IF (cell_i(1).EQ.0) THEN
           cell_i(1)=nncell(1)
         ENDIF
         cell_i(3)=INT(icell_max/(nncell(1)*nncell(2)+1.0d-12))+1
         cell_i(2)=INT((icell_max-(cell_i(3)-1)*nncell(1)*nncell(2))/
     &   (nncell(1)+1.0d-12))+1
         node_cell_max(NPD_i)=cell_i(NPD_i)
         node_cell_max(NPD_j)=cell_i(NPD_j)
         cell_i(1)=MOD(icell_min,nncell(1))
         IF (cell_i(1).EQ.0) THEN
           cell_i(1)=nncell(1)
         ENDIF
         cell_i(3)=INT(icell_min/(nncell(1)*nncell(2)+1.0d-12))+1
         cell_i(2)=INT((icell_min-(cell_i(3)-1)*nncell(1)*nncell(2))/
     &   (nncell(1)+1.0d-12))+1
         node_cell_min(NPD_i)=cell_i(NPD_i)
         node_cell_min(NPD_j)=cell_i(NPD_j)
       ENDIF
      WRITE(*,*) mynode,'node_cell x',node_cell_min(1),node_cell_max(1)
      WRITE(*,*) mynode,'node_cell y',node_cell_min(2),node_cell_max(2)
      WRITE(*,*) mynode,'node_cell z',node_cell_min(3),node_cell_max(3)

       RETURN
       END
     

