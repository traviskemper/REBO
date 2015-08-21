       SUBROUTINE ljlink
C
C This subroutine divide the atoms into different cells
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      USE dyn_array
      USE MPIvars
      USE STRUCTR

      IMPLICIT none
c
      INTEGER  icell,i,j,n,NTCELL
      REAL*8  ncell(3)

       DO i=1,3
C         nncell_lj(i)=INT(CUBE(i)/RLJMAX)
         ncell(i)=float(nncell_lj(i))
       ENDDO

       ntcell=nncell_lj(1)*nncell_lj(2)*nncell_lj(3)
C       WRITE(*,*) 'total # of cells=',ntcell
C       IF (ntcell.GT.NPMAX) THEN
C         WRITE(*,*)'total # of cells=',ntcell
C         WRITE(*,*) nncell_lj(1),nncell_lj(2),nncell_lj(3)
C         WRITE(*,*) 'increase NPMAX in parameters.inc'
C         CALL write_error
C       ENDIF
C       WRITE(*,*) 'ljcell',nncell_lj(1),nncell_lj(2),nncell_lj(3)
       DO icell=1,ntcell
         ljhead_node(icell)=0
       ENDDO

       DO i=1,NPtot4_node
         icell=1+INT(((r0_node(1,i)/(CUBE(1)+1.0d-12))+0.5)*ncell(1))
     &          +INT(((r0_node(2,i)/(CUBE(2)+1.0d-12))+0.5)*ncell(2))
     &          *nncell_lj(1)
     &          +INT(((r0_node(3,i)/(CUBE(3)+1.0d-12))+0.5)*ncell(3))
     &          *nncell_lj(2)*nncell_lj(1)
         ljlist_node(i)=ljhead_node(icell)
         ljhead_node(icell)=i
         n_icell_lj(i)=icell
       ENDDO
C       WRITE(*,*) 'normal execution ljlink complete:CALL write_errorping'

       RETURN
       END 
