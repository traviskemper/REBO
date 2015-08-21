       SUBROUTINE clink
C
C This subroutine divide the atoms into different cells
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      IMPLICIT REAL*8(A-H,O-Z)
c
      include 'common_files.inc'
      INCLUDE 'mpif.h'

C       INTEGER  ::icell,i,j,n
        REAL*8  ncell(3)
C       WRITE(*,*) 'start cagut_linklist'

       DO i=1,3
         ncell(i)=float(nncell(i))
       ENDDO

       ntcell=nncell(1)*nncell(2)*nncell(3)
C       WRITE(*,*) 'total # of cells=',ntcell
       IF (ntcell.GT.NCMAX) THEN
         WRITE(*,*) 'total # of cells=',ntcell
         WRITE(*,*) nncell(1),nncell(2),nncell(3)
         WRITE(*,*) 'increase NCMAX in parameter.inc'
         CALL write_error
       ENDIF
C       WRITE(*,*) nncell(1),nncell(2),nncell(3)
       DO icell=1,ntcell
         nchead_node(icell)=0
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

C       WRITE(*,*) 'normal execution clink complete'

       RETURN
       END 
