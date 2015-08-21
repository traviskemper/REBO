       SUBROUTINE ljgut_nabor_cell
       
C This subroutine first devide system into cells by maximum L-J cut-off
C Then find the neighbor cells for every cell
C save neighbor cell list in an array "ljnbc(27)"
C the 14th number in that array is the self cell

      USE dyn_array
      USE MPIvars

      IMPLICIT none

C local variables
       REAL*8  ncell(3)
       INTEGER  :: icell,cellx,celly,cellz
     & ,x,y,z,n,i,ntcell
       INTEGER  newcx,newcy,newcz

C       DO i=1,3
C         nncell_lj(i)=INT(CUBE(i)/RLJMAX)
C         IF (nncell_lj(i).LT.3) THEN
C           WRITE(*,*) 'System size is too small'
C           WRITE(*,*) 'Size of the system should be',3*RLJMAX
C           CALL write_error
C         ENDIF
C       ENDDO
C
       ntcell=nncell_lj(1)*nncell_lj(2)*nncell_lj(3)
CC       WRITE(*,*) 'total # of cells=',ntcell
C       IF (ntcell.GT.NPMAX) THEN
C         WRITE(6,*) ntcell,'increase NbufMAX in parameters.inc'
C         CALL write_error
C       ENDIF
C       WRITE(*,*) nncell_lj(1),nncell_lj(2),nncell_lj(3)

C find neighbor cells

       DO x=1,3
         ncell(x)=float(nncell_lj(x))
       ENDDO

       DO i=1,ntcell
         icell=i
         cellx=MOD(icell,nncell_lj(1))
         IF (cellx.EQ.0) THEN
           cellx=nncell_lj(1)
         ENDIF
         cellz=INT(icell/(nncell_lj(1)*nncell_lj(2)+1.0d-12))+1
         celly=INT((icell-(cellz-1)*nncell_lj(1)*nncell_lj(2))/
     &  (nncell_lj(1)+1.0d-12))+1

         n=0

         DO x=-1,1
           DO y=-1,1
             DO z=-1,1
               newcx=cellx+x
               IF (newcx.LT.1) THEN
                 newcx=newcx+nncell_lj(1)
               ENDIF
               IF (newcx.GT.nncell_lj(1)) THEN
                 newcx=newcx-nncell_lj(1)
               ENDIF
               newcy=celly+y
               IF (newcy.LT.1) THEN
                 newcy=newcy+nncell_lj(2)
               ENDIF
               IF (newcy.GT.nncell_lj(2)) THEN
                 newcy=newcy-nncell_lj(2)
               ENDIF
               newcz=cellz+z
               IF (newcz.LT.1) THEN
                 newcz=newcz+nncell_lj(3)
               ENDIF
               IF (newcz.GT.nncell_lj(3)) THEN
                 newcz=newcz-nncell_lj(3)
               ENDIF
               n=n+1
               ljnbc((i-1)*27+n)=((newcz-1)*nncell_lj(1)*nncell_lj(2))+
     &         (newcy-1)*nncell_lj(1)+newcx
C             WRITE(*,*)ljnbc(n)
             ENDDO
           ENDDO
         ENDDO

       ENDDO

C       WRITE(*,*) 'normal execution ljnabor complete:CALL write_errorping'

       RETURN
       END
