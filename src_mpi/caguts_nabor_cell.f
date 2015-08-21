       SUBROUTINE cagut_nabor_cell
       
C This subroutine first devide system into cells by maximum REBO cut-off
C Then find the neighbor cells for every cell
C save neighbor cell list in an array "cneigh(27)"
C the 14th number in that array is the self cell

      USE dyn_array
      USE MPIvars

      IMPLICIT none
c

C local variables
       REAL*8  ncell(3)
       INTEGER  icell,cellx,celly,cellz,ntcell
     & ,newcx,newcy,newcz,i
     & ,x,y,z,n

       ntcell=nncell(1)*nncell(2)*nncell(3)
C       WRITE(*,*) 'total # of cells=',ntcell
C       IF (ntcell.GT.NCMAX) THEN
C         WRITE(*,*) 'total # of cells=',ntcell
C         WRITE(6,*) 'increase NCMAX in parameters.inc'
C         CALL write_error
C       ENDIF
C       WRITE(*,*) nncell(1),nncell(2),nncell(3)

C find neighbor cells

C This does not seem to be used and x is undifined 
C -travisk
c$$$       DO i=1,3
c$$$         ncell(x)=float(nncell(i))
c$$$       ENDDO

       DO i=1,ntcell
         icell=i
         cellx=MOD(icell,nncell(1))
         IF (cellx.EQ.0) THEN
           cellx=nncell(1)
         ENDIF
         cellz=INT(icell/(nncell(1)*nncell(2)+1.0d-12))+1
         celly=INT((icell-(cellz-1)*nncell(1)*nncell(2))/
     &  (nncell(1)+1.0d-12))+1

         n=0

         DO x=-1,1
           DO y=-1,1
             DO z=-1,1
               newcx=cellx+x
               IF (newcx.LT.1) THEN
                 newcx=newcx+nncell(1)
               ENDIF
               IF (newcx.GT.nncell(1)) THEN
                 newcx=newcx-nncell(1)
               ENDIF
               newcy=celly+y
               IF (newcy.LT.1) THEN
                 newcy=newcy+nncell(2)
               ENDIF
               IF (newcy.GT.nncell(2)) THEN
                 newcy=newcy-nncell(2)
               ENDIF
               newcz=cellz+z
               IF (newcz.LT.1) THEN
                 newcz=newcz+nncell(3)
               ENDIF
               IF (newcz.GT.nncell(3)) THEN
                 newcz=newcz-nncell(3)
               ENDIF
               n=n+1
               icneigh((i-1)*27+n)=((newcz-1)*nncell(1)*nncell(2))+
     &         (newcy-1)*nncell(1)+newcx
C             WRITE(*,*)cneigh(n)
             ENDDO
           ENDDO
         ENDDO

       ENDDO

C       WRITE(*,*) 'normal execution cnabor complete:CALL write_errorping'

       RETURN
       END
