       SUBROUTINE check_center(NN,i)
       
C This subroutine check received atoms in read_coord
C if the atoms is belong to the correct
C if correct save as RO_node, R1_node, R2_node, R3_node
C if not then get rid of them
      USE dyn_array
      USE STRUCTR
      USE MPIvars

      IMPLICIT none

      INCLUDE 'mpif.h'
C Local variables
      INTEGER icheck
     & ,i,ii,iii,NN
C      icheck=0
      DO ii=1,np_send(i)
        IF ((R0(1,ii).GE.pmin_node(1)).AND.
     &      (R0(1,ii).LT.pmax_node(1))) THEN
          IF ((R0(2,ii).GE.pmin_node(2)).AND.
     &        (R0(2,ii).LT.pmax_node(2))) THEN
            IF ((R0(3,ii).GE.pmin_node(3)).AND.
     &          (R0(3,ii).LT.pmax_node(3))) THEN
              NN=NN+1
              NA(NN)=KA(ii)
              KTYPE_node(NN)=KTYPE(ii)
              itr_node(NN)=itr(ii)
              DO iii=1,3
                R0_node(iii,NN)=R0(iii,ii)
                R1_node(iii,NN)=R1(iii,ii)
                R2_node(iii,NN)=R2(iii,ii)
                R3_node(iii,NN)=R3(iii,ii)
              ENDDO
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      
C      IF (icheck.EQ.1) THEN
C        WRITE(*,*) 'there are missing atoms'
C        CALL write_error
C      ENDIF
C      WRITE(*,*) 'finish chcek_center(i)'
      RETURN
      END

