
      SUBROUTINE close_file

      USE MPIvars
      USE SPECIF 

      IMPLICIT none
c
      INCLUDE 'mpif.h'

C Local variables
      CHARACTER*20 filename
      INTEGER  nfile
c standard input/output files
      IF (mynode.EQ.0) THEN
       CLOSE(1)
       close(13)
       close(10)
       close(9)
       close(21)
       IF(PRNLD)CLOSE(51)
      ENDIF
     
c$$$      nfile = 11
c$$$      filename = 'coord_node_'
c$$$      if (nprocs+1.lt.10) then
c$$$          nfile=12
c$$$          filename='coord_node_0'
c$$$          nfile = nfile+1
c$$$C      elseif (nprocs+1.lt.100) then
c$$$C          filename(nfile+1:) ='0'
c$$$C          nfile = nfile+1
c$$$C          write (filename(nfile+1:),'(i2)') in_beam
c$$$C          nfile = nfile+2
c$$$      else
c$$$          write (filename(nfile+1:),'(i2)') mynode
c$$$          nfile = nfile+2
c$$$      endif
c$$$      filename(nfile+1:) = '.d'
c$$$      CLOSE(mynode+20)
c
c data files for REBO potential
C
      IF (mynode.EQ.0) THEN
        CLOSE(14)
        CLOSE(15)
        CLOSE(16)
        CLOSE(17)
        CLOSE(18)
      ENDIF

      RETURN
      END
