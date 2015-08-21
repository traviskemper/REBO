
      SUBROUTINE open_file

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
        open(10,file='coord.d',status='old')
        open(11,file='out.d',status='unknown')
        open(12,file='input.d',status='old')
        OPEN(13,FILE='out.warn',STATUS='UNKNOWN')
        open(21 ,file='out.dat',status='unknown')
        IF(PRNLD)open(51,file='out.load',status='unknown')
      ENDIF
C     
c$$$      nfile = 11
c$$$      filename = 'coord_node_'
c$$$      if (nprocs+1.lt.10) then
c$$$          nfile=12
c$$$          filename='coord_node_0'
c$$$          write (filename(nfile+1:),'(i1)') mynode
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
c$$$      OPEN(mynode+20,file=filename,STATUS='UNKNOWN')
c
c data files for REBO potential
C
      IF (mynode.EQ.0) THEN
        open(14,file='../Spline/inter3d_iv_new.d',
     &  status='old')
        open(15,file='../Spline/inter3d_iv_CHF.d',
     &  status='old')
        open(16,file='../Spline/inter3dtors.d',
     &  status='old')
        open(17,file='../Spline/inter3d_h.d',
     &  status='old')
        open(18,file='../Spline/inter3d_ch.d',
     &  status='old')
      ENDIF
c
c optional input/output files
c
C      open(55,file='max_ke.d',status='unknown')
C      open(85,file='pair_energy.d',status='unknown')
C      open(50,file='overwrite.d',status='unknown')
C
C output files for tight binding potential
c
C       open(85,file='eigenvectors.d',status='unknown')
C       open(86,file='eigenenergies.d',status='unknown')
C       open(87,file='dos.d',status='unknown')
C       open(88,file='ldos.d',status='unknown')

C       WRITE(*,*) 'finish open_files'
       RETURN
       END
