! Subroutine to open input files
!
      SUBROUTINE OPENIN
!
      USE SPECIF
!
c standard input/output files 
      open(10,file='coord.d',status='old') 
      OPEN(12,file='input.d',status='old')
      OPEN(13,FILE='out.warn',STATUS='UNKNOWN')
c 
c data files for REBO potential 
C 
      open(14,file='../Spline/inter3d_iv_new.d',
     &status='old')
      open(15,file='../Spline/inter2d_iv.d',
     &status='old')
      open(16,file='../Spline/inter3dtors.d',
     &status='old')
      open(17,file='../Spline/inter3d_h.d',
     &status='old')
      open(18,file='../Spline/inter3d_ch.d',
     &status='old')
      open(19,file='../Spline/inter3d_iv_CHF.d',
     &  status='old')


!

      RETURN
      END SUBROUTINE OPENIN

! Subroutine to open output files
!
      SUBROUTINE OPENOUT
!
      USE SPECIF
!
c 
c optional input/output files
c 
      IF(PXMOL) open(1,file='out.xmol',status='unknown') 
      IF(WSTRI ) OPEN(UNIT=13,FILE='initial.d',STATUS='UNKNOWN')
      OPEN(UNIT=11,FILE='out.coord',STATUS='unknown')
      OPEN(UNIT=21,FILE='out.data',STATUS='unknown')
      WRITE(21,210)
 210  FORMAT( '#  Step ; Time ; Total energy ; Kinetic energy '
     & ,'; Potential energy , Temperature ; Time step ') 
C
C output files for tight binding potential 
c 
      RETURN
      END SUBROUTINE OPENOUT

