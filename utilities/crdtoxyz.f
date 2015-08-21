c This program is to write the coordinates from the new code into old code, more
c specifically, into xmol format.    
C  v1.1 prints max and min values of x,y,z
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Programer:
C     Travis Kemper
C     Department of Materials SCience and ENgineering
C     University of FLorida
C     traviskemper@ufl.edu
C
C     Version 1.0 7/15/09 T. W. Kemper
C     Version 1.0 7/27/09 T. W. Kemper
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Input files:
C       [file].d    
C     Output files:
C       [file].xyz
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      PROGRAM crdtoxyz

      IMPLICIT none

      INTEGER :: AMX,J,N
      PARAMETER ( AMX = 150000 ) ! Maxiumum number of atoms insingle file
      INTEGER NA,PARA,NUM,ITYPE(AMX)
      REAL*8 :: R(AMX,3),RMAX(3),RMIN(3)

C Set max and min values for R 
      RMAX(:)=0
      RMIN(:)=100

C  Read in values from coord.d
      READ(5,*)
      READ(5,*) NA
      READ(5,*)
      READ(5,*)
 
      WRITE(6,*)NA
      Do J=1,NA
        READ(5,*)  NUM,ITYPE(J), R(J,1),R(J,2),R(J,3),PARA
        DO N=1,3
          IF ( R(J,N) .GT. RMAX(N) ) RMAX(N)=R(J,N) 
          IF ( R(J,N) .LT. RMIN(N) ) RMIN(N)=R(J,N) 
        ENDDO
      ENDDO
      WRITE(6,162) 'Max:',RMAX(1),RMAX(2),RMAX(3),
     &             'Min',RMIN(1),RMIN(2),RMIN(3)

      DO J=1,NA
        WRITE(6,161) ITYPE(J), R(J,1),R(J,2),R(J,3)
      ENDDO

  161 FORMAT(I7,3e16.8)     
  162 FORMAT(A5,3f12.8,A5,3f12.8)     

      STOP  
      END    


 
 

