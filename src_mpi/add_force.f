!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Programer:
!     Travis Kemper
!     Department of Materials SCience and ENgineering
!     University of FLorida
!     traviskemper@ufl.edu
!
!     Version 1.0 04/04/11 T. W. Kemper
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Adds a force to surface atoms for compressition study
!
      subroutine ADD_FORCE

      USE dyn_array
      USE MPIvars
      USE STRUCTR
      USE SPECIF

      IMPLICIT none

      INTEGER :: i,j
C
      do i = 1,NP_node
           if(itr_node(i).eq.3) then
                 RNP(SDIR,I)=RNP(SDIR,I) + ADFB
           endif
           if(itr_node(i).eq.4) then
                 RNP(SDIR,I)=RNP(SDIR,I) + ADFT
           endif
      enddo
c
c add velocity to carbon tip for first step
c      if(ttime.le.delta) then
c           do i=1,np
c                if((ktype(i).le.2).and.(itr(i).ne.2)) then
c                     r1(i,2) = -0.0005d0
c                endif
c           enddo
c      endif
      return
      END SUBROUTINE ADD_FORCE

