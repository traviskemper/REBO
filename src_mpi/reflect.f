      SUBROUTINE REFLECT
C
      IMPLICIT REAL*8(A-H,O-Z)
c
      include 'common_files.inc'
C
C use reflecting boundary in y
C
      do i=1,np
           if((r0(i,2).gt.38.50d0).and.(r1(i,2).gt.0.0d0))
     &      r1(i,2) = -r1(i,2)
      enddo
      return
      end

