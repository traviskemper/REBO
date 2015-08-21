      SUBROUTINE LJCSET

      USE CONTM
      USE POTS
      USE PARAMS

      IMPLICIT none

      INTEGER :: I,J
      REAL*8 :: CONV1,CONV2

C     ALLOCATE(
C    & TAU(RTYPES)
C     &,EPSS(RTYPES,RTYPES)
C     &,SIGS(RTYPES,RTYPES)
C     &)


C set up for continuum
c
C  Atomic densities must be custom made for each substrate
C
C tau: No. atoms per unit volume (in A**3) for each atom type
C      It is 1/t from Goodwin's notes
C
C surf: location of surface
C
C ndir: direction of boundary
C
C note that epss(atom in continuum, atom above surface)
C
      ndir = 2
      surf = 0.00d0
c
      do i=1,NTYPES
           tau(i) = 0.0d0
      enddo

c polyethylene
      tau(1) = 64.0d0/1385.1466d0
      tau(2) = 128.0d0/1385.1466d0
c diamond (111)
c      tau(1) = 96.0d0/(10.085d0*8.734d0*6.176d0)

      conv1 = (0.40d0)**(1.0d0/6.0d0)
      conv2 = 10.0d0*pi/9.0d0/8.0d0
c last divide above accounts for *4 and 1/2 in class notes
      do i=1,NTYPES
           do j=1,NTYPES
                sigs(i,j) = conv1*sqrt(sig(i,j))
                epss(i,j) = conv2*eps(i,j)*tau(i)*(sigs(i,j)**3)
           enddo
      enddo
      return
      end

