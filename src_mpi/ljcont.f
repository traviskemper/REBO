      subroutine ljcont
      USE dyn_array      
      USE CONTM
      USE POTS
      USE PARAMS
      USE STRUCTR
      USE MDSTP

      IMPLICIT none

      INTEGER :: I,J,kj
      REAL*8 :: CONV1,CONV2,ds,vsurf,dv
C
c use continuum boundary with LJ potential. See Goodman class
C notes for derivation of form.
C
C sigs and epss depend on atomic volume so need to be custom made
C for each new substrate
c
      do 10 i=1,RTYPES
           if(tau(i).lt.1.0d-05) go to 10
           do j=1,np
                kj = ktype(j)
                ds = r0(j,ndir) - surf
                s3 = (sigs(i,kj)/ds)**3
                vsurf = epss(i,kj)*s3*(s3*s3 - 3.0d0)
                dv = epss(i,kj)*s3*9.0d0/ds*(s3*s3 - 1.0d0)
                tote = tote + vsurf
                rnp(j,ndir) = rnp(j,ndir) + dv
           enddo
10     continue
       return
       end

