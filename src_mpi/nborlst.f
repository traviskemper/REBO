C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Programer:
C     Travis Kemper
C     Department of Materials SCience and ENgineering
C     University of FLorida
C     traviskemper@ufl.edu
C
C     Version 1.0 10/06/09 T. W. Kemper
C     Version 2.0 10/10/09 T. W. Kemper
C       - add number of non-hydrogen neighbors
C          -1 - hydrogen
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  create nieghbor list based on rebo cut offs 

      SUBROUTINE nborlst 

      USE MPIvars
      USE STRUCTR
      USE POTS

      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
      INTEGER ::  I,J,N
      REAL*8 :: RIJSQ,DR,
     &    RCUTSQ,RHCUTSQ,HHCUTSQ


      IF ( mynode.EQ.0 ) THEN
C Set squares of cut offs
      RCUTSQ = RCUT*RCUT
      RHCUTSQ = RHCUT*RHCUT
      HHCUTSQ = HHCUT*HHCUT
 
C Find the nearest nieghbors for each atom 
      DO I=1,NP
        NCTi(I) = 0 
        CNNi(I) = 0
      ENDDO

      IF (mynode.EQ.0) WRITE(*,*)'Atoms=',NP
      IF (mynode.EQ.0) WRITE(*,*)'box=',CUBE(:)
      DO I=1,NP-1
        DO J=I+1,NP
          RIJSQ = 0d0
          DO N=1,3
            DR = R0i(N,J) - R0i(N,I)  
            DR = DR - 
     &         CUBE(N)*ANINT(DR/CUBE(N))
            RIJSQ =  RIJSQ + DR*DR
          ENDDO
!         IF(I.EQ.34)  WRITE(46,*) I,J,RIJSQ
          IF ( RIJSQ.LT.RCUTSQ)  THEN
            IF ( KTYPEi(J).eq.ihyd .and. KTYPEi(I).eq.ihyd
     &           .and. RIJSQ.lt.HHCUTSQ ) THEN
              NCTi(I) = NCTi(I) + 1
              NCTi(J) = NCTi(J) + 1
              NLISTi(I,NCTi(I)) = J  
              NLISTi(J,NCTi(J)) = I
            ELSEIF ( KTYPEi(J).eq.ihyd.and.KTYPEi(I).ne.ihyd
     &             .and. RIJSQ.lt.RHCUTSQ ) THEN
              NCTi(I) = NCTi(I) + 1
              NCTi(J) = NCTi(J) + 1
              NLISTi(I,NCTi(I)) = J  
              NLISTi(J,NCTi(J)) = I
              CNNi(J) = CNNi(J) + 1
            ELSEIF ( KTYPEi(I).eq.ihyd .and. KTYPEi(J).ne.ihyd
     &             .and. RIJSQ.lt.RHCUTSQ ) THEN
              NCTi(I) = NCTi(I) + 1
              NCTi(J) = NCTi(J) + 1
              NLISTi(I,NCTi(I)) = J  
              NLISTi(J,NCTi(J)) = I
              CNNi(I) = CNNi(I) + 1
            ELSEIF ( KTYPEi(J).ne.ihyd.and.KTYPEi(I).ne.ihyd )  THEN
              NCTi(I) = NCTi(I) + 1
              NCTi(J) = NCTi(J) + 1
              NLISTi(I,NCTi(I)) = J  
              NLISTi(J,NCTi(J)) = I
              CNNi(I) = CNNi(I) + 1
              CNNi(J) = CNNi(J) + 1
            ENDIF
          ENDIF
          IF ( NCTi(I).eq.NMX.or.NCTi(J).eq.NMX ) THEN
            WRITE(*,*) 'The number of neighbors of',I,J
            WRITE(*,*) ' is ',NCTi(I)
            WRITE(*,*) 'equal to the max defined in modual_analysis.f'
            WRITE(*,*) R0i(:,I)
            WRITE(*,*) R0i(:,J)
            DO N = 1,NCTi(I)
              WRITE(*,*) N,NLISTi(I,N),KTYPEi(NLISTi(I,N))
            ENDDO
            STOP
          ENDIF
        ENDDO
      ENDDO
      ENDIF

      END SUBROUTINE nborlst
