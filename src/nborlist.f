C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Programer:
C     Travis Kemper
C     Department of Materials SCience and ENgineering
C     University of FLorida
C     traviskemper@ufl.edu
C
C     Version 1.0 10/06/09 T. W. Kemper
C     Version 2.0 10/10/09 T. W. Kemper
C     Version 3.0 01/13/11 T. W. Kemper
!        - same as was in caguts 
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  create nieghbor list based on cut offs 
!
      SUBROUTINE nborlist 
!
      USE STRUCTR
      USE POTS
      USE MDSTP
      USE PARAMS 
!
      IMPLICIT none
!
      INTEGER :: I,J,K,KI,KJ,L
      REAL*8 :: RSQ,RLIS,RR(3),RI(3)

      IF(LCHK.EQ.1) THEN
C
C Set up neighbor list
C
           K=0
           DO 302 I=1,NP
                NABORS(I)=K+1
                DO 299 L=1,3
                     RI(L)=R0(I,L)
299             CONTINUE
                KI=KTYPE(I)
C
C
                if(ki.gt.RTYPEs) go to 302
C
                DO 301 J=1,NP
C
                     IF(I.EQ.J) GO TO 301
C
                     KJ=KTYPE(J)
C
c cuts out all but C,H,Si, and Ge
C
                     if(kj.gt.RTYPES) go to 301
                     RLIS=RLIST(KI,KJ)
C
                     RSQ=0.0D0
                     DO 298 L=1,3
                          RR(L)=RI(L)-R0(J,L)
                          RR(L)=RR(L) -
     &                          CUBE(L)*ANINT(RR(L)/CUBE(L))
                          RSQ=RSQ+RR(L)*RR(L)
                          IF(RSQ.GT.RLIS) GO TO 301
298                  CONTINUE
C
405                  CONTINUE
                     K=K+1
                     LIST(K)=J
                     IVCT2B(K)=I
                     JVCT2B(K)=J
C
301             CONTINUE
302        CONTINUE
C
           NABORS(NP+1)=K+1
           KEND=K
           if(kend.gt.nlmax) then
                 write(*,*) 'kend exceeds nlmax'
                 write(*,*) 'kend,nlmax = ',kend,nlmax
                 write(*,*) 'increase nlmax and recompile' 
                 include 'close.inc' 
                 stop 
           endif
C           write(*,*) 'kend= ',kend
      ENDIF
!
      RETURN
      END SUBROUTINE nborlist
