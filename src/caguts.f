C***WARNING***WARNING***WARNING***WARNING***WARNING***WARNING***
C
C This version of the hydrocarbon potential is currently
C (8/3/95) unpublished. There may still be some small changes
C forthcoming, particularly with respect to interstitial
C hydrogen in diamond. For updates, e-mail D. Brenner at
C dwb@ripley.mte.ncsu.edu.
C
C input coordinate file coord.d is written over; remember to
C keep a backup copy!!!!!
C
C***WARNING***WARNING***WARNING***WARNING***WARNING***WARNING***

      SUBROUTINE CAGUTS
C
C CALCULATE TWO-BODY FORCES AND NEIGHBOR LIST for hydrocarbons
C
      USE SPECIF
      USE MDSTP
      USE STRUCTR
      USE POTS
      USE PARAMS
      USE TABS

      IMPLICIT none

      INTEGER :: I,J,K,L,M,N,KI,KJ,IT,JBEGIN,JEND,NBJ
      REAL*8 :: RPP(nlmax,3),RR(3),RI(3),RLIS,RSQ,RC,RT,VV,RP

c
      DO 320 K=1,KEND
C
           I=IVCT2B(K)
           J=JVCT2B(K)

           KI=KTYPE(I)
           KJ=KTYPE(J)
C
           LCHECK(K)=0
           RSQ=0.0D0
           DO L=1,3
                RR(L)=R0(I,L)-R0(J,L)
                RR(L)=RR(L) - CUBE(L)*ANINT(RR(L)/CUBE(L))
                RSQ=RSQ+RR(L)*RR(L)
                COR(K,L)=RR(L)
           ENDDO
c
           IF(RSQ.GT.RMAX(KI,KJ)) GOTO 320
C Need to add other types to REBO interaction
C Change tight binding to >10
C  -travisk
           if((kj.le.10).and.(ki.le.10)) LCHECK(K)=1
           if((kj.ge.10).and.(ki.ge.10)) LCHECK(K)=2
C temperarly disable for test of CH only
C           if((kj.le.2).and.(ki.le.2)) LCHECK(K)=1
C           if((kj.ge.3).and.(ki.ge.3)) LCHECK(K)=2
           RC=SQRT(RSQ)
           rt = rc/ddtab(ki,kj)
           it = min(int(rt) + 1,ntab-1)

           RCOR(K)=RC
C
           WW(K)=TABFC(ki,kj,it)
     &            +(TABFC(ki,kj,it+1)-TABFC(ki,kj,it))*(rt-it+1)
           DWW(K)=TABDFC(ki,kj,it)
     &           +(TABDFC(ki,kj,it+1)-TABDFC(ki,kj,it))*(rt-it+1)
C
           EXX1(K) = atable(ki,kj,it)
     &              +(atable(ki,kj,it+1)-atable(ki,kj,it))*(rt-it+1)
           DEXX1(K) = datable(ki,kj,it) +
     &        (datable(ki,kj,it+1)-datable(ki,kj,it))*(rt - it +1)
C
           IF(I.GE.J) GO TO 320
C
           vv = rtable(ki,kj,it)
     &              +(rtable(ki,kj,it+1)-rtable(ki,kj,it))*(rt-it+1)
           rp = drtable(ki,kj,it)
     &              +(drtable(ki,kj,it+1)-drtable(ki,kj,it))*(rt-it+1)
           tote = tote + vv
           eatom(i) = eatom(i) + vv/2.0d0
           eatom(j) = eatom(j) + vv/2.0d0
!
! Print pair interaction
! -travisk
          IF (PRNPR.and.I.EQ.PR1.AND.J.EQ.PR2) THEN
            WRITE(23,2301) I,J,vv,EXX1(K)*2.0d0
            PREN = VV
            PRDIST = RC
          ENDIF

          DO 318 L=1,3
               RPP(K,L)=RP*RR(L)
318       CONTINUE
320   CONTINUE
C
      DO 321 K=1,KEND
           if(lcheck(k).eq.0) go to 321
           I=IVCT2B(K)
           J=JVCT2B(K)
           IF(I.GE.J) GO TO 321
            DO 322 L=1,3
                RNP(I,L)=RNP(I,L) + RPP(K,L)
                RNP(J,L)=RNP(J,L) - RPP(K,L)
322        CONTINUE
321   CONTINUE
c      if(noa(11)+noa(12).ne.0) call sili_germ
      CALL  cutfunc 
C
      RETURN
 2301 FORMAT(2I8,2F14.6)
C
      END
c

