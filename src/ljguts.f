      SUBROUTINE LJGUTS

      USE SPECIF
      USE MDSTP
      USE STRUCTR
      USE POTS
      USE PARAMS
      USE TABS

      IMPLICIT none

      INTEGER :: I,J,K,L,M,N,KI,KJ,IT,kli,llj,II
      REAL*8 :: RPP(nlmax,3),RR,RI(3)
      REAL*8 :: RLIS,RSQ,RC,RT,VV,RP
     &,R,vdw,dvdw,RSQS,dr

C
c Find the LJ neighbor list
c
      IF(LCHK.EQ.1) THEN
      kli=0
      do i=1,np-1
          do 311 j=i+1,np
               ki = ktype(i)
               kj = ktype(j)
               if(RSLJ(ki,kj).eq.0.0d0) go to 311
               RSQS=0.0D0
c
               DO L=1,3
                    RR=R0(I,L)-R0(J,L)
                    RR = RR - CUBE(L)*ANINT(RR/CUBE(L))
                    RSQS=RSQS+rr*rr
                    if(rsqs.gt.RSLJ(ki,kj)) go to 311
               enddo
               IF(RSQS.lt.xmms(ki,kj)) GO TO 311
               kli=kli+1
               iv(kli)=i
               jv(kli)=j
 311       continue
      enddo

      kliend=kli
      if(kliend.gt.nmabig) then
               write(*,*) 'LJ neighbor array= ',kliend
               write(*,*) 'increase nmabig and recompile'
               stop 
      endif
      endif
c
300   format(2i10,3f15.5)
      pvdw=0.D0
      do 323 llj=1,kliend
           i=iv(llj)
           j=jv(llj)
           ki=ktype(i)
           kj=ktype(j)
c
           RSQS=0.0D0

           DO L=1,3
               RRS(L)=R0(I,L)-R0(J,L)
               RRS(L) = RRS(L) - CUBE(L)*ANINT(RRS(L)/CUBE(L))
               RSQS=RSQS+RRS(L)*RRS(L)
               IF(RSQS.GT.RMAXLJ(ki,kj)) go to 323
           ENDDO
C
           if(RSQS.LT.XMM(KI,kj)) go to 323
c table look up with linear interpolation
                 r = sqrt(rsqs) - xm(ki,kj)
                 rt = r/dellj
                 ii = int(rt)+1
                 vdw  = vlook(ii,ki,kj) +
     &                  (vlook(ii+1,ki,kj)-vlook(ii,ki,kj))*(rt-ii+1)
                 dvdw = dlook(ii,ki,kj) +
     &                  (dlook(ii+1,ki,kj)-dlook(ii,ki,kj))*(rt-ii+1)


c$            endif
c
            pvdw=pvdw+vdw
c
            do l=1,3
                 dr = dvdw*rrs(l)
                 rnp(i,l)=rnp(i,l) + dr
                 rnp(j,l)=rnp(j,l) - dr
            enddo
 323  continue

      tote = tote + pvdw
      return
      end

