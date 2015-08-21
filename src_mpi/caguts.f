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

      SUBROUTINE colist
C
C CALCULATE TWO-BODY FORCES AND NEIGHBOR LIST for hydrocarbons
C
      USE MPIvars
      USE dyn_array
      USE STRUCTR
      USE MDSTP
      USE PARAMs
      USE POTS

      IMPLICIT none
c
      INCLUDE 'mpif.h'
c
      REAL*8 ncell(3),RR(3),RI(3),RSQ,RLIS
      INTEGER  :: I,J,K,L,KI,KJ,KK
     & ,cell,icell
C
c	write(273,*) lchk
C       IF(LCHK.EQ.1) THEN
C
C Set up neighbor list
C
C Sort atoms into the cell
       DO icell=1,nncell(1)*nncell(2)*nncell(3)
         nchead_node(icell)=0
       ENDDO
       DO i=1,3
         ncell(i)=float(nncell(i))
       ENDDO
       DO i=1,NPtot3_node
         icell=1+INT(((r0_node(1,i)/(CUBE(1)+1.0d-12))+0.5)*ncell(1))
     &          +INT(((r0_node(2,i)/(CUBE(2)+1.0d-12))+0.5)*ncell(2))
     &          *nncell(1)
     &          +INT(((r0_node(3,i)/(CUBE(3)+1.0d-12))+0.5)*ncell(3))
     &          *nncell(2)*nncell(1)
         nclist_node(i)=nchead_node(icell)
         nchead_node(icell)=i
         n_icell_node(i)=icell
       ENDDO
C Start to make the neighbor list by link cell technique
           DO I=1,NLMAX
             LIST(I)=0
           ENDDO
           K=0
           DO 302 I=1,NPtot_node
             NABORS(I)=K+1
             DO 299 L=1,3
               RI(L)=R0_node(L,I)
299          CONTINUE
             KI=KTYPE_node(I)
C
c cuts out all but C,H,Si, and Ge
C
             IF(ki.GE.6) go to 302
             icell=n_icell_node(I)
             DO kk=1,27
               cell=icneigh((icell-1)*27+kk)
               j=nchead_node(cell)
               DO while (j.NE.0)
                 IF(I.EQ.J) GO TO 301
                 KJ=KTYPE_node(J)
C
c cuts out all but C,H,Si, and Ge
C
                 IF(kj.GE.6) go to 301
                 RLIS=RLIST(KI,KJ)
C	write(274,*)'rlist=',ki,kj,rlist(ki,kj)
C
                 RSQ=0.0D0
                 DO 298 L=1,3
                   RR(L)=RI(L)-R0_node(L,J)
                   RR(L)=RR(L)-CUBE(L)*ANINT(RR(L)/CUBE(L))
                   RSQ=RSQ+RR(L)*RR(L)

c	write(274,*)'rsq, rlis',rsq,rlis,ki,kj

                   IF(RSQ.GT.RLIS) GO TO 301
298              CONTINUE
                 IF (RSQ.EQ.0) THEN
                   WRITE(*,*) 'Zero distance in colist',NA(I),NA(J)
                   CALL write_error
                 ENDIF
                 K=K+1
                 LIST(K)=J

C                 IVCT2B(K)=I
C                 JVCT2B(K)=J
301              j=nclist_node(j)
               ENDDO
             ENDDO
302        CONTINUE
           NABORS(NPtot_node+1)=K+1
C           WRITE(*,*) k
           KEND=K
C           CALL MPI_REDUCE(kend,KENDTO,1,MPI_INTEGER,MPI_SUM,0,
C     &     MPI_COMM_WORLD,ierr)
C           WRITE(*,*) mynode,'NPtot_node',NPtot_node
C           WRITE(*,*) mynode,'LSTEP=',LSTEP,' kend=',kend
           if(kend.gt.nlmax) then
                 write(*,*) 'kend exceeds nlmax'
                 write(*,*) 'kend,nlmax = ',kend,nlmax
                 write(*,*) 'increase NLMAX in allocate_arrays.f'
                 CALL write_error
           endif

      RETURN
C
      END
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
      USE dyn_array
      USE MPIvars
      USE SPECIF
      USE MDSTP
      USE STRUCTR
      USE POTS
      USE TABS

      IMPLICIT none

      INTEGER :: I,J,K,L,M,N,KI,KJ,IT,JMAX,JMIN
     &,ierr
      REAL*8 :: RR(3),RI(3),RLIS,RSQ,RC,RT,VV,RP

!
      INCLUDE 'mpif.h'
!
C      do i=1,NP
C           eatom(i) = 0.0d0
C      enddo
c
C      WRITE(*,*) 'LSTEP=',LSTEP,' kend=',kend
      DO 320 I=1,NPtot_node   !K=1,KEND
C
C           I=IVCT2B(K)
        jmin=nabors(i)
        jmax=nabors(i+1)-1
        DO 319 k=jmin,jmax
          J=list(k)
C           J=JVCT2B(K)
           KI=KTYPE_node(I)
           KJ=KTYPE_node(J)
C
           LCHECK(K)=0
           RSQ=0.0D0
           DO L=1,3
                RR(L)=(R0_node(L,I)-R0_node(L,J))-
     &          CUBE(L)*ANINT((R0_node(L,I)-R0_node(L,J))/CUBE(L))
                RSQ=RSQ+RR(L)*RR(L)
                COR(K,L)=RR(L)
           ENDDO
           IF (RSQ.EQ.0) THEN
             WRITE(*,*) 'Zero distance in caguts',NA(I),NA(J)
             CALL write_error
           ENDIF
c
           IF(RSQ.GT.RMAX(KI,KJ)) GOTO 319
C Add oxygen
C  -travisk
           if((kj.le.10).and.(ki.le.10)) LCHECK(K)=1
           if((kj.ge.11).and.(ki.ge.11)) LCHECK(K)=2

C Just CHF for debuging
!           if((kj.le.3).and.(ki.le.3)) LCHECK(K)=1
!           if((kj.ge.4).and.(ki.ge.4)) LCHECK(K)=2

c	write(273,*)lcheck(k),k

c	write(274,*)'ki,kj=',ki,kj

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
C          IF (I.GT.NP_node) GO TO 319
          IF(I.GT.NP_node) GO TO 319
          IF(I.GE.J) GOTO 319
C
          vv = rtable(ki,kj,it)
     &              +(rtable(ki,kj,it+1)-rtable(ki,kj,it))*(rt-it+1)
          rp = drtable(ki,kj,it)
     &              +(drtable(ki,kj,it+1)-drtable(ki,kj,it))*(rt-it+1)

          IF (J.GT.NP_node) THEN
           tote = tote + 0.5*vv
          ELSE
           tote=tote+vv
          ENDIF

c	write(273,*) vv,ki,kj,rtable(ki,kj,it)
c	write(273,*) exx1(k),k,atable(ki,kj,it)

C          eatom(NA(i)) = eatom(NA(i)) + vv/2.0d0
C          eatom(NA(j)) = eatom(NA(j)) + vv/2.0d0
C Print pair interaction
C -travisk
C          IF ( mynode.EQ.0) THEN
C            IF (PRNPR.and.I.EQ.PR1.AND.J.EQ.PR2) THEN
            IF (PRNPR.and.NA(I).EQ.PR1.AND.NA(J).EQ.PR2) THEN
              WRITE(23,2301) NA(I),NA(J),vv,EXX1(K)*2.0d0
              PREN = VV
              PRDIST = RC
C      CALL MPI_BCAST(PREN,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
C      CALL MPI_BCAST(PRDIST,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
            ENDIF
C          ENDIF

          DO 318 L=1,3
               RPP(L,K)=RP*RR(L)
318       CONTINUE
319     CONTINUE
320   CONTINUE
C
      DO 321 I=1,NP_node !K=1,KEND
        jmin=nabors(i)
        jmax=nabors(i+1)-1
        DO 322 k=jmin,jmax
          J=list(k)
          if(lcheck(k).eq.0) go to 322
          IF (I.GE.J) GOTO 322
C           I=IVCT2B(K)
C           J=JVCT2B(K)
C          IF(NA(I).GE.NA(J)) GO TO 322
          DO 323 L=1,3
C            RNP_node(L,I)=RNP_node(L,I)+RPP(L,K)
C            RNP_node(L,J)=RNP_node(L,J)-RPP(L,K)
            RNP(L,NA(I))=RNP(L,NA(I)) + RPP(L,K)
            RNP(L,NA(J))=RNP(L,NA(J)) - RPP(L,K)
323       CONTINUE
322     CONTINUE
321   CONTINUE

C      IF (mynode.EQ.0 ) WRITE(*,*) 'Caguts finished',RNP(2,2),tote

      RETURN
 2301 FORMAT(2I8,2F14.6)
C
      END
c

