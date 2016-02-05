      SUBROUTINE LJVLIST

      USE MPIvars
      USE dyn_array
      USE STRUCTR
      USE TABS
      USE MDSTP

      IMPLICIT none

      INCLUDE 'mpif.h'

      INTEGER :: i,j,k,L,KI,KJ
     &,cell,kli,icell
      REAL*8 ncell(3),RR,RSQS
C
c Find the LJ neighbor list
c


C      INTEGER  Kli_cell(3000)
      CALL ljlink
C      WRITE(*,*) mynode,'after ljlink'
C      WRITE(*,*)'ljvlist',mynode,NP_node,NPtot4_node
      ljvlit=0
      kli=0
C -travisk
      OPEN(71,file='lj.dat',status='unknown')
      WRITE(71,*) ljhead_node(:)
      WRITE(*,*)mynode,'NP_node',NP_node
      do i=1,NP_node
C          IF (mynode.EQ.0) THEN
C            WRITE(22,*) mynode,LSTEP,i
C          ELSE
C            WRITE(23,*) mynode,LSTEP,i
C          ENDIF
          ljnabr(I)=kli+1
          icell=n_icell_lj(I)
           DO k=1,27
             cell=ljnbc((icell-1)*27+k)
            j=ljhead_node(cell)
            DO WHILE (j.NE.0)
              IF (I.GE.J) GOTO 311
C           do 311 j=i+1,np
                KI=KTYPE_node(I)
                KJ=KTYPE_node(J)
C               ki = ktype(NA(i))
C               kj = ktype(NA(j))
               if(RSLJ(ki,kj).eq.0.0d0) go to 311
               RSQS=0.0D0
c
               DO L=1,3
                    RR=R0_node(L,I)-R0_node(L,J)
                    RR = RR - CUBE(L)*ANINT(RR/CUBE(L))
                    RSQS=RSQS+rr*rr
                    if(rsqs.gt.RSLJ(ki,kj)) go to 311
               enddo
               IF(RSQS.lt.xmms(ki,kj)) GO TO 311
               kli=kli+1
               WRITE(*,*) kli
               ljvlit(kli)=J
C               iv(kli)=i
C               jv(kli)=j
C 311       continue
311            j=ljlist_node(j)
            ENDDO
          ENDDO
      enddo
      ljnabr(NP_node+1)=kli+1
C     CALL MPI_REDUCE(kli,KLITO,1,MPI_INTEGER,MPI_SUM,0,
C    &MPI_COMM_WORLD,ierr)
      kliend=kli
C      WRITE(*,*) mynode,'LSTEP=',LSTEP,'  kli=',kli
      if(kliend.gt.(3*nmabig)) then
               write(*,*) 'LJ neighbor array= ',kliend
               write(*,*) 'increase NMABIG in allocate_arrays.f'
               CALL write_error
      endif
C      endif
C      INTEGER  Kli_cell(3000)
      CALL ljlink
C      WRITE(*,*) mynode,'after ljlink'
C      WRITE(*,*)'ljvlist',mynode,NP_node,NPtot4_node
      ljvlit=0
      kli=0
      do i=1,NP_node
C          IF (mynode.EQ.0) THEN
C            WRITE(22,*) mynode,LSTEP,i
C          ELSE
C            WRITE(23,*) mynode,LSTEP,i
C          ENDIF
          ljnabr(I)=kli+1
          icell=n_icell_lj(I)
          DO k=1,27
            cell=ljnbc((icell-1)*27+k)
            j=ljhead_node(cell)
            DO WHILE (j.NE.0)
              IF (I.GE.J) GOTO 311
C           do 311 j=i+1,np
                KI=KTYPE_node(I)
                KJ=KTYPE_node(J)
C               ki = ktype(NA(i))
C               kj = ktype(NA(j))
               if(RSLJ(ki,kj).eq.0.0d0) go to 311
               RSQS=0.0D0
c
               DO L=1,3
                    RR=R0_node(L,I)-R0_node(L,J)
                    RR = RR - CUBE(L)*ANINT(RR/CUBE(L))
                    RSQS=RSQS+rr*rr
                    if(rsqs.gt.RSLJ(ki,kj)) go to 311
               enddo
               IF(RSQS.lt.xmms(ki,kj)) GO TO 311
               kli=kli+1
               ljvlit(kli)=J
C               iv(kli)=i
C               jv(kli)=j
C 311       continue
311            j=ljlist_node(j)
            ENDDO
          ENDDO
      enddo
      ljnabr(NP_node+1)=kli+1
C     CALL MPI_REDUCE(kli,KLITO,1,MPI_INTEGER,MPI_SUM,0,
C    &MPI_COMM_WORLD,ierr)
      kliend=kli
      WRITE(*,*) mynode,'LSTEP=',LSTEP,'  kli=',kli
      if(kliend.gt.(3*nmabig)) then
               write(*,*) 'LJ neighbor array= ',kliend
               write(*,*) 'increase NMABIG in allocate_arrays.f'
               CALL write_error
      endif
C      endif


      RETURN
      END

      SUBROUTINE LJGUTS

      USE dyn_array
      USE MDSTP
      USE MPIvars
      USE TABS
      USE POTS
      USE PARAMS
      USE STRUCTR

      IMPLICIT none

      INCLUDE 'mpif.h'

      INTEGER :: I,J,K,L,M,N,KI,KJ,IT,kli,llj,II
     &,cell,jmin,jmax,ljlist,ljhead
      REAL*8 :: RR,RI(3)
      REAL*8 :: RLIS,RSQ,RC,RT,VV,RP
     &,R,vdw,dvdw,RSQS,dr

c
300   format(2i10,3f15.5)
      pvdw=0.D0
      RNP_lj=0.0d0
      do 323 I=1,NP_node !llj=1,kliend
        jmin=ljnabr(i)
        jmax=ljnabr(i+1)-1
        DO 322 k=jmin,jmax
           J=ljvlit(k)
C           i=iv(llj)
C           j=jv(llj)
           KI=KTYPE_node(I)
           KJ=KTYPE_node(J)
c
           RSQS=0.0D0

           DO L=1,3
               RRS(L)=R0_node(L,I)-R0_node(L,J)
               RRS(L) = RRS(L) - CUBE(L)*ANINT(RRS(L)/CUBE(L))
               RSQS=RSQS+RRS(L)*RRS(L)
               IF(RSQS.GT.RMAXLJ(ki,kj)) go to 322
           ENDDO
C
           if(RSQS.LT.XMM(KI,kj)) go to 322
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
            IF (j.GT.NP_node) THEN
              pvdw=pvdw+0.5d0*vdw
            ELSE
              pvdw=pvdw+vdw
            ENDIF
c
            do l=1,3
                 dr = dvdw*rrs(l)
                 RNP_lj(l,I)=RNP_lj(l,I) + dr
                 RNP_lj(l,J)=RNP_lj(l,J) - dr
            enddo
 322    CONTINUE
 323  continue

      tote = tote + pvdw

C       IF (mynode.EQ.0 ) WRITE(*,*) 'LJGUTS finished',RNP(2,2),tote
C      WRITE(*,*) 'lj',TOTE

      return
      end

