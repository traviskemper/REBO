C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C This program move all gas phase molecules with non depositionC
C velocities under the substrate so they can continue to used  C
C in analysis, but will not effect further depositions         C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C     Programer:                                               C
C      Travis Kemper                                           C
C      traviskemper@ufl.ed                                     C
C    Advisor:                                                  C
C      Dr. Susan Sinnott                                       C
C      ssinn@mse.ufl.edu                                       C
C      Department of Materials Science and Engineering         C
C      University of FLorida                                   C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C     Version 1.0 06/29/2010 T. W. Kemper                      C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C

      SUBROUTINE cleargas

      USE ANALYSIS
      USE STRUCTR
      USE BEAM
      USE POTS
      USE MPIvars
      USE MDSTP
      USE SPECIF
!
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
      INTEGER :: I,M,Io,Jo,Oo,IM,N
     &,ECNT,ELIST(STHS),ELCNT(ELMX)
      REAL*8 :: R0MX,R0MN,MSHIFT,GBUF
     &,COM(3),MOLM,MSZ,MXCG
     &,SUBMXi,SMN(2,3),MVEL(3)
!
      IF ( mynode.EQ.0 ) THEN
      GBUF = RMAX(icarb,icarb)
!
!     Get max and min position of subtrate atoms and
      SMN(1,:) = SUBMN
      SMN(2,:) = SUBMX
      DO I=1,NP
         IF ( ASUBi(I) ) THEN
          DO N=1,3
            IF ( R0i(N,I) .GT. SMN(2,N) ) SMN(2,N)=R0i(N,I)
            IF ( R0i(N,I) .LT. SMN(1,N) ) SMN(1,N)=R0i(N,I)
          ENDDO
         ENDIF
        IF ( R0i(DEPD,I).LT. GASD ) GASD=R0i(DEPD,I)
      ENDDO
!
      SUBMXi = SMN(2,DEPD)
      GASD = GASD - GBUF
      MXCG = -CUBE(DEPD)/2.0d0
      WRITE(*,*) 'Sub max',SUBMXi
      WRITE(*,*) 'Gas depth=',GASD
      WRITE(*,*) 'Dep V=',VDEP
      WRITE(*,*) 'Dep dimension',DEPD
      WRITE(*,*) 'G buffer',GBUF
C
      DO 1001 M=1,MOLCNTi
        Io = MPNTi(M)
        Jo = MPNTi(M+1)-1
        ECNT = 0
        R0MX = - CUBE( DEPD )
        R0MN =   CUBE( DEPD )
        DO Oo=Io,Jo
          I=MOLSTi(Oo)
          IF ( ASUBi(I) ) GOTO 1001
          ECNT = ECNT + 1
          ELIST(ECNT) = I
          IF( R0i(DEPD,I).GT.R0MX ) R0MX=R0i(DEPD,I)
          IF( R0i(DEPD,I).LT.R0MN ) R0MN=R0i(DEPD,I)
        ENDDO
        CALL COMAS(ECNT,COM,MOLM,MVEL)
        WRITE(*,*) 'Molecule M=',M,COM(DEPD),Jo-Io+1
        WRITE(*,*) I,R0i(DEPD,I),R1(DEPD,I)
        IF (R0MN.GT.SMN(2,DEPD)+R0MX-R0MN ) THEN !.AND.R1(DEPD,I).NE.VDEP ) THEN
          MSHIFT = R0MX - GASD + GBUF
          DO Oo=Io,Jo
            IM=MOLSTi(Oo)
            R0i(DEPD,IM) = R0i(DEPD,IM) - MSHIFT
            R1i(:,IM) = (/0.0,0.0,0.0/)
            R2i(:,IM) = (/0.0,0.0,0.0/)
            R3i(:,IM) = (/0.0,0.0,0.0/)
          ENDDO
          MSZ = R0MX - R0MN + GBUF
          GASD = GASD - MSZ
          WRITE(*,*) 'Molecule',M,' with last atom',IM
          WRITE(*,*) ' has been cleared with: New GASD,by,with size'
          WRITE(*,*) GASD,MSHIFT,MSZ
          WRITE(*,*) 'Max Min',R0MX,R0MN
          WRITE(*,*) 'New pos',R0i(DEPD,IM)
          IF ( GASD.LT.MXCG) THEN
            WRITE(*,*) 'Warning gas molecules have extended beyond box
     &     Please reoganize them more compactly under substrate and
     &      reset GASD accordingly'            
            WRITE(*,*) GASD,MXCG
          ENDIF
        ENDIF
 1001 CONTINUE 
!     Print out coord.d file 
        OPEN(UNIT=31,FILE='sub.d',STATUS='unknown')
        OPEN(UNIT=32,FILE='out.rclist',STATUS='unknown')
        WRITE(31,1800) CTITLE
        WRITE(31,100) NP
        WRITE(31,300) TTIME,DELTA
        WRITE(31,300) (CUBE_r(N),N=1,3)
        WRITE(31,300) (CUBE(N),N=1,3)
        DO I=1,NP
           WRITE(31,350) I,KT2(KTYPEi(I)),
     &     (R0i(N,I),N=1,3),ITRi(i),RCLIST(I)
           WRITE(31,360) I,((R1i(N,I)/DELTA),N=1,3)
           WRITE(31,360) I,((R2i(N,I)/DELTSQ),N=1,3)
           WRITE(31,360) I,(R3i(N,I),N=1,3)
        ENDDO
        CLOSE(31)
      ENDIF
      RETURN
  100 FORMAT(I7)
  300 FORMAT(3E20.11)
  350 FORMAT(2I5,3E20.11,I3,I8)
  360 FORMAT(I7,3E20.11)
  110 FORMAT(I7,3E20.11,I8)
 1800 FORMAT(20A2)
      END SUBROUTINE cleargas
