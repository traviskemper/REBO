C Set variables for analysis 
C
      SUBROUTINE set_anly
C
      USE ANALYSIS
      USE PARAMS
      USE STRUCTR
      USE POTS   
      USE MDSTP
      USE BEAM
      USE SPECIF
      USE MPIvars
C
      IMPLICIT none
C
      INTEGER :: I,J,N,A,T,D,M,NMB,NMU,NM,ADUM,BDUM,DP
     &,TYPC(12),MT,SCNT,NPi,K,NATOM
     & ,NBJ,JBEGIN,JEND,IOS
      CHARACTER(30) :: CDUM,DDUM
      LOGICAL ::  depexists
!
      INCLUDE 'mpif.h'
!
C
      IF(mynode.EQ.0) THEN
      CALL molcheck
C     Structure
!     Initialize previous step information
      NNBo       = NNBi
      NABORSo(:) = NABORSi(:)
      LISTo(:)   = LISTi(:)
      MOLCNTo    = MOLCNTi
      MPNTo(:)   = MPNTi(:)
      MOLSTo(:)  = MOLSTi(:) 
      ASUBo(:)   = ASUBi(:)
      IMOLo(:)   = IMOLi(:)
      VDEP = R1i(DEPD,NP)
!
! Get max and min position of subtrate atoms 
      SMNo(1,:) = (/1.d16,1.0d16,1.0d16/)
      SMNo(2,:) = (/-1.0d16,-1.0d16,-1.0d16/)
      SCNT = 0
      MTOT = 0.0d0 
      DO I=1,NP
        IF ( ASUBo(I) ) THEN
          SCNT = SCNT + 1
          MTOT = MTOT + XMASS(KTYPEi(I) )
          DO N=1,3
           IF ( R0i(N,I) .LT. SMNo(1,N) ) SMNo(1,N)=R0i(N,I)
           IF ( R0i(N,I) .GT. SMNo(2,N) ) SMNo(2,N)=R0i(N,I)
          ENDDO 
        ENDIF
      ENDDO
      VOLMXN = (SMNo(2,1) - SMNo(1,1))*(SMNo(2,3) - SMNo(1,3))
     &         *(SMNo(2,3) - SMNo(1,3))
      MXNDENS = MTOT/VOLMXN*DCONV
!      
! Read in information from previous run
      IF ( RECAN ) THEN
        INQUIRE(file='out.react',exist=depexists)
!       Read in information about equilibrated substrate
        IF ( depexists) THEN
          WRITE(*,*) 'Read in out.react'
          OPEN(UNIT=27,FILE='out.react',STATUS='unknown')
          READ(27,'(A70)') DTITLE
          READ(27,*) 
          READ(27,*) SUBNA,SUBMN,SUBMX,SUBTPEN,SUBTKEN,SUBTEMP
          REC = 0
          DO
            READ(27,*,IOSTAT=ios) DDUM,CDUM,REC 
            IF ( ios .gt. 0 ) THEN
              WRITE(*,*) 'Something is wrong'
              EXIT
            ELSE IF (ios .lt. 0 ) THEN
              EXIT
            ENDIF
          ENDDO
          CLOSE(27)
!         Re initialize reaction file so it does not get too long
          OPEN(UNIT=27,FILE='out.react',STATUS='unknown')
          WRITE(27,'(A70)') DTITLE
          WRITE(27,*) 'Substrate info:,# atoms, min, max, PE, KE, TEMP'
          WRITE(27,2701) SUBNA,SUBMN,SUBMX,SUBTPEN,SUBTKEN,SUBTEMP
        ELSE
          WRITE(*,*) 'out.react does not appear to exist'
          WRITE(*,*) 'please run a single step of the 
     &                 equilibrated substrate with '
          WRITE(*,*) '       PRINTSUBINFO = .TRUE. '
          WRITE(*,*) 'to recorded substrate information'
          STOP 
        ENDIF
      ENDIF
      IF ( DEPPROF ) THEN
        OPEN(UNIT=30,FILE='out.depth',status='unknown')
        SUBCT = SUBMX + HBUF 
        SDEP = SUBCT - SUBMN
        SDIV = INT( SDEP/DPDIV  )
      ENDIF 
      ALLOCATE( 
     &   BRB(BTYPES),FRB(BTYPES),EXB(BTYPES)
     &   ,BRBDP(SDIV,BTYPES),FRBDP(SDIV,BTYPES),EXBDP(SDIV,BTYPES)
     &   ,BRBGS(BTYPES),FRBGS(BTYPES),EXBGS(BTYPES)
     &   ,BID(NTYPES,NTYPES)
     &   ,ELCOMP(2,NTYPES),GASKE(NTYPES)
     &   )
! debug 
      WRITE(*,*) mynode,'BTYPES',BTYPES
      WRITE(*,*) mynode,'SDIV',SDIV
      WRITE(*,*) mynode,'BRB',BRB
!     end debug
!
!     Set bond types for depth analysis
      BID(:,:)          =  1
      BID(icarb,icarb)  =  2
      BID(icarb,ihyd)   =  3
      BID(ihyd,icarb)   =  3
      BID(ihyd,ihyd)    =  4
      BID(isulfur,isulfur) = 5
      BID(icarb,isulfur)   = 6
      BID(isulfur,icarb)   = 6
      BID(isulfur,ihyd)    = 7
      BID(ihyd,isulfur)    = 7
      BID(ioxygen,ioxygen) = 8
      BID(icarb,ioxygen)   = 9
      BID(ioxygen,icarb)   = 9
      BID(ioxygen,ihyd)    = 10
      BID(ihyd,ioxygen)    = 10 
!
      ALLOCATE( GSA(NTYPES),SBA(NTYPES),DSBA(SDIV,NTYPES)
     &            ,MXDEP(NTYPES) )
C
!
      ENDIF
      RETURN
 2701 FORMAT(I10,2F16.8,3F14.3)
      END SUBROUTINE set_anly
!
!     Get information about substrate 
!     Get information about substrate 
      SUBROUTINE PRNSUBE 
!
      USE dyn_array
      USE ANALYSIS
      USE PARAMS
      USE STRUCTR
      USE MDSTP
      USE POTS
      USE BEAM
      USE SPECIF
      USE MPIvars
!
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
      INTEGER :: I,J,IOS,ierr 
      REAL*8 ::  XX,XXtmp
!
!  CALCULATE KINETIC ENERGY
!
      XX=0.0d0
      DO I=1,NP_node
           DO J=1,3
                XX=XX+(R1_node(J,I)**2)*XMASS(KTYPE_node(I))
                XXtmp=XX
           enddo
      enddo
           CALL MPI_REDUCE(XXtmp,XX,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
C      WRITE(*,*) mynode,'TOTE=',TOTE
           CALL MPI_REDUCE(TOTE,TTOTE,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)


!     Store substrate information 
      IF (mynode.EQ.0) THEN
          SUBNA = NP
          SUBMN = RMN(1,DEPD)
          SUBMX = RMN(2,DEPD)
          SUBTPEN = TOTE
          SUBTKEN = XX/(4.0d0*DELTSQ)*ECONV
          SUBTEMP = ENPR*XX/6.0d0/DELTSQ*ECONV
          DTITLE = CTITLE 
          OPEN(UNIT=27,FILE='out.react',STATUS='unknown')
          WRITE(27,'(A70)') DTITLE
          WRITE(27,*) 'Substrate info:,# atoms, min, max, PE, KE, TEMP'
          WRITE(27,2701) SUBNA,SUBMN,SUBMX,SUBTPEN,SUBTKEN,SUBTEMP
          CLOSE(27)
      ENDIF
!
      RETURN
 2701 FORMAT(I10,2F16.8,3F14.3)
      END SUBROUTINE PRNSUBE
