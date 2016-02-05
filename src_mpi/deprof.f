!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Programer:
!     Travis Kemper
!     Department of Materials SCience and ENgineering
!     University of FLorida
!     traviskemper@ufl.edu
!
!     Version 1.0 04/04/11 T. W. Kemper
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Print depth profile of deposited atoms
!
      SUBROUTINE DEPTHPROF
!
      USE MDSTP
      USE MPIvars
      USE SPECIF
      USE ANALYSIS
      USE POTS
      USE STRUCTR
      USE BEAM 
!
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
      INTEGER :: i,j,D,KI
      REAL*8  :: DP,DPTH,DPTHDIV
      CHARACTER(12*SDIV) :: DS
!
      IF ( mynode.EQ.0 )THEN   
        GSA(:)    = 0
        SBA(:)    = 0
        DSBA(:,:) = 0
        MXDEP(:) = 0.0d0
        DO I=SUBNA+1,NP
          KI = KTYPEi(I)
          DP = R0i(DEPD,I)
          IF ( DP.GT.SUBMN.AND. DP.LT.SUBCT ) THEN
            DPTH =  SUBCT - DP
            DPTHDIV = DPTH/DPDIV + 1.0d0
            D = INT(DPTHDIV) 
            SBA(KI) = SBA(KI) + 1
            DSBA(D,KI) = DSBA(D,KI) + 1
            IF( DPTH.GT.MXDEP(KI) ) MXDEP(KI) = DPTH
          ELSE
            GSA(KI) = GSA(KI) + 1
          ENDIF
        ENDDO
        DO KI = 1,NTYPES
          WRITE(DS,*) DSBA(:,KI)
          WRITE(30,301) TTIME,LSTEP,KI,MXDEP(KI),GSA(KI),SBA(KI),DS
        ENDDO 
      ENDIF
 301  FORMAT(F16.2,I10,I4,F12.4,2I8,A) 
      END SUBROUTINE DEPTHPROF
