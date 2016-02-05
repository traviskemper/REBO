!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Programer:
!     Travis Kemper
!     Department of Materials SCience and ENgineering
!     University of FLorida
!     traviskemper@ufl.edu
!
!     Version 1.0 9/15/09 T. W. Kemper
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Merge a beam of ions in xyz format with a substrate in coord.d format 
!
      SUBROUTINE mkbeam 
!
      USE STRUCTR
      USE beam
      USE POTS
      USE ANALYSIS
      USE SPECIF
      USE MPIvars
      USE dyn_array
! 
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
      INTEGER :: NAI,M,I,J,ETMP,K,D,II,NN
     & ,n_pass,nma_node,IERR,NPtot
      DOUBLE PRECISION :: MASS,RMNACT(2,3)
      CHARACTER(7) :: fin,fout,STITLE,ENR
      LOGICAL :: PAR
!
      IF (mynode.EQ.0) THEN
       WRITE(6,*) 'starting beam build'
! Active area to be deposited on
!     Find rigid  amd thermostat box sizes 
c$$$      DO D=1,3
c$$$        M=1 !min
c$$$        RIG(M,D) = RMN(M,D) + RBOX(M,D)
c$$$        IF( RBOX(M,D).LT.0.1d0) RIG(M,D)=RMN(M,D)-1.0d0 !Account for round off
c$$$        THM(M,D) = RIG(M,D) + TBOX(M,D)
c$$$        IF( TBOX(M,D).LT.0.1d0) THM(M,D)=RMN(M,D)-1.0d0 !Account for round off
c$$$        M=2 !max
c$$$        RIG(M,D) = RMN(M,D) - RBOX(M,D)
c$$$        IF( RBOX(M,D).LT.0.1d0) RIG(M,D)=RMN(M,D)+1.0d0 !Account for round off
c$$$        THM(M,D) = RIG(M,D) - TBOX(M,D)
c$$$        IF( TBOX(M,D).LT.0.1d0) THM(M,D)=RMN(M,D)+1.0d0 !Account for round off
c$$$      ENDDO
!     Set to cube for now no
      WRITE(6,*) '!!!!Warning!!!'
      WRITE(6,*) 'no boundry for beam!'
      DEPA(1) = CUBE(1) !- SZTHM(1) - 2*THBUF
      DEPA(2) = CUBE(2) !- SZTHM(2) - 2*THBUF
      DEPA(3) = CUBE(3) !- SZTHM(3) - 2*THBUF
!     Set 
      IF( CENTR ) THEN 
        BMHT = 0.50d0*RMN(2,DEPD) - 0.50d0*RMN(1,DEPD) + DBUF
      ELSE
        BMHT = RMN(2,DEPD) + DBUF
      ENDIF
! 
! Build beam
      CALL bb
!
      WRITE(6,*) 'beam built'
      WRITE(6,*)' beam atoms',NABM
      WRITE(6,*)' with velocity of',VB(DEPD),' A/fs'
! Read in coordinates of beam
      DO J=NP+1,NABM+NP
        K = BTYP(J-NP)
        KTYPEi(J) = KT(K)
        ITRi(J)  = 0
        DO M=1,3
          R0i(M,J) = RB(J-NP,M)
        ENDDO
        IF(ktypei(J).eq.0.OR.ktypei(J).GT.NTYPES) then
                write(6,*) 'unknown beam atom type',ktypei(J)
     &                     ,' for atom ',k
                include 'close.inc'
                stop
        ENDIF 
        DO M=1,3
          R1i(M,J) = VB(M)*DELTA
        ENDDO
        R2i(:,J) = (/0.0,0.0,0.0/)
        R3i(:,J) = (/0.0,0.0,0.0/)
       ENDDO
!
! Update total number of atoms NP
       NP = NP + NABM
       WRITE(6,*) 'Beam has been added total atoms:',NP
!! Calculate the number of steps need to have impact
       CALL beamtime
       WRITE(6,*) 'Steps to completion: ',STPS 
! Modify box size to accommodate beam
       IF( CUBE(DEPD).LT.BMH*2.0d0 ) THEN
!         CUBE(DEPD) = 2.0d0*BMH + 2.0d0*THBUF
        WRITE(6,*)'WARNING!'
        WRITE(6,*)'Box size in dimension',DEPD,
     &             ' is to small to accommodate beam'
        WRITE(6,*)'It should be modified to ',2.0d0*BMH + 2.0d0*THBUF
       ENDIF
      ENDIF !mynode=0

      RETURN
!
   61 FORMAT(4A)
  101 FORMAT('Dep of ',A5,'on ',A7,'at ',F6.1,'eV')
!
      END SUBROUTINE mkbeam 
