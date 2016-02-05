      subroutine setin
!
      USE MPIvars
      USE dyn_array
      USE PARAMS
      USE SPECIF 
      USE STRUCTR
      USE POTS
      USE TABS
      USE CONTM
      USE MDSTP
!
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
      INTEGER :: I,J,NTYPSQ
!
      ALLOCATE(AD(RTYPES,RTYPES)
     &,AXL(RTYPES,RTYPES)
     &,BD(RTYPES,RTYPES)
     &,BXL(RTYPES,RTYPES)
     &,CD(RTYPES,RTYPES)
     &,CXL(RTYPES,RTYPES)
     &,DD(RTYPES,RTYPES)
     &,DXL(RTYPES,RTYPES)
     &,ED(RTYPES,RTYPES)
     &,RB1(RTYPES,RTYPES)
     &,RB2(NTYPES,NTYPES)
     &,PID(RTYPES,RTYPES)
     &,RMAX(RTYPES,RTYPES),RMIN(RTYPES,RTYPES)
     &,CHI(RTYPES,RTYPES)
     &,RLIST(RTYPES,RTYPES)
     & )
      ALLOCATE(
     & DDTAB(RTYPES,RTYPES)
     &,tabfc(RTYPES,RTYPES,NTAB)
     &,tabdfc(RTYPES,RTYPES,NTAB)
     &,atable(RTYPES,RTYPES,NTAB)
     &,datable(RTYPES,RTYPES,NTAB)
     &,rtable(RTYPES,RTYPES,NTAB)
     &,drtable(RTYPES,RTYPES,NTAB)
     &)

C Allocate noa since it is initialized here
C -travisk
      ALLOCATE( NOA(NTYPES) )

C LJ arrays
      ALLOCATE(
     & XM(NTYPES,NTYPES)
     &,XMM(NTYPES,NTYPES)
     &,XMMS(NTYPES,NTYPES)
     &,RSPL(NTYPES,NTYPES)
     &,RMAXLJ(NTYPES,NTYPES)
     &,RSLJ(NTYPES,NTYPES)
     &,RSPLS(NTYPES,NTYPES)
     &,C2(NTYPES,NTYPES),C3(NTYPES,NTYPES)
     &,vlook(10000,NTYPES,NTYPES)
     &,dlook(10000,NTYPES,NTYPES)
     &,SIG(NTYPES,NTYPES)
     &,EPS(NTYPES,NTYPES)
     &,iv(NMABIG)
     &,jv(NMABIG)
     &)
      
      ALLOCATE(
     & TAU(NTYPES)
     &,EPSS(NTYPES,NTYPES)
     &,SIGS(NTYPES,NTYPES)
     &)


      KT(6)  = icarb
      KT(1)  = ihyd
      KT(8)  = ioxygen
      KT(16) = isulfur
      KT(9)  = iflor
      KT(10) = ineon
      KT(18) = iarg
      KT(36) = ikryp

C      KT(14) = 11
C      KT(32) = 12

      KT2(icarb) = 6
      KT2(ihyd) = 1 
      KT2(ioxygen) = 8
      KT2(isulfur) = 16
      KT2(iflor) = 9
      KT2(iarg) = 18
      KT2(ineon) = 10
      KT2(ikryp) = 36
!
      ATYPE(icarb)   = ' C'
      ATYPE(ihyd)    = ' H'
      ATYPE(ioxygen) = ' O'
      ATYPE(isulfur) = ' S'
      ATYPE(iflor)   = ' F'
      ATYPE(iarg)    = 'Ar'
      ATYPE(ineon)   = 'Ne'
      ATYPE(ikryp)   = 'Kr'
C      KT2(11) = 14
C      KT2(12) = 32


      xmass(:) = 0.0d0
      noa(:) = 0
      do i=1,NTYPES
          do j=1,NTYPES
               sig(i,j) = 0.0d0
               eps(i,j) = 0.0d0
          enddo
      enddo


C  Set the molecular masses
      XMASS(icarb)   = 12.0107  !Carbon
      XMASS(ihyd)    = 1.00794  !Hydrogen
      XMASS(ioxygen) = 15.9994  !Oxygen
      XMASS(iflor)   = 18.99840 !Flourine
      XMASS(isulfur) = 32.065   !Sulfur
      XMASS(ineon) = 20.1797d0  
      XMASS(iarg) =  39.948d0 
      XMASS(ikryp) = 83.80d0
C Set LJ parameters
      EPS(icarb,icarb) = 51.2d0 
      SIG(icarb,icarb) = 3.35d0 
      EPS(ihyd,ihyd) = 15.0d0 
      SIG(ihyd,ihyd) = 2.81d0
      EPS(ioxygen,ioxygen) = 80.5d0
      SIG(ioxygen,ioxygen) =  3.03d0
      EPS(isulfur,isulfur) = 80.4d0
      SIG(isulfur,isulfur) =   3.13d0
        EPS(iarg,iarg)       =  119.8d0
        SIG(iarg,iarg)       =  3.41d0 
!        EPS(ineon,ineon) =   47.0 
!        SIG(ineon,ineon) =   2.72 
!        EPS(ikryp,ikryp) =   164.0
!        SIG(ikryp,ikryp) =   3.83 
      KTMAX =  NTYPES 
!      NTYPSQ =RTYPES*RTYPES

      return 
      end 

      subroutine setpp


      USE MPIvars
      USE POTS
      USE SPECIF
      
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
      INTEGER :: I,J,K,L,M,N,I2D,IC,I3D,ITD,ierr
!     
      if(ipot.eq.1) then
           CALL PARAM
           call mtable
      endif       
!
      if(ILJ) THEN
        call ljparam
        call ljcset
      ENDIF
!
      return
      end

      subroutine setmd 

      USE MPIvars
      USE STRUCTR
      USE PARAMS
      USE SPECIF
      USE MDSTP

      IMPLICIT none 


      INCLUDE 'mpif.h'
c

            IF(NMA.ne.0) THEN
                  ENPR=EPSI/FLOAT(NMA)
            ELSE
                 ENPR=0.0D0
            ENDIF
C
            xkea = 0.0d0
            time = 0.0d0
            lchk = 1
      return 
      end 

! Create cfg files for atomeye
!
      SUBROUTINE setcfg
!
      USE MPIvars
      USE STRUCTR
      USE SPECIF
      USE ANALYSIS
      USE BEAM
!
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
      INTEGER :: N  
      REAL*8  :: DM
!
      IF(mynode.EQ.0) THEN
         N = KVC/NSTR
!     
         IF ( N.GT.999999) THEN
            WRITE(6,*) 'Warning more than 999999 cfg files are needed'
            WRITE(6,*) ' change filename format in cfg.f'
            STOP
         ENDIF
!       Set box size 
         DO N = 1,3 
            CBOX(N) = CUBE(N)
            DM = RMN(2,N) - RMN(1,N) 
            WRITE(6,*) 'act dimension',N,DM
            IF(CBOX(N).GT.DM+50.d0 ) THEN
               cbox(N) = DM*2.0d0
            ENDIF
            IF ( CBOX(N).LT.10.d0 ) THEN
               CBOX(N) = 10.d0
            ENDIF
         ENDDO
!        Create Cfgs folder
         CALL SYSTEM('sh -c "if [ ! -d Cfgs ]; then mkdir Cfgs ; fi"')
      ENDIF
!      
      RETURN
      END SUBROUTINE setcfg 
