!
      subroutine setin
!
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
      INTEGER :: I,J
!     REBO arrays 
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
     &,PID(RTYPES,RTYPES)
     &,RMAX(RTYPES,RTYPES)
     &,CHI(RTYPES,RTYPES)
     &,RLIST(RTYPES,RTYPES)
     &,XHC(NPMAX,RTYPES)  !Number of neighbors of each type
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
C     LJ arrays
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
     &,iv(NMABIG)
     &,jv(NMABIG)
     &,RB2(NTYPES,NTYPES)
     &,SIG(NTYPES,NTYPES)
     &,EPS(NTYPES,NTYPES)
     &) 
      ALLOCATE(
     & TAU(NTYPES)
     &,EPSS(NTYPES,NTYPES)
     &,SIGS(NTYPES,NTYPES)
     &,NOA(NTYPES)
     &)
!    not sure about this allocation may cause problems
      ALLOCATE(
     & GL(3*NPMAX)
     &)
!
      ALLOCATE( EATOM(NPMAX)
     &,NABORS(NLMAX)          !Neighbor list
     &,LIST(NLMAX)            !N list pointer
     &,NABORSi(NLMAX),LISTi(NLMAX) !current step NB list
     &,RCLIST(NPMAX)            !Reaction #
     & ,IMOLo(NPMAX),IMOLi(NPMAX)
     & ,SIMOLo(NPMAX),SIMOLi(NPMAX)
     &,IVCT2B(NLMAX)
     &,JVCT2B(NLMAX)
     &,LCHECK(NLMAX)
     &,COR(NLMAX,DMS)
     &,RCOR(NLMAX)
     &,WW(NLMAX)
     &,DWW(NLMAX)
     &,EXX1(NLMAX)
     &,DEXX1(NLMAX)
     &)
!
!
!
!
      ALLOCATE(
     & KTYPE(NPMAX)
     &,R0(NPMAX,DMS) 
     &,R1(NPMAX,DMS)
     &,R2(NPMAX,DMS)
     &,R3(NPMAX,DMS)
     &,R4(NPMAX,DMS)
     &,ITR(NPMAX)
     &,MLIST(NPMAX)
     &,NLIST(NPMAX)
     &,R0L(NPMAX,DMS)
     &)
      ALLOCATE( RNP(NPMAX,DMS) 
     & )
 !
      KT(6)  = icarb
      KT(1)  = ihyd
      KT(8)  = ioxygen
      KT(16) = isulfur
      KT(9)  = iflor
      KT(10) = ineon
      KT(18) = iarg
      KT(36) = ikryp
!
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
!
      xmass(:) = 0.0d0
      noa(:) = 0
      do i=1,NTYPES
          do j=1,NTYPES
               sig(i,j) = 0.0d0
               eps(i,j) = 0.0d0
          enddo
      enddo
!
!  Set the molecular masses
      XMASS(icarb)   = 12.0107d0  !Carbon
      XMASS(ihyd)    = 1.00794d0  !Hydrogen
      XMASS(ioxygen) = 15.9994d0  !Oxygen
      XMASS(iflor)   = 18.99840d0 !Flourine
      XMASS(isulfur) = 32.065d0  !Sulfur
      XMASS(ineon) = 20.1797d0  
      XMASS(iarg) =  39.948d0 
      XMASS(ikryp) = 83.80d0
!
C      xmass(21) = 28.0d0
C      xmass(22) = 72.0d0
!     Initialize reaction # to 0 if used or not
      RCLIST(:) = 0
      IMOLo(:) = 0
      
!
      return 
      end 

      subroutine setpp

      USE POTS
      USE SPECIF
      USE PARAMS

      IMPLICIT none
 
      INTEGER :: I,J,K,L,M,N,I2D,IC,I3D,ITD
       
      if(ipot.eq.1) then
           CALL PARAM
           call mtable
      endif
      if(ilj) THEN
!       Set LJ parameters
        EPS(icarb,icarb)     = 51.2d0 
        SIG(icarb,icarb)     =  3.35d0 
        EPS(ihyd,ihyd)       = 15.0d0 
        SIG(ihyd,ihyd)       =  2.81d0
        EPS(ioxygen,ioxygen) = 80.5d0
        SIG(ioxygen,ioxygen) =  3.03d0
        EPS(isulfur,isulfur) = 80.4d0
        SIG(isulfur,isulfur) =  3.13d0
        EPS(iarg,iarg)       =  119.8d0
        SIG(iarg,iarg)       =  3.41d0 
        EPS(iflor,iflor)       =  62.4d0
        SIG(iflor,iflor)       =  2.81d0
!        EPS(ineon,ineon) =   47.0 
!        SIG(ineon,ineon) =   2.72 
!        EPS(ikryp,ikryp) =   164.0
!        SIG(ikryp,ikryp) =   3.83 
!       Set the number of lj types REBO types
        ktmax = NTYPES
        call ljparam
        call ljcset
      ENDIF
      return
      end 
!
      subroutine setmd 
!
      USE STRUCTR
      USE PARAMS
      USE SPECIF
      USE MDSTP
!
      IMPLICIT none 
!
!     Set heeting rate
      DELTMP = DT/1000.0d0*DELTA
!
            IF(NMA.ne.0) THEN
                  ENPR=EPSI/FLOAT(NMA)
            ELSE
                 ENPR=0.0D0
            ENDIF
!
            xkea = 0.0d0
            time = 0.0d0
            lchk = 1
            return 
            end 
