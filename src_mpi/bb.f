C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* This code built by John D. Zahrt, 28 July 1999
* Last Modified:  JAN, 2002(by IJang)
* Last Modified June 30, 2004 by S Pregler to include Ar
* Last Modified June 6,  2009 by Travis Kemper
*
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* This code, bb for beam build, requires an input file of
* type and coordinates for each atom in a (thermalized)
* cluster.  The user must define 7 parameters and 17 variables
* defined in two data statements.  These numbers define the beam
* and the ranges in which random variables are to be chosen.
* The input cluster has its center of mass determined and then
* it is set down on the z axis, normal to and in the center of 
* the crystal surface, a given distance above the surface.
* Copies of the cluster are made at prescribed periodic
* distances along the z axis.  Then each cluster is randomly
* rotated (by angles theta and phi) and translated by random
* small amounts from its original position.  Finally the beam(s)
* is oriented by rotations and translations.  The output file 
* is then written.  A document exists which descibes the code 
* in more detail.
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Input files:
C    mol.xyz    !Molecule xyz 
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  Subroutine to build beam 

      SUBROUTINE bb

      USE beam

c The value of "max" shouldn't be less than "nta".
      integer   na,nm,nb,nta
      real      xc,yc,tmp

C     Prameters:
*       na      the number of atoms in a molecule      
*	nm	the number of molecules per cluster
*	nc	the number of clusters per beam
*	nb	the number of beams
*       nta     the number of total atoms
*       xc      the x coordinate of the center of the crystal surface
*       yc      the y coordinate of the center of the crystal surface

        parameter(
     &   nm    = 1,
     &   nb    = 1,
     &   xc    = 0.0,
     &   yc    = 0.0)
      integer   ntyp(BAMX),i,j
      real      x(BAMX,3),com(3),pb(10),rc(7),m,tm

*  Beam line   lo  lr    thet1     phi1  dx1 dy1 thet2 phi2 dx2 dy2   
      data pb/10.0,20.0, 0.0,   0.0,  0.0,0.0, 0.0,  0.0,0.0,0.0/
*  Rand. clus. pl      pu     tl    tu        xv    yv    zv
      data rc/0.0,6.11, 0.0, 6.11,0.0,0.0,0.0/ 
! Switch beam direction to z for building (will be switched back and end)
      TMP = DEPA(DEPD)
      DEPA(DEPD) = DEPA(3)
      DEPA(3) = TMP
      pb(1) = BMHT
  
      call rdclus(na,nm,nta,ntyp,x)

      call cenofmass(1,na*nm,nta,ntyp,x,com,tm)

      CALL mvelocity(tm,pb)

      call bldbeam(na,nm,nb,nta,ntyp,x,com,pb,xc,yc)

      call rndconb(na,nm,nb,nta,ntyp,x,rc)

      call beamorient(na,nm,nb,nta,x,pb,xc,yc)

      call wrtcoord(nta,ntyp,x,tm)

      END SUBROUTINE bb 

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C RDCLUS reads the coordinates and the type of atoms of the starting cluster. 

      subroutine rdclus(na,nm,nta,ntyp,x)

      USE beam

      IMPLICIT none

      integer    na,nm,nta,ntyp(BAMX)
      real       x(BAMX,3),rijsq,dr
      integer    i,J,N
      REAL :: rijmsq
!
! Set locals should be eliminated 
      NA = NPM
      DO I =1,NA
        ntyp(i) = ITYPEM(i)
        x(I,:) = RM(I,:)
      ENDDO

C Calculate total number of atoms 
      NTA = NA*NC

C Find the maximum bond length in the molecule 
      rijmsq= 0.0
      DO I=1,NA-1
        DO J=I+1,NA
          rijsq =0.0
          DO N=1,3
            dr = x(I,N) - x(J,N)
            rijsq = rijsq + dr*dr 
          ENDDO
          IF (rijsq.gt.rijmsq) rijmsq = rijsq
        ENDDO
      ENDDO
      MOLMXR = SQRT(rijmsq)

C Modify the deposition area accordingly
      DO N=1,3
        DEPA(N) = DEPA(N) - MOLMXR
      ENDDO

      return
      end

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* CENOFMASS determines the center of mass of a cluster (both
* in MAIN and in BLDBEAM).  Only atoms of carbon and hydro-
* are currently allowed although the coding should be trans-
* parent enough to allow the user to add other atoms.

      subroutine cenofmass(ll,lu,nta,ntyp,x,com,tm)

      USE beam
      USE POTS

      integer    ll,lu,nta,ntyp(BAMX)
      real       x(BAMX,3),com(3)
      integer    i
      real       m,tx,ty,tz,tm

      tx = 0.0
      ty = 0.0
      tz = 0.0
      tm = 0.0
      do 100 i=ll,lu
         m = XMASS( KT( ntyp(i) ))
         if ( m .LT. 0.00001d0 ) THEN 
           WRITE(*,*) 'No mass for ',ntyp(i)
           STOP
	 end if
         tx = tx + m*x(i,1)
         ty = ty + m*x(i,2)
         tz = tz + m*x(i,3)
         tm = tm + m
  100 continue
      com(1) = tx/tm
      com(2) = ty/tm
      com(3) = tz/tm
C      print *,com(1), com(2), com(3)

      return
      end
      
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Calculate the velocity of molecule based on energy given in
C  mol.xyz. As well as the inter molecular spacing determend by
C  the buffer time between impacts

      SUBROUTINE mvelocity(tm,pb)

      USE beam

      INTEGER :: M
      REAL :: pb(10)

C Calculate beam velocity
      DO M=1,3
        VB(M)=-sqrt(2*EB(M)/ECON/tm)
      ENDDO
      IF ( INT(tm*100) .NE.INT(MOLMASS*100)) THEN      
        WRITE(*,*) 'Total molecular mass is off'
        WRITE(*,*) MOLMASS,tm
        STOP
      ENDIF 

C Calculate distance between molecules
      pb(2) = -1d0*TBUF*VB(DEPD) 
!
      END SUBROUTINE mvelocity
!
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* BEAMORIENT orients the beam(s) from the default 
* orientation of being normal to the crystal surface.
* The rotation angles and the translation vectors are
* contained in the data statement for PB. 
!
      subroutine beamorient(na,nm,nb,nta,x,pb,xc,yc)

      USE beam

      integer    na,nm,nb,nta
      real       x(BAMX,3),pb(10),xc,yc
      integer    napb,i,j
      real       ct,st,cp,sp,xd,yd,zd
!
      napb = na*nm*nc
* For each beam
      do 100 i=1,nb
         ct = cos(((pb(3+(i-1)*4)/180.0)*3.141593))
         st = sin(((pb(3+(i-1)*4)/180.0)*3.141593))
         cp = cos(((pb(4+(i-1)*4)/180.0)*3.141593))
         sp = sin(((pb(4+(i-1)*4)/180.0)*3.141593))
* The 2 rotations and the translation are done in 1 step.
	 do 90 j=1,napb
            xd = x(j+(i-1)*napb,1) - xc
            yd = x(j+(i-1)*napb,2) - yc
            zd = x(j+(i-1)*napb,3)
            x(j+(i-1)*napb,1)= cp*ct*xd-sp*yd-cp*st*zd+xc
            x(j+(i-1)*napb,2)= sp*ct*xd+cp*yd-sp*st*zd+yc
            x(j+(i-1)*napb,3)= st*xd         +ct*zd
   90    continue
  100 continue

      return
      end

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* BLDBEAM builds the beam from a given cluster.  The center
* of mass of the cluster is among the input variables and
* is used as the center of the cluster.  This point is dropped
* down on the z axis a distance pb(1) from the center of the
* crystal surface (xc,yc) and then repeated at the prescribed distance,
* pb(2).  If a second beam is present, the second beam is an
* exact copy of the first beam since they will later be oriented
* differently.

      subroutine bldbeam(na,nm,nb,nta,ntyp,x,com,pb,xc,yc)

       USE beam

      integer    na,nm,nb,nta,ntyp(BAMX)
      real       x(BAMX,3), com(3), pb(10), xc, yc
      integer    i,j,k,napc,napb

*  First translate first cluster to starting point

      do 100 i=1,na*nm
         x(i,1) = x(i,1) - com(1) + xc
	 x(i,2) = x(i,2) - com(2) + yc
         x(i,3) = x(i,3) - com(3) + pb(1)
  100 continue

*  Now repeat cluster along beam

      napc = na*nm          !Number of atoms in Beam
      do 200 i=1,nc-1       
         do 190 j=1,napc
	    ntyp(i*napc+j) = ntyp(j)
	    do 180 k=1,2
	       x(i*napc+j,k) = x(j,k)
  180       continue
            x(i*napc+j,3) = x(j,3) + i*pb(2)
  190    continue
  200 continue
 
      if (nb .eq. 2) then

         napb = na*nm*nc
	 do 300 i=1,napb
	    ntyp(i+napb) = ntyp(i)
	    do 290 j=1,3
	       x(i+napb,j) = x(i,j)
  290       continue
  300    continue

      end if
*       na      the number of atoms in a molecule      
*	nm	the number of molecules per cluster
*	nc	the number of clusters per beam
*	nb	the number of beams
*       nta     the number of total atoms
*       xc      the x coordinate of the center of the crystal surface
*       yc      the y coordinate of the center of the crystal surface
      return
      end

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C RNDCONB is designed to randomize the orientation and the 
C position of each cluster on the beam(s).  The center of mass
C coordinates of each cluster are determined in CENOFMASS.
C Random numbers (uniformly distributed from 0 to 1) then
C allow random rotation angles, dp and dt, and random transla-
C tions from the center of the crystal surface to be assigned.
C The two separate rotions and the translation are all done in
C one DO loop. 
      
      subroutine rndconb(na,nm,nb,nta,ntyp,x,rc)
  
      USE beam
 
C      IMPLICIT none

      integer    na,nm,nb,nta,ntyp(BAMX)
      real       x(BAMX,3),rc(7),m
      integer    i,j,im1nc,napci,napc,seed,sec(3),days(3)
      INTEGER*4 :: today(3),now(3)
      real       dp,dt,cp,sp,ct,st,dx,dy,dz,xd,yd,zd,com(3)
     
C Get seed from time 
C      CALL time(seed)
      CALL itime(sec)
      CALL idate(days)
C      CALL time_and_date(day)
C      CALL CPU_Time(sec)
      seed = sec(1)*3600 + sec(2)*60 + sec(3) + days(1)*30 + days(2)
      DO I=1,100
       M = ran(seed)
      ENDDO
!      
      napc = na*nm
      do 100 i=1,nc*nb
         im1nc = (i-1)*napc
         call cenofmass(1+im1nc,i*napc,nta,ntyp,x,com,m)
         dp = ran(seed)*(rc(2)-rc(1))
	 dt = ran(seed)*(rc(4)-rc(3))
	 cp = cos(dp)
	 sp = sin(dp)
	 ct = cos(dt)
	 st = sin(dt)
C Area of random beam
	 dx = (ran(seed)-0.5)*DEPA(1)
	 dy = (ran(seed)-0.5)*DEPA(2)  ! y and z still swapt need to fix
	 dz = 0.0
	 do 90 j=1,napc
	    xd = x(im1nc+j,1) - com(1)
	    yd = x(im1nc+j,2) - com(2)
	    zd = x(im1nc+j,3) - com(3)
	    x(im1nc+j,1)=cp*ct*xd-sp*yd-cp*st*zd + dx + com(1)
	    x(im1nc+j,2)=sp*ct*xd+cp*yd-sp*st*zd + dy + com(2)
	    x(im1nc+j,3)=  st*xd         +ct*zd + dz + com(3)
   90    continue
  100 continue

      return
      end

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C WRTCOORD writes the atom type (1=hydrogen, 6=carbon)
C and the x,y,z coordinates of each atom in the beam
C to a file.  The routine is straightforward.

      subroutine wrtcoord(nta,ntyp,x,tm)

      USE beam

      integer    nta,ntyp(BAMX)
      real       x(BAMX,3),tm
      INTeger    i,j

      DO i=1,nta
          BTYP(I) = ntyp(I)
          RB(i,1) = x(i,1)
          RB(i,2) = x(i,2) 
          RB(i,3) = x(i,3)
      ENDDO
! Switch z to desired beam direction
      DO I=1,NTA
        TMP        = RB(I,DEPD)
        RB(I,DEPD) = RB(I,3)
        RB(I,3)    = TMP
        IF( TMP.GT.BMH) BMH=TMP
      ENDDO

C Set local variables to global variables 
      MMS = tm
      NABM = nta
!
 1001 format(i3, 3e15.6)
 3301 FORMAT(A30,'with',I6,'clusters')
!
      return
      end
