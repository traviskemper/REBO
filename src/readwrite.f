C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C     Version 1.0 01/16/2010 T. W. Kemper                      C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C  Subroutines involved in reading and writing data
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  Write header of output file
!
      SUBROUTINE WRTHEAD
!
      USE SPECIF
      USE TM
      USE PARAMS
!
      IMPLICIT none
!
      INTEGER :: sec(3),day1,day2,day3
!
      WRITE(6,*) ''
     &,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     &,'+                  REBO  (CH/O/F/S)                         +'
     &,'+ THIS PROGRAMS PERFORM MOLECULAR DYNAMICS SIMULATIONS WITH +'
     &,'+ SHORT RANGE BOND-ORDER POTENTIALS AND LONG RANGE          +'
     &,'+ LENNARD-J0NES POTENTIALS FOR HYDROGEN, CARBON, OXYGEN     +'
     &,'+ SULFUR AND FLUORIEN                                       +'
     &,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     &,'+  Units: mass = AMU s, length = Angstroms, time = fs       +'
     &,'+         energy = eV s                                     +'
     &,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     &,'+                                                           +'
     &,'+  Details about development and properties of potential    +'
     &,'+  is described in :                                        +'
     &,'+  1. CH                                                    +'
     &,'+     D.W. Brenner+, O.A. Shenderova,                       +'
     &,'+     J.A. Harrison, S.J. Stewart, B. Ni, S.B. Sinnott,     +'
     &,'+     Journal of Physics, Condensed Matter 14,              +'
     &,'+     783-802 (2002).                                       +'
     &,'+  2. CHO                                                   +'
     &,'+     Boris Ni, Ki-Ho Lee and Susan B. Sinnott              +'
     &,'+     J. Phys. Condens. Matter 16 (2004) 7261-7275          +'
     &,'+  3. CHF                                                   +'
     &,'+     Inkook Jang and Susan B. Sinnott*                     +'
     &,'+     J. Phys. Chem. B, 2004, 108 (49), pp 18993-19001      +'
     &,'+  4. CHS                                                   +'
     &,'+     Travis Kemper and Susan B. Sinnott                    +'
     &,'+     J. Phys. Chem. C, 2011, 115 (48), pp 23936-23945      +'
     &,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     &,'+    Version Serial 1.00                                    +'
     &,'+       07/11/2011 T. W. Kemper: traviskemper@ufl.ed        +'
     &,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

C Calculate current date 
      !call idate(day1,day2,day3)   ! today(1)=day, (2)=month, (3)=year
      !call itime(sec)     ! now(1)=hour, (2)=minute, (3)=second
      !WRITE(6,102) day1,day2,2000+day3,sec
      WRITE(6,*)
     &,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE(6,*)''
!
 101  FORMAT(a16,3(i3,a8,1x),f6.1,1x,a8)
 102  FORMAT( ' + Date: ', i2.2, '/', i2.2, '/', i4.4, '; time ',
     & i2.2, ':', i2.2, ':', i2.2, '                           +' )
!
      END SUBROUTINE WRTHEAD
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Read in settings for simultion from input file
!  input file must specify variables with the following format:
!   [VAR] = [VALUE]
!  and is blank line terminated
!
      SUBROUTINE READIN 
!
      USE SPECIF
      USE STRUCTR
      USE POTS
      USE STRUCTR
!   
      IMPLICIT none
!
      INTEGER :: IOS,N
      CHARACTER(20) :: var
      CHARACTER(2) :: dumb
      REAL*8 :: val
!
C     Read in input file
      WRITE(13,*)'Reading input.d:'
      READ(12,*,IOSTAT=ios) var
      BACKSPACE(12) 
      DO WHILE( IOS.EQ.0)
C     Find variable then go back a space and read it in
        IF ( var(1:6).eq.'SYSTEM' ) THEN
         READ(12,'(A,A1,A70)') var,dumb,CTITLE
        ELSEIF  ( var(1:12).eq. 'LENNARDJONES' ) THEN
         READ(12,*) var,dumb,ILJ
        ELSEIF ( var(1:5).eq.'STEPS' ) THEN
         READ(12,*) var,dumb,KVC
        ELSEIF  ( VAR(1:7) .eq. 'TEMPMAX' ) THEN
          READ(12,*) var,dumb,TMAX
        ELSEIF  ( var(1:4).eq. 'TEMP' ) THEN
          READ(12,*) var,dumb,TEM
        ELSEIF  ( VAR(1:4) .eq. 'WOUT' ) THEN
          READ(12,*) var,dumb,MAXKB
        ELSEIF  ( VAR(1:11) .eq. 'THERMOSTATE' ) THEN
          READ(12,*) var,dumb,KFLAG
        ELSEIF  ( VAR(1:4) .eq. 'SEED' ) THEN
          READ(12,*) var,dumb,PSEED
        ELSEIF  ( VAR(1:10) .eq. 'INTEGRATOR' ) THEN
          READ(12,*) var,dumb,INTG
        ELSEIF  ( VAR(1:4) .eq. 'HEAT' ) THEN
          READ(12,*) var,dumb,HT
        ELSEIF  ( VAR(1:6) .eq. 'DELTAT' ) THEN
          READ(12,*) var,dumb,DT
        ELSEIF  ( VAR(1:13) .eq. 'SPACIALDECOMP' ) THEN
          READ(12,*) var,dumb
          WRITE(13,*) 'Not parallel code '
        ELSEIF  ( VAR(1:5) .eq. 'MDRUN' ) THEN
          READ(12,*) var,dumb,MDRUN
        ELSEIF  ( VAR(1:9) .eq. 'PRINTXMOL' ) THEN
          READ(12,*) var,dumb,PXMOL
        ELSEIF  ( VAR(1:8) .eq. 'PRINTCFG' ) THEN
          READ(12,*) var,dumb,PCFG
        ELSEIF  ( VAR(1:5) .eq. 'WSTR' ) THEN
          READ(12,*) var,dumb,NSTR
        ELSEIF  ( VAR(1:7) .eq. 'RETHERM' ) THEN
          READ(12,*) var,dumb,RTHRM
        ELSEIF  ( VAR(1:8) .eq. 'BOXTHERM' ) THEN
          READ(12,*) var,dumb,BXTHRM
        ELSEIF  ( VAR(1:11) .eq. 'RIGIDBOXMIN' ) THEN
          READ(12,*) var,dumb,RBOX(1,:)
        ELSEIF  ( VAR(1:11) .eq. 'RIGIDBOXMAX' ) THEN
          READ(12,*) var,dumb,RBOX(2,:)
        ELSEIF  ( VAR(1:11) .eq. 'THERMBOXMIN' ) THEN
          READ(12,*) var,dumb,TBOX(1,:)
        ELSEIF  ( VAR(1:11) .eq. 'THERMBOXMAX' ) THEN
          READ(12,*) var,dumb,TBOX(2,:)
        ELSEIF  ( VAR(1:11) .eq. 'CUSTOMTHERM' ) THEN
          READ(12,*) var,dumb,CUSTHRM
        ELSEIF ( var(1:6).eq.'MVTIP' ) THEN
         READ(12,*) var,dumb,MVTIP
        ELSEIF ( var(1:7).eq.'MVSTEPS' ) THEN
         READ(12,*) var,dumb,ISTEPS
        ELSEIF ( var(1:11).eq.'TIPVELOCITY' ) THEN
         READ(12,*) var,dumb,tpvel(:)
        ELSEIF ( var(1:4).eq.'CALC' ) THEN
         READ(12,'(A,A1,A70)') var,dumb,ITITLE
        ELSEIF  ( VAR(1:1) .eq. '!' ) THEN
          READ(12,*) 
          WRITE(13,*) 'Comment ',var
        ELSEIF  ( VAR(1:2) .eq. 'C ' ) THEN
          READ(12,*) 
          WRITE(13,*) 'Comment ',var
        ELSE
          WRITE(13,*) var,' unknown'
          READ(12,*) 
        ENDIF
        READ(12,*,IOSTAT=ios) var
        BACKSPACE(12) 
      ENDDO
      CLOSE(12)
      RETURN
      END SUBROUTINE READIN
!
      SUBROUTINE WRITEOPT
!     Write out simulation specifications
      USE SPECIF
      USE STRUCTR
      USE POTS
      USE STRUCTR
      USE MDSTP
      USE PARAMS 
!
      IMPLICIT none
!
      INTEGER :: N,KI,KJ
!
      WRITE(6,*)' Simulation parameters read in.'
      WRITE(6,601)' System    : ',CTITLE
      WRITE(6,601)' Operation : ',ITITLE  
      WRITE(6,*)' Units: '
      WRITE(6,602)'   Length = ', 1.0,' A' 
      WRITE(6,602)'   Mass   = ',1.0,' AMU'
      WRITE(6,602)'   Time   = ',1.0,' fs'
      WRITE(6,602)'   Energy = ',ECONV,' eV'
!
      WRITE(6,607)'  Number of atom types',RTYPES
      WRITE(6,604)'  Temperature:',TEM
      WRITE(6,604)'  Start time:',TTIME
      WRITE(6,604)'  Time Step:',DELTA
      WRITE(6,*)'  Seed :',PSEED

! Potential used 
      WRITE(6,*)' Potential: '
      WRITE(6,*)'  REBO'
      WRITE(6,*)'  Max cutoffs:'
      DO KI=1,RTYPES-1
        DO KJ=KI,RTYPES
          WRITE(6,605)'   ',ATYPE(KI),ATYPE(KJ),SQRT(RMAX(KI,KJ))
        ENDDO
      ENDDO
      WRITE(6,604)'  Nieghbor list buffer: cutoff max +',RLL
      IF (LTORS) WRITE(6,*)'  Torsion'
      IF (ILJ  ) WRITE(6,*)'  Lenard Jones'
      IF (INTG.EQ.1) WRITE(6,*)' Integrator : 3rd order PC'
      IF (INTG.EQ.2) WRITE(6,*)' Integrator : Velocity Verlet'
      IF ( KFLAG.EQ.-1 ) THEN
        WRITE(6,*)'  Thermostat      : Landervan'
      ENDIF
! Options
      WRITE(6,*)' Options   : '
      IF (MDRUN) WRITE(6,*)'  MD simulation'
      IF (WSTRI) WRITE(6,*)'  The initial structure will be writen'
      IF (RTHRM) WRITE(6,*)'  The systems thermostates ',
     &                     ' will be reset to ', NTHRM
      IF(BXTHRM) THEN
        WRITE(6,*)'  Box thermostates will be applied'
        WRITE(6,610)
        WRITE(6,611)' x: ',RMN(1,1),RIG(1,1),THM(1,1)
     &             ,THM(2,1),RIG(2,1),RMN(2,1)
        WRITE(6,611)' y: ',RMN(1,2),RIG(1,2),THM(1,2)
     &             ,THM(2,2),RIG(2,2),RMN(2,2)
        WRITE(6,611)' z: ',RMN(1,3),RIG(1,3),THM(1,3)
     &             ,THM(2,3),RIG(2,3),RMN(2,3)
        WRITE(6,*)'    '
      ENDIF
      IF (HT)  WRITE(6,608)'  System will be heated by',DT,'K/fs to '
     & ,TMAX,' K'
      IF (IRSVDT)WRITE(6,*)'  Variable time step'
      IF (IRSVDT)WRITE(6,604)'    with a max displacement of',DISMAX
      IF (IRSVDT)WRITE(6,*)'  Variable time step'
      IF (IRSVDT)WRITE(6,*)'    with a max displacement of',DISMAX
! Detials
      WRITE(6,*)' Simulation detials'
      WRITE(6,607)'  Total number of atoms',NP
      WRITE(6,612)'  Number of active atoms',nma-nta
     &   ,FLOAT(nma-nta)/FLOAT(NP)
      WRITE(6,612)'  Number of thermostat atoms',nta
     &   ,FLOAT(nta)/FLOAT(NP)
      WRITE(6,612)'  Number of rigid atoms',NP-nma
     &   ,FLOAT(NP-nma)/FLOAT(NP)
      WRITE(6,607)'  Total number of atoms',NP
      WRITE(6,607)'  Number of steps : ', KVC
      WRITE(6,607)'  Number of steps between output writes:',MAXKB
      WRITE(6,607)'  Number of steps between structure writes :',NSTR 
      IF ( PXMOL ) THEN
        WRITE(6,*)'  Positions will be printed out in xmol format'
      ENDIF
      IF ( PCFG ) THEN
        WRITE(6,*)'  Positions will be printed out in cfg format'
        WRITE(6,609)'    Box size for cfg output will be:',cbox(:)
      ENDIF
      WRITE(6,609)'  Max:',RMN(2,:)
      WRITE(6,609)'  Min:',RMN(1,:)
      WRITE(6,*) "Box Density",CDENS,"g/cm^3"
      WRITE(6,*) "Max/Min Density",MXNDENS,"g/cm^3"
      WRITE(6,*)'  Data ouput format for out.dat'
      WRITE(6,2102)
!
      RETURN
 601  FORMAT(2A)
 602  FORMAT(A,F10.3,A) 
 603  FORMAT(A,I3,A)
 604  FORMAT(A,F10.3)
 605  FORMAT(A,2A3,F10.3)
 606  FORMAT(A,A,A,F12.3,A)
 607  FORMAT(A,I12)
 608  FORMAT(A,F10.3,A,F10.3,A)
 609  FORMAT(A,3F10.3)
 610  FORMAT('   |     Rigid    |  Therm   |  Active  |'
     &,'  Therm   |   Rigid     |')
 611  FORMAT(A,F6.2,' <-> ',F6.2,' <-> ',F6.2,' <-> '
     &        ,F6.2,' <-> ',F6.2,' <-> ',F6.2)
 612  FORMAT(A,I8,F10.3)
2102  FORMAT(2x,'STEP',2x,'TIME(fs)',
     &7x,'POTENTIAL',9x,'KINETIC',7x,'TOTAL(eV)',
     &13x,'TEMP (K)',1x,'TSTEP(fs)') 
      END SUBROUTINE WRITEOPT
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Read input structure and set things accordingly 
!
      subroutine read_coord
!
      USE SPECIF
      USE STRUCTR
      USE SPECIF
      USE POTS
      USE PARAMS
      USE MDSTP
!
      IMPLICIT none
!
      INTEGER :: K,KATOM,DM,I,NATOM
      INTEGER :: IOS,N,M
      CHARACTER(6) :: var
      CHARACTER(30) :: dumb
      REAL*8 :: val
C
      INTEGER :: J,IDUM,NRA,NLA
      REAL*8 :: xma,epst,sigt
C
C Read structure 
C
      READ(10,1800) CTITLE 
      READ(10,*) NP
      READ(10,*) TTIME,DELTA
      READ(10,*) (CUBE_r(N),N=1,3)
      READ(10,*) (CUBE(N),N=1,3)
! Put allocations here
      if(np.gt.npmax) then
           write(*,*) 'np= ',np,' greater than npmax= ',npmax
           write(*,*) 'increase npmax and recompile'
           include 'close.inc'
           stop
      endif
!
      DELTSQ=DELTA*DELTA/2.0D0
!
      MTOT = 0.0d0
      RMN(1,:) = (/1.d16,1.0d16,1.0d16/)
      RMN(2,:) = (/-1.0d16,-1.0d16,-1.0d16/)
      DO I=1,NP
        READ(10,*) K,NATOM,(R0(K,N),N=1,3),itr(k)
        do n=1,3
              if (R0(K,N).gt.cube(n)/2.0d0) R0(K,N)=R0(K,N)-cube(n)
              IF ( R0(K,N) .LT. RMN(1,N) ) RMN(1,N)=R0(K,N)
              IF ( R0(K,N) .GT. RMN(2,N) ) RMN(2,N)=R0(K,N)
        enddo
        KTYPE(K)=KT(NATOM)
        MTOT = MTOT + XMASS(KTYPE(K) )
        if(ktype(k).eq.0.OR.ktype(k).GT.NTYPES) then
                write(6,*) 'unknown atom type for atom read in ',k
                include 'close.inc'
                stop
        endif
        noa(ktype(k)) = noa(ktype(k)) + 1
!
!       Reset thermostate value if specified
        IF( RTHRM ) ITR(K) = NTHRM
c
        READ(10,*) K,(R1(K,N),N=1,3)
        READ(10,*) K,(R2(K,N),N=1,3)
        READ(10,*) K,(R3(K,N),N=1,3)
        DO N=1,3
                R1(k,N)=R1(k,N)*DELTA
        ENDDO
        DO N=1,3
                R2(k,N)=R2(k,N)*DELTSQ
        ENDDO
      ENDDO
      REWIND(10)
C
C ESTABLISH INITIAL POSITIONS FOR NEIGHBOR LIST UPDATE
C
      nma = 0   ! number of non fixed atoms (thermostated/active
      nta = 0
      DO I=1,NP
          R0L(I,:)=R0(I,:)
!         Set up list of nonrigid and thermostated atoms
          if((itr(I).eq.0).or.(itr(I).eq.1)) then
            nma = nma + 1
            mlist(nma) = I
            if(itr(I).eq.1) then
              nta = nta + 1
              nlist(nta) = I
            endif
          endif
      ENDDO
C     
      VOL=CUBE(1)*CUBE(2)*CUBE(3)
!     Calculate density
      VOLMXN = (RMN(2,1) - RMN(1,1))*(RMN(2,3) - RMN(1,3))
     &         *(RMN(2,3) - RMN(1,3))
      MXNDENS = MTOT/VOLMXN*DCONV
      CDENS = MTOT/VOL*DCONV
C
      DO 6 I=1,3
           CUBE2(I)=CUBE(I)/2.0D0
6     CONTINUE
C
       TTCONV=2.0d0/3.0d0/FLOAT(NP)
C
C Determine system type
       WRITE(13,*)'Checking atom types to determine which 
     & pibond to use:'
       IF(noa(isulfur).EQ.0.and.noa(ioxygen).EQ.0
     &    .and.noa(iflor).EQ.0) THEN
          SYTP = 1
          WRITE(13,*) ' Using CHF pibond'! IF(mynode.EQ.0)
        ELSEIF ( noa(isulfur).GT.0.OR.noa(ioxygen).GT.0
     &    .OR.noa(iflor).GT.0) THEN
          IF  ( noa(ioxygen).GT.0.and.noa(isulfur).EQ.0
     &    .AND.noa(iflor).EQ.0) THEN
            SYTP = 2
            WRITE(13,*) ' Using CHO pibond'! IF(mynode.EQ.0)
          ELSEIF( noa(ioxygen).EQ.0.and.noa(isulfur).GT.0
     &      .AND.noa(iflor).EQ.0) THEN
            SYTP = 3
            WRITE(13,*) ' Using CHS pibond'! IF(mynode.EQ.0)
          ELSEIF( noa(ioxygen).EQ.0.and.noa(isulfur).EQ.0
     &      .AND.noa(iflor).GT.0) THEN
            SYTP = 4
            WRITE(13,*) ' Using CHF pibond'! IF(mynode.EQ.0)
          ELSE 
            WRITE(13,*) ' System contianes oxygen,sulfur and fluorine'
            WRITE(13,*) ' The S-O,S-F,O-F interactions is not  currently 
     &                 parameterized'
            STOP 
          ENDIF
      ENDIF

!
      return
  100 FORMAT(4I6)
  200 FORMAT(4F12.6)
  300 FORMAT(3E20.11)
  350 FORMAT(2I5,3E20.11,I3)
  360 FORMAT(I5,3E20.11)
 1800 FORMAT(A70)

      END SUBROUTINE read_coord
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Write structure in coord.d format

      subroutine writeinit 
!
      USE STRUCTR
      USE POTS
      USE PARAMS
      USE MDSTP
      USE SPECIF
!
      IMPLICIT none
!
      INTEGER :: I,J,N
      INTEGER :: IDUM,NRA,NLA
!
      IDUM = 1
      NRA = 1
      NLA = 1
!
      WRITE(6,*) ' Writing initial structure'
      WRITE(13,110) CTITLE
      WRITE(13,111) NP,IDUM,NRA,NLA
      WRITE(13,112) TTIME,DELTA
      WRITE(13,113) (CUBE_r(N),N=1,3)
      WRITE(13,113) (CUBE(N),N=1,3)
C
      DO 11 I=1,NP
           WRITE(13,114)I,KT2(KTYPE(I)),(R0(I,N),N=1,3),itr(i),RCLIST(I)
           WRITE(13,115) I,((R1(I,N)/DELTA),N=1,3)
           WRITE(13,115) I,((R2(I,N)/DELTSQ),N=1,3)
           WRITE(13,115) I,((R3(I,N)),N=1,3)
11    CONTINUE
C
      CLOSE(13)
!
      return 
 110  FORMAT(20A2)
 111  FORMAT(4I6)
 112  FORMAT(2E20.11)
 113  FORMAT(3E20.11)
 114  FORMAT(2I5,3E20.11,I3,I8)
 115  FORMAT(I5,3E20.11)
      END SUBROUTINE WRITEINIT
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Write data for current step
!
      subroutine write_data

      USE STRUCTR
      USE POTS
      USE PARAMS
      USE MDSTP
      USE SPECIF

      IMPLICIT none
 
      INTEGER :: I,J,N
      REAL*8 ::  XX,TOTEK,EETOT,KEPA

c
C  CALCULATE KINETIC ENERGY
C
      XX=0.0d0
      DO J=1,3
           DO I=1,NP
                XX=XX+(R1(I,J)**2)*XMASS(KTYPE(I)) 
           enddo
      enddo
      
      TOTEK = XX/(4.0d0*DELTSQ)*ECONV
      EETOT = TOTE + TOTEK
      TEMPK = ENPR*XX/6.0d0/DELTSQ*ECONV
      KEPA = TOTEK/FLOAT(NP)
C
      WRITE(21,2101) LSTEP,TIME,TOTE/FLOAT(NP),
     &   KEPA,EETOT/FLOAT(NP),
     &   TEMPK,DELTA
      call flush(9)
C
C      IF((KFLAG.NE.6).AND.(KFLAG.NE.8)) THEN
C           WRITE(95,*) lstep, TOTE/FLOAT(NP)
C      ENDIF
C
      WRITE(6,601)' PE= ',TOTE,' KE= ',TOTEK,' TOTE= ',EETOT

      IF (QUENCH) THEN
        IF( KEPA.GT.KEMX) KFLAG = 2
        IF( KEPA.LE.KEMX) KFLAG = 1
      ENDIF

c write one-half pair energy for each atom
C        do i=1,np
C             write(95,*) i,ktype(i),eatom(i)
C        enddo
      call flush(85)

      return
  601 FORMAT(3(A,F16.8,1x))
 2101 FORMAT(I10,1x,F16.3,4F16.8,F14.3,F6.3)
      END SUBROUTINE WRITE_DATA
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Write coordinates
!
      subroutine write_coord

      USE STRUCTR
      USE POTS
      USE PARAMS
      USE MDSTP
      USE SPECIF

      IMPLICIT none
 
      INTEGER :: I,J,N
      INTEGER :: IDUM,NRA,NLA

      IDUM = 1
      NRA = 1
      NLA = 1

      WRITE(11,1800) CTITLE
      WRITE(11,100) NP,IDUM,NRA,NLA
      WRITE(11,300) TTIME,DELTA
      WRITE(11,300) (CUBE_r(N),N=1,3)
      WRITE(11,300) (CUBE(N),N=1,3)
C
      DO 11 I=1,NP
        WRITE(11,350) I,KT2(KTYPE(I)),(R0(I,N),N=1,3)
     &                 ,itr(i)
        WRITE(11,360) I,((R1(I,N)/DELTA),N=1,3)
        WRITE(11,360) I,((R2(I,N)/DELTSQ),N=1,3)
        WRITE(11,360) I,((R3(I,N)),N=1,3)
11    CONTINUE
C
      CALL FLUSH(11)
      REWIND(11)
      return 

  100 FORMAT(4I6)
  200 FORMAT(4F12.6)
  300 FORMAT(3E20.11)
  350 FORMAT(I8,I5,3E20.11,I3)
  360 FORMAT(I8,3E20.11)
  500 FORMAT('*CLASSICAL DYNAMICS SIMULATION OF ',20A2)
  600 FORMAT(/,'TOTAL STEPS= ',I6,' DATA WRITTEN EVERY ',I4,
C    &' STEPS WITH ',F7.4,' fs/STEP')
     &' STEPS')
  700 FORMAT(/,'UNITS OF LENGTH, MASS, TIME AND ENERGY:',I2,' A ',I2,
     &' AMU ',I2,' fs ',F8.4,' eV')
  800 FORMAT(/,'NEIGHBOR LIST PARAMETERS: ',F12.3,' A ')
  900 FORMAT(/,'LANGEVIN PARAMETERS: ',F12.3,' k PSEED= ',F9.6,/)
 1200 FORMAT('NEIGHBOR LIST UPDATES: ',/)
 1300 FORMAT(10I4,/)
 1400 FORMAT(8X,'ENERGY(eV)',5X,'T',6X,'TSTEP(fs)  TIME(fs)',/)
 1500 FORMAT(I6,F14.8,F10.3,2F9.3)
 1800 FORMAT(20A2)

      END SUBROUTINE WRITE_COORD
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
