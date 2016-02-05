C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C     Version 1.0 01/16/2010 T. W. Kemper                      C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C  Subroutines involved in reading and writing data
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  Write header of output file
!
      SUBROUTINE WRTHEAD
!
      USE MPIvars 
      USE SPECIF
      USE TM
!
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
      INTEGER :: sec(3),day1,day2,day3
!
      IF ( mynode.EQ.0 ) THEN
        WRITE(6,*) ''
     &,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     &,'+                  REBO  (CHF/CHO/CHS)                      +'
     &,'+ THIS PROGRAMS PERFORM MOLECULAR DYNAMICS SIMULATIONS WITH +'
     &,'+ SHORT RANGE BOND-ORDER POTENTIALS AND LONG RANGE          +'
     &,'+ LENNARD-J0NES POTENTIALS FOR HYDROGEN, CARBON AND OXYGEN  +'
     &,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     &,'+  Units: mass = AMU s, length = Angstroms, time = fs       +'
     &,'+         energy = eV s                                     +'
     &,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     &,'+                                                           +'
     &,'+  Details about development and properties of potential    +'
     &,'+  is described at :                                        +'
     &,'+  1. CH                                                    +'
     &,'+     D.W. Brenner+, O.A. Shenderova,                       +'
     &,'+     J.A. Harrison, S.J. Stewart, B. Ni, S.B. Sinnott,     +'
     &,'+     Journal of Physics: Condensed Matter 14,              +'
     &,'+     783-802 (2002).                                       +'
     &,'+  2. CHO                                                   +'
     &,'+     Boris Ni, Ki-Ho Lee and Susan B. Sinnott              +'
     &,'+     J. Phys.: Condens. Matter 16 (2004) 7261-7275         +'
     &,'+  3. CHS                                                   +'
     &,'+     Travis Kemper                                         +'
     &,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     &,'+    Version Parellel 10.01                                  +' 
     &,'+       05/14/2011 T. W. Kemper: traviskemper@ufl.ed        +'
     &,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

C Calculate current date 
!        call idate(day1,day2,day3)   ! today(1)=day, (2)=month, (3)=year
!        call itime(sec)     ! now(1)=hour, (2)=minute, (3)=second
C        WRITE(6,102) day1,day2,2000+day3,sec
        WRITE(6,*)
     &,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        WRITE(6,*)''
      ENDIF
!
 101  FORMAT(a16,3(i3,a8,1x),f6.1,1x,a8)
 102  FORMAT( ' + Date: ', i2.2, '/', i2.2, '/', i4.4, '; time ',
     & i2.2, ':', i2.2, ':', i2.2, '                           +' )
!
      END SUBROUTINE WRTHEAD
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Read in molecular structure from mol.xyz
!
      SUBROUTINE READMOL
!
      USE MPIvars
      USE BEAM
      USE ANALYSIS !only 4 DEN
      USE POTS
!
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
      INTEGER :: I,N
!
      IF( mynode.EQ.0 ) THEN
        OPEN(UNIT=20,file='mol.xyz',status='old')
        READ(20,*) NPM,NC,EB(DEPD),TBUF  !# atoms in mol.,# of cluster/beam
        READ(20,*) MTITLE
!     Allocate need arrays
        ALLOCATE(
     &  ITYPEM(NPM)
     & ,RM(NPM,3)
     & )
        DO I=1,NPM 
          READ(20,*) N,RM(i,1),RM(i,2),RM(i,3)
          ITYPEM(i) = N
          MOLMASS = MOLMASS + XMASS(KT(N))
        ENDDO
        CLOSE(20)
      ENDIF
!
      RETURN
      END SUBROUTINE READMOL
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Read in settings for simulation from input file
!  input file must specify variables with the following format:
!   [VAR] = [VALUE]
!  and is blank line terminated
!
      SUBROUTINE READIN 
!
      USE MPIvars
      USE SPECIF
      USE STRUCTR
      USE POTS
      USE ANALYSIS
      USE STRUCTR
      USE BEAM
!   
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
      INTEGER :: IOS,N,J,ierr
      CHARACTER(20) :: var
      CHARACTER(2) :: dumb
      REAL*8 :: val
!
! Read in input file
      IF ( mynode.EQ.0) THEN 
       WRITE(13,*)'Reading input.d:'
       READ(12,*,IOSTAT=ios) var
       BACKSPACE(12) 
       DO WHILE( IOS.EQ.0)
C       BACKSPACE 12  !Find variable then go back a space and read it in
        IF ( var(1:4).eq.'CALC' ) THEN
         READ(12,'(A,A1,A70)') var,dumb,ITITLE
        ELSEIF ( var(1:6).eq.'SYSTEM' ) THEN
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
        ELSEIF  ( VAR(1:8) .eq. 'NABORCUT' ) THEN
          READ(12,*) var,dumb,RLL
        ELSEIF  ( VAR(1:13) .eq. 'SPACIALDECOMP' ) THEN
          READ(12,*) var,dumb,NDP(:)
        ELSEIF  ( VAR(1:7) .eq. 'PRINTPR' ) THEN
          READ(12,*) var,dumb,PRNPR
        ELSEIF  ( VAR(1:6) .eq. 'PRINTF' ) THEN
          READ(12,*) var,dumb,PRNF
        ELSEIF  ( VAR(1:5) .eq. 'PAIR1' ) THEN
          READ(12,*) var,dumb,PR1
        ELSEIF  ( VAR(1:5) .eq. 'PAIR2' ) THEN
          READ(12,*) var,dumb,PR2
        ELSEIF  ( VAR(1:10) .eq. 'INTEGRATOR' ) THEN
          READ(12,*) var,dumb,INTG
        ELSEIF  ( VAR(1:4) .eq. 'HEAT' ) THEN
          READ(12,*) var,dumb,HT
        ELSEIF  ( VAR(1:6) .eq. 'DELTAT' ) THEN
          READ(12,*) var,dumb,DT
        ELSEIF  ( VAR(1:11) .eq. 'MOLANALYSIS' ) THEN
          READ(12,*) var,dumb,MOLAN
        ELSEIF  ( VAR(1:11) .eq. 'DEPANALYSIS' ) THEN
          READ(12,*) var,dumb,DEPAN
        ELSEIF  ( VAR(1:12) .eq. 'DEPTHPROFILE' ) THEN
          READ(12,*) var,dumb,DEPPROF
        ELSEIF  ( VAR(1:12) .eq. 'DEPTHPDIV' ) THEN
          READ(12,*) var,dumb,DPDIV
        ELSEIF  ( VAR(1:13) .eq. 'REACTANALYSIS' ) THEN
          READ(12,*) var,dumb,RECAN
        ELSEIF  ( VAR(1:15) .eq. 'SURFACEANALYSIS' ) THEN
          READ(12,*) var,dumb,SMOLAN
        ELSEIF  ( VAR(1:12) .eq. 'PRINTSUBINFO' ) THEN
          READ(12,*) var,dumb,RECSUB
        ELSEIF  ( VAR(1:7) .eq. 'READREF' ) THEN
          READ(12,*) var,dumb,REFSTR
        ELSEIF  ( VAR(1:6) .eq. 'DEPDIM' ) THEN
          READ(12,*) var,dumb,DEPD
        ELSEIF  ( VAR(1:6) .eq. 'MKBEAM' ) THEN
          READ(12,*) var,dumb,MKBM
        ELSEIF  ( VAR(1:6) .eq. 'CENTER' ) THEN
          READ(12,*) var,dumb,CENTR
        ELSEIF  ( VAR(1:5) .eq. 'SHIFT' ) THEN
          READ(12,*) var,dumb,SHFT
        ELSEIF  ( VAR(1:8) .eq. 'CLEARGAS' ) THEN
          READ(12,*) var,dumb,CLRGAS
        ELSEIF  ( VAR(1:5) .eq. 'MDRUN' ) THEN
          READ(12,*) var,dumb,MDRUN
        ELSEIF  ( VAR(1:9) .eq. 'PRINTXMOL' ) THEN
          READ(12,*) var,dumb,PXMOL
        ELSEIF  ( VAR(1:5) .eq. 'WSTR' ) THEN
          READ(12,*) var,dumb,NSTR
        ELSEIF  ( VAR(1:8) .eq. 'PRINTCFG' ) THEN
          READ(12,*) var,dumb,PCFG
        ELSEIF  ( VAR(1:6) .eq. 'SUBBUF' ) THEN
          READ(12,*) var,dumb,HBUF
        ELSEIF  ( VAR(1:6) .eq. 'DEPBUF' ) THEN
          READ(12,*) var,dumb,DBUF
        ELSEIF  ( VAR(1:8) .eq. 'THERMBUF' ) THEN
          READ(12,*) var,dumb,THBUF
        ELSEIF  ( VAR(1:6) .eq. 'NTHERM' ) THEN
          READ(12,*) var,dumb,NTHRM
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
        ELSEIF  ( VAR(1:12) .eq. 'WRITEINITIAL' ) THEN
          READ(12,*) var,dumb,WSTRI
        ELSEIF ( var(1:6).eq.'MVTIP' ) THEN
         READ(12,*) var,dumb,MVTIP
        ELSEIF ( var(1:7).eq.'MVSTEPS' ) THEN
         READ(12,*) var,dumb,ISTEPS
        ELSEIF ( var(1:11).eq.'TIPVELOCITY' ) THEN
         READ(12,*) var,dumb,tpvel(:)
        ELSEIF  ( VAR(1:8) .eq. 'ADDFORCE' ) THEN
          READ(12,*) var,dumb,ADDFC
        ELSEIF  ( VAR(1:9) .eq. 'FORCESTEP' ) THEN
          READ(12,*) var,dumb,FSTEP
        ELSEIF  ( VAR(1:9) .eq. 'FPRESSURE' ) THEN
          READ(12,*) var,dumb,PRADF
        ELSEIF  ( VAR(1:13) .eq. 'SURFDIMENSION' ) THEN
          READ(12,*) var,dumb,SDIR
        ELSEIF  ( VAR(1:13) .eq. 'PRESSURETHERM' ) THEN
          READ(12,*) var,dumb,PTHRM
        ELSEIF  ( VAR(1:11) .eq. 'PRESSUREBUF' ) THEN
          READ(12,*) var,dumb,PRBUF
        ELSEIF  ( VAR(1:11) .eq. 'PRINTLOAD' ) THEN
          READ(12,*) var,dumb,PRNLD
        ELSEIF  ( VAR(1:17) .eq. 'CROSSLINKANALYSIS' ) THEN
          READ(12,*) var,dumb,CRSAN
        ELSEIF  ( VAR(1:1) .eq. '!' ) THEN
          READ(12,*) 
        ELSEIF  ( VAR(1:2) .eq. 'C ' ) THEN
          READ(12,*) 
        ELSE
          WRITE(13,*) var,' unknown'
          READ(12,*) 
        ENDIF
        READ(12,*,IOSTAT=ios) var
        BACKSPACE(12) 
       ENDDO
       CLOSE(12)
!      Open optional files 
        IF (mynode.EQ.0) THEN
         IF(PXMOL) open(1,file='out.xmol',status='unknown')
         IF(MOLAN) OPEN(UNIT=28,FILE='out.gasmol',STATUS='UNKNOWN')
         IF(SMOLAN) OPEN(UNIT=29,FILE='out.surfmol',STATUS='UNKNOWN')
         IF(PRNLD) OPEN(UNIT=51,FILE='out.load',STATUS='UNKNOWN')
         IF(RECAN)OPEN(UNIT=31,FILE='out.element',STATUS='UNKNOWN')
         IF(CRSAN)OPEN(UNIT=32,FILE='out.crslink',STATUS='UNKNOWN')
         IF(RECAN)OPEN(UNIT=33,FILE='out.bonds',STATUS='UNKNOWN')
        ENDIF
       ENDIF
       CALL MPI_BCAST(ITITLE,1,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(CTITLE,1,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(ILJ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(KVC,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(TMAX,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(TEM,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(MAXKB,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(KFLAG,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(PSEED,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(RLL,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(PRNPR,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(PRNF,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(PR1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(PR2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(INTG,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(DT,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(HT,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(MOLAN,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(DEPAN,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(SMOLAN,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(RECAN,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(RECSUB,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(REFSTR,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(DEPD,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(MKBM,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(CENTR,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(CLRGAS,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(PXMOL,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(NSTR,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(PCFG,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(HBUF,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(DBUF,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(THBUF,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(NTHRM,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(RTHRM,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(WSTRI,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(F02,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(F12,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(F32,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(NDP,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(IPOT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(MVTIP,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(ISTEPS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(disp,3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(tpvel,3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(PRNLD,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
!      CALL MPI_BCAST(NLJ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
!
!
      RETURN
      END SUBROUTINE READIN
!
      SUBROUTINE WRITEOPT
!     Write out simulation specifications
      USE MPIvars
      USE SPECIF
      USE STRUCTR
      USE POTS
      USE ANALYSIS
      USE STRUCTR
      USE BEAM
!
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
      INTEGER :: N,KI,KJ
!
      IF( mynode.EQ.0) THEN
!
      WRITE(6,*)' Simulation parameters read in.'
      WRITE(6,601)' System    : ',CTITLE
      WRITE(6,601)' Operation : ',ITITLE  
      WRITE(6,*)' Units: '
      WRITE(6,602)'   Length = ', 1.0,' A' 
      WRITE(6,602)'   Mass   = ',1.0,' AMU'
      WRITE(6,602)'   Time   = ',1.0,' fs'
      WRITE(6,602)'   Energy = ',ECONV,' eV'
      WRITE(6,*)' Spacial decomposition',NDP(:)
!
      WRITE(6,607)'  Number of atom types',RTYPES
      WRITE(6,604)'  Temperature:',TEM
      WRITE(6,604)'  Start time:',TTIME 
      WRITE(6,604)'  Time Step:',DELTA
      WRITE(6,*)'  Seed :',PSEED
! Potential used 
      WRITE(6,*)' Potential: '
      WRITE(6,*)'  REBO'
      WRITE(6,*)'    Max cutoffs:'
      DO KI=1,RTYPES-1
        DO KJ=KI,RTYPES
          WRITE(6,605)'    ',ATYPE(KI),ATYPE(KJ),SQRT(RMAX(KI,KJ))
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
      IF (ENMIN) WRITE(6,*)'  Energy minimization'
      IF (MDRUN) WRITE(6,*)'  MD simulation'
      IF (CENTR) WRITE(6,*)'  The system will be centered'
      IF (WSTRI) WRITE(6,*)'  The initial structure will be writen'
      IF (RTHRM) WRITE(6,*)'  The systems thermostates ',
     &                     ' will be reset to ', NTHRM
      IF (MKBM)  THEN
        WRITE(6,606)'  A beam of ',MTITLE,' at ',EB(DEPD)
     &           ,'eV will be added '
        WRITE(6,603)'    In the ',DEPD,' indexed direction'
        WRITE(6,604)'    Starting at :',BMHT
        WRITE(6,602)'    Which should be ',DBUF,' A'
        WRITE(6,602)'    above the max point of the input geom '
     &           ,RMN(2,DEPD),' A'
      ENDIF
      IF (PRNPR) THEN
        WRITE(6,*)'  Information about',
     &                           PR1,PR2
     &                      ,'will be printed'
      ENDIF
      IF (PRNF) WRITE(6,*)'  Forces for each atom will printed'
      IF (HT)  WRITE(6,608)'  System will be heated by',DT,'K/fs to '
     & ,TMAX,' K'
      IF (IRSVDT)WRITE(6,*)'  Variable time step'
      IF (IRSVDT)WRITE(6,604)'    with a max displacement of',DISMAX
      WRITE(6,*)'  Analysis:'
      IF(MOLAN)WRITE(6,*)'    Molecular composition'
      IF(SMOLAN)WRITE(6,*)'    Surface species analysis'
      IF(RECAN)WRITE(6,*)'    Reactions
     & starting from previous reaction ',REC
      IF(DEPAN)WRITE(6,*)'    Deposition'
      IF(CLRGAS)WRITE(6,*)
     &  '    Gas molecules will be placed under the substrate'
      IF(ADDFC ) THEN
        WRITE(6,*)'    A presure of ',PRADF,'MP will be added'
        WRITE(6,*)'     Conversion from MPa to AMU/Afs2',FCONV
        WRITE(6,*)'     Surface area ',SAREA,' A2'
        WRITE(6,*)'     Top atoms ',TATOMS,' will have a force of',ADFT
      WRITE(6,*)'     Bottom atoms ',BATOMS,' will have a force of',ADFB
      ENDIF 
      IF ( MVTIP ) THEN
        WRITE(6,*)'     ITR = 3 atoms will be moved'
        WRITE(6,*)'       at velocity',tpvel,' m/s'
        WRITE(6,*)'       which was found to be ',disp(:),' A'
        WRITE(6,*)'       every ',isteps,' steps'
      ENDIF 
      IF(DEPPROF ) THEN
        WRITE(6,*) '    Depth profile will be tracked'
        WRITE(6,*) '     The maximum of the substrate will be',SUBCT
        WRITE(6,*) '     Substrate will be divided by',SDIV
     &             ,'divisions of ',DPDIV,' A each'
      ENDIF 
! Detials
      WRITE(6,*)' Simulation detials'
      WRITE(6,607)'  Total number of atoms',NP
      WRITE(6,612)'  Number of active atoms',nma-nta
     &   ,FLOAT(nma-nta)/FLOAT(NP)
      WRITE(6,612)'  Number of thermostat atoms',nta
     &   ,FLOAT(nta)/FLOAT(NP)
      WRITE(6,612)'  Number of rigid atoms',NP-nma
     &   ,FLOAT(NP-nma)/FLOAT(NP)
      WRITE(6,607)'  Number of steps : ', KVC
      WRITE(6,607)'  Number of steps between energy output writes:'
     & ,MAXKB
      WRITE(6,607)'  Number of steps between structure writes :',NSTR
      IF ( PXMOL ) THEN
        WRITE(6,*)'  Positions will be printed out in xmol format'
      ENDIF
      IF ( PCFG ) THEN
        WRITE(6,*)'  Positions will be printed out in cfg format'
        WRITE(6,609)'    Box size for cfg output will be:',CUBE(:)
      ENDIF
      WRITE(6,609)'  Max:',RMN(2,:)
      WRITE(6,609)'  Min:',RMN(1,:)
      WRITE(6,*) "Box Density",CDENS,"g/cm^3"
      WRITE(6,*) "Max/Min Density",MXNDENS,"g/cm^3"
      WRITE(6,*)'  Data ouput format for out.dat'
      WRITE(6,2102)
!
      ENDIF
      RETURN
 601  FORMAT(2A)
 602  FORMAT(A,F10.3,A) 
 603  FORMAT(A,I3,A)
 604  FORMAT(A,F16.3)
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
      USE MPIvars
      USE dyn_array
      USE STRUCTR
      USE SPECIF
      USE PARAMS
      USE POTS
      USE ANALYSIS
      USE BEAM
      USE MDSTP
!
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
C Local variables
      INTEGER :: N,I,II,NN,NATOM,K
     &,ierr,M
      INTEGER n_pass,nma_node,NPtot
C
C RESTART INFORMATION FROM PREVIOUS RUN
C
!     Read in structure information 
      IF (mynode.EQ.0) THEN
         READ(10,1800) CTITLE
         READ(10,*) NP
         READ(10,*) TTIME,DELTA
         READ(10,*) (CUBE_r(N),N=1,3)
         READ(10,*) (CUBE(N),N=1,3)
         IF( MKBM ) THEN
          CALL READMOL
          WRITE(*,*) 'Molecule read in number of atoms changed from'
     &    ,NP
          NP = NP + NPM*NC
          WRITE(*,*) ' To ',NP
        ENDIF
      ENDIF
      CALL MPI_BCAST(NP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(TTIME,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(DELTA,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(CUBE_r,3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(CUBE,3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
!
!     R
      DELTSQ=DELTA*DELTA/2.0D0
      CALL arrange_cell
      CALL allocate_array
      IF( MKBM ) NP = NP - NPM*NC
!
!     Read in initial possitions
      IF ( mynode.EQ.0) THEN
        MTOT = 0.0d0 
        RMN(1,:) = (/1.d16,1.0d16,1.0d16/)
        RMN(2,:) = (/-1.0d16,-1.0d16,-1.0d16/)
        DO I=1,NP
          IF( RECAN ) THEN
            READ(10,*)K,NATOM,R0i(:,K),ITRi(k),IMOLo(K),RCLIST(K)
          ELSE
            READ(10,*) K,NATOM,(R0i(N,K),N=1,3),ITRi(K)
            RCLIST(K) = 0
            IMOLo(K) = 0
          ENDIF
          DO N=1,3
            R0i(N,K)=R0i(N,K)-CUBE(N)*ANINT(R0i(N,K)/CUBE(N))
!           Calculate r_min and r_max
            IF ( R0i(N,K) .GT. RMN(2,N) ) RMN(2,N)=R0i(N,K)
            IF ( R0i(N,K) .LT. RMN(1,N) ) RMN(1,N)=R0i(N,K)
          ENDDO
          KTYPEi(K)=KT(NATOM)
          MTOT = MTOT + XMASS(KTYPEi(K) )
          READ(10,*) K,(R1i(N,K),N=1,3)
          DO N=1,3
            R1i(N,K)=R1i(N,K)*DELTA
          ENDDO
          READ(10,*) K,(R2i(N,K),N=1,3)
          DO N=1,3
            R2i(N,K)=R2i(N,K)*DELTSQ
          ENDDO
          READ(10,*) K,(R3i(N,K),N=1,3)
!         Reset thermostate value if specified
          IF( RTHRM ) ITRi(K) = NTHRM
          IF (KTYPEi(K).EQ.0) THEN
            WRITE(13,*) 'UNKNOWN atom type for atom ',K
            WRITE(13,*) K,NATOM,R0(:,K),ITR(K),IMOLo(K),RCLIST(K)
            CALL write_error
          ENDIF
        ENDDO
        IF( BXTHRM ) CALL SETTHRMBOX
        IF( CUSTHRM ) CALL SETTHRMCUST
        IF( PTHRM ) CALL prestherm
! Center structure 
        IF ( CENTR ) CALL centstr
! I shift structure 
        IF ( SHFT) CALL shift 
! Make super cell
!        IF ( SPCELL ) CALL supcell
! Add beam molecules
        IF( MKBM ) CALL MKBEAM
      ENDIF
!
!     Pass coordinate information to nodes
      nma = 0
      nta = 0
      nlist=0      
      NN=0
      nta=0
      n_pass=INT(NP/NSendMAX)
      np_remain=0
      np_remain(n_pass)=MOD(NP,NSendMAX)
      DO i=1,n_pass
        np_send(i)=NSendMAX+np_remain(i)
C        WRITE(*,*)mynode,i,np_send(i)
      ENDDO
      IF (mynode.EQ.0) THEN 
        K=0
        TATOMS = 0 
        BATOMS = 0
      ENDIF 
      DO i=1,n_pass
        IF (mynode.EQ.0) THEN
          DO ii=1,np_send(i)
            K = K + 1
            KA(ii) = K
            KTYPE(ii)= KTYPEi(K)
            ITR(ii) = ITRi(k)
            R0(:,ii) = R0i(:,k)
            R1(:,ii) = R1i(:,k)
            R2(:,ii) = R2i(:,k)
            R3(:,ii) = R3i(:,k)
            IF (itr(ii).NE.2) THEN
              nma = nma + 1
              IF (itr(ii).EQ.1) THEN
                nta=nta+1
                nlist(KA(ii))=nta
              ENDIF
            ENDIF
            IF ( ADDFC ) THEN
              IF( ITR(ii).EQ.3 ) TATOMS = TATOMS + 1
              IF( ITR(ii).EQ.4 ) BATOMS = BATOMS + 1
            ENDIF 
          ENDDO
        ENDIF
!
        CALL MPI_BCAST(KA,np_send(i),MPI_INTEGER,0,MPI_COMM_WORLD,
     &       ierr)
        CALL MPI_BCAST(KTYPE,np_send(i),MPI_INTEGER,0,MPI_COMM_WORLD,
     &       ierr)
        CALL MPI_BCAST(itr,np_send(i),MPI_INTEGER,0,MPI_COMM_WORLD,
     &       ierr)
        CALL MPI_BCAST(R0,3*np_send(i),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(R1,3*np_send(i),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(R2,3*np_send(i),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(R3,3*np_send(i),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
!
        CALL check_center(NN,i)
      ENDDO
!     Make sure all atoms have been distributed
      IF ( mynode.EQ.0) THEN
        IF (K.NE.NP) THEN
          WRITE(13,*)'error in distributing atoms',
     &  'total distributed', K,'not equal to # read in',NP
          CALL write_error
        ENDIF
      ENDIF
!
      CALL MPI_BCAST(nta,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nlist,NP,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      NP_node=NN
      CALL MPI_REDUCE(NP_node,NPtot,1,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(NPtot,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      IF (NPtot.NE.NP) THEN
        WRITE(*,*)'error in distributing atom',NPtot,NP,ierr
        CALL write_error
      ENDIF
      nma_node=0
      noa_node=0
      DO i=1,NP_node
        noa_node(ktype_node(i)) = noa_node(ktype_node(i)) + 1
        IF(itr_node(i).LT.2) THEN
          nma_node = nma_node + 1
        endif
      ENDDO
      CALL MPI_REDUCE(noa_node,noa,NTYPES,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(noa,NTYPES,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(nma_node,nma,1,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nma,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
C
      VOL=CUBE(1)*CUBE(2)*CUBE(3)
      VOLMXN = (RMN(2,1) - RMN(1,1))*(RMN(2,3) - RMN(1,3))
     &         *(RMN(2,3) - RMN(1,3))
      MXNDENS = MTOT/VOLMXN*DCONV
      CDENS = MTOT/VOL*DCONV
      IF ( ADDFC ) THEN
!       Calculate surface area
        SAREA = 1.0d0 
        DO M=1,3
          IF( M.NE.SDIR) SAREA= SAREA*CUBE(M)
        ENDDO
        ADFT = PRADF*FCONV/DELTSQ/FLOAT(TATOMS)*SAREA
        ADFB = -PRADF*FCONV/DELTSQ/FLOAT(BATOMS)*SAREA
      ENDIF
!     Calculate initial displacement
      IF ( MVTIP ) THEN 
         DO I = 1,3
          disp(I) = tpvel(I) * float(ISTEPS) * VCONV
         ENDDO
       CALL MPI_BCAST(disp,3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(tpvel,3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      ENDIF 
!
      DO 6 I=1,3
           CUBE2(I)=CUBE(I)/2.0D0
6     CONTINUE 
C
      TTCONV=2.0d0/3.0d0/FLOAT(NP)
!
      IF (mynode.EQ.0) THEN
C Determine system type
       WRITE(13,*)'Checking atom types to determine which 
     & pibond to use:'
        IF(noa(isulfur).EQ.0.and.noa(ioxygen).EQ.0
     &  .and.noa(iflor).EQ.0) THEN
          SYTP = 1
          WRITE(13,*) ' Using CHF pibond'! IF(mynode.EQ.0)
        ELSEIF ( noa(isulfur).GT.0.OR.noa(ioxygen).GT.0
     &  .OR.noa(iflor).GT.0) THEN
          IF  ( noa(ioxygen).GT.0.and.noa(isulfur).EQ.0
     &    .AND.noa(iflor).EQ.0) THEN
            SYTP = 2
            WRITE(13,*) ' Using CHO pibond'! IF(mynode.EQ.0)
          ELSEIF( noa(ioxygen).EQ.0.and.noa(isulfur).GT.0
     &    .AND.noa(iflor).EQ.0) THEN
            SYTP = 3
            WRITE(13,*) ' Using CHS pibond'! IF(mynode.EQ.0)
          ELSEIF( noa(ioxygen).EQ.0.and.noa(isulfur).EQ.0
     &    .AND.noa(iflor).GT.0) THEN
            SYTP = 1
            WRITE(13,*) ' Using CHF pibond'! IF(mynode.EQ.0)
          ELSE 
            WRITE(13,*) ' System contianes oxygen,sulfur and fluorine'
            WRITE(13,*) ' The S-O,S-F,O-F interactions is not  currently 
     &                 parameterized'
            WRITE(13,*) ' Or System contianes flourine which 
     &                is not included yet'
            STOP 
          ENDIF
        ENDIF
      ENDIF
!
      CALL MPI_BCAST(SYTP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!
      return
  100 FORMAT(4I6)
  200 FORMAT(4F12.6)
  300 FORMAT(3E20.11)
  350 FORMAT(2I5,3E20.11,I3)
  360 FORMAT(I5,3E20.11)
 1800 FORMAT(20A2)
      END SUBROUTINE READ_COORD
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Write data for current step
      subroutine write_data
!
      USE dyn_array
      USE MPIvars
      USE STRUCTR
      USE POTS
      USE PARAMS
      USE MDSTP
      USE SPECIF

      IMPLICIT none
 
      INTEGER :: I,J,N
     &,ierr
      REAL*8 ::  XX,TOTEK,EETOT,XXtmp,TEMPK
!
      INCLUDE 'mpif.h'
c
C  CALCULATE KINETIC ENERGY
C
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

      IF (mynode.EQ.0) THEN
C        WRITE(*,*) 'XX',XX
        TOTEK = XX/(4.0d0*DELTSQ)*ECONV
        EETOT = TTOTE + TOTEK
        TEMPK = ENPR*XX/6.0d0/DELTSQ*ECONV
C
        WRITE(21,2101) LSTEP,TIME,TTOTE/FLOAT(NP),TOTEK/FLOAT(NP),
     &   EETOT/FLOAT(NP),TEMPK,DELTA
        call flush(9)
        IF((KFLAG.NE.6).AND.(KFLAG.NE.8)) THEN
           WRITE(6,601) 'PE= ',TTOTE,' KE= ',TOTEK
     &            ,'  TOTE= ',EETOT
        ENDIF
      ENDIF
C      WRITE(*,*) 'finish write 2'
      return
  601 FORMAT(3(A,F16.8,1x))
 2101 FORMAT(I10,1x,F16.3,4F16.8,F14.3,F6.3)
      end 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Reduce the coordinates on all node to head node r*i
!
      SUBROUTINE reduce_coord
!      
      USE dyn_array
      USE MPIvars
      USE STRUCTR
      USE POTS
      USE PARAMS
      USE MDSTP
      USE SPECIF
      USE ANALYSIS
      USE beam
! 
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
C local variable
      INTEGER  :: i,ii,NN,KK,N,K,ierr,TP
     & ,NP_xmol,n_source
     & ,LP,D,J
      INTEGER istat(MPI_STATUS_SIZE)
      REAL*8  r0_xmol(3),r1_xmol(3),r2_xmol(3),r3_xmol(3)   
!
!     Save coordinates to buffer
      NN=0
      DO i=1,NP_node
        NN=NN+1
        buf_cor(NN)=float(NA(i))
        buf_cor(NN+1)=float(ktype_node(i))
        buf_cor(NN+2)=float(itr_node(i))
        buf_cor(NN+3)=r0_node(1,i)
        buf_cor(NN+4)=r0_node(2,i)
        buf_cor(NN+5)=r0_node(3,i)
        buf_cor(NN+6)=r1_node(1,i)
        buf_cor(NN+7)=r1_node(2,i)
        buf_cor(NN+8)=r1_node(3,i)
        buf_cor(NN+9)=r2_node(1,i)
        buf_cor(NN+10)=r2_node(2,i)
        buf_cor(NN+11)=r2_node(3,i)
        buf_cor(NN+12)=r3_node(1,i)
        buf_cor(NN+13)=r3_node(2,i)
        buf_cor(NN+14)=r3_node(3,i)
        NN=NN+14
      ENDDO
!     Save coordinates to local value on node 0
      IF (mynode.EQ.0) THEN
        DO I=1,NP_node
               K = NA(I)
               KTYPEi(K)= ktype_node(I)
               ITRi(K)  = itr_node(i)
               R0i(1,K) = R0_node(1,I)
               R0i(2,K) = R0_node(2,I)
               R0i(3,K) = R0_node(3,I)
               R1i(:,K) = R1_node(:,I)
               R2i(:,K) = R2_node(:,I)
               R3i(:,K) = R3_node(:,I)
        ENDDO
        DO n_source=1,nprocs-1
          CALL MPI_RECV(NP_xmol,1,MPI_INTEGER,
     &    n_source,1,MPI_COMM_WORLD,istat,ierr)
          CALL MPI_RECV(buf_cor_recv,15*NP_xmol,MPI_REAL8,
     &    n_source,2,MPI_COMM_WORLD,istat,ierr)
          NN=0
          DO ii=1,NP_xmol
            NN=NN+1
            K=INT(buf_cor_recv(NN))
            KK=INT(buf_cor_recv(NN+1))
            TP=INT(buf_cor_recv(NN+2))
            r0_xmol(1)=buf_cor_recv(NN+3)
            r0_xmol(2)=buf_cor_recv(NN+4)
            r0_xmol(3)=buf_cor_recv(NN+5)
            r1_xmol(1)=buf_cor_recv(NN+6)
            r1_xmol(2)=buf_cor_recv(NN+7)
            r1_xmol(3)=buf_cor_recv(NN+8)
            r2_xmol(1)=buf_cor_recv(NN+9)
            r2_xmol(2)=buf_cor_recv(NN+10)
            r2_xmol(3)=buf_cor_recv(NN+11)
            r3_xmol(1)=buf_cor_recv(NN+12)
            r3_xmol(2)=buf_cor_recv(NN+13)
            r3_xmol(3)=buf_cor_recv(NN+14)
            NN=NN+14
            KTYPEi(K)=KK
            ITRi(k) = TP
            R0i(1,K) = R0_xmol(1)
            R0i(2,K) = R0_xmol(2)
            R0i(3,K) = R0_xmol(3)
            R1i(:,K) = R1_xmol(:)
            R2i(:,K) = R2_xmol(:)
            R3i(:,K) = R3_xmol(:)
          ENDDO
        ENDDO
!     Pass coordinates to node 0
      ELSE
        CALL MPI_SEND(NP_node,1,MPI_INTEGER,0
     &  ,1,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(buf_cor,15*NP_node,MPI_REAL8,0
     &  ,2,MPI_COMM_WORLD,ierr)
      ENDIF
      end 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Write coordinates
!
      SUBROUTINE write_coord
!      
      USE dyn_array
      USE MPIvars
      USE STRUCTR
      USE POTS
      USE PARAMS
      USE MDSTP
      USE SPECIF
      USE ANALYSIS
      USE beam
! 
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
C local variable
      INTEGER  :: i,ii,NN,KK,N,K,ierr,TP
     & ,NP_xmol,n_source
     & ,LP,D,J
      REAL*8, DIMENSION(3) :: R0F,R1V,R2A
      REAL*8 :: RNPSQ,PRNRNP,VSQ,ASQ,XX,TEMPKi
      CHARACTER *79 FILENAME
!
!     Print out coord.d file 
      IF ( mynode.EQ.0 ) THEN
        WRITE(11,1800) CTITLE
        WRITE(11,100) NP
        WRITE(11,300) TTIME,DELTA
        WRITE(11,300) (CUBE_r(N),N=1,3)
        WRITE(11,300) (CUBE(N),N=1,3)
        DO I=1,NP
           WRITE(11,350) I,KT2(KTYPEi(I)),
     &     (R0i(N,I),N=1,3),ITRi(i),IMOLi(I),RCLIST(I)
           WRITE(11,360) I,((R1i(N,I)/DELTA),N=1,3)
           WRITE(11,360) I,((R2i(N,I)/DELTSQ),N=1,3)
           WRITE(11,360) I,(R3i(N,I),N=1,3)
        ENDDO
        REWIND(11)
!       Print xmol output
        IF ( PXMOL ) THEN
          WRITE(1,100) NP
          WRITE(1,1800) CTITLE
          DO I=1,NP
           WRITE(1,110) KT2(KTYPEi(I)),(R0i(N,I),N=1,3),RCLIST(I)
          ENDDO     
        ENDIF 
!       Print cfg files for atomeye 
        IF ( PCFG ) THEN
           N = LSTEP / NSTR + 1
           LP = 19
!          Write output
           WRITE (FILENAME, '("Cfgs/", I6.6, ".cfg")'),N
           OPEN(UNIT=LP,FILE=FILENAME,STATUS='UNKNOWN')
           write(LP,10)'Number of particles = ',NP
           write(LP,11)'A = 1.0'
           write(LP,12)'H0(1,1) =', cbox(1)
           write(LP,12)'H0(1,2) =', 0.
           write(LP,12)'H0(1,3) =', 0.
           write(LP,12)'H0(2,1) =', 0.
           write(LP,12)'H0(2,2) =', cbox(2)
           write(LP,12)'H0(2,3) =', 0.
           write(LP,12)'H0(3,1) =', 0.
           write(LP,12)'H0(3,2) =', 0.
           write(LP,12)'H0(3,3) =', cbox(3)
           write(LP,11)'.NO_VELOCITY.'
           write(LP,11)'entry_count = 13'
           write(LP,11)'auxiliary[0]= Atom #'
           write(LP,11)'auxiliary[1]= Thermostat'
           write(LP,11)'auxiliary[2]= Temperature [K]'
           write(LP,11)'auxiliary[3]= Molecule number [#]'
           write(LP,11)'auxiliary[4]= Reaction number [#]'
           write(LP,11)'auxiliary[5]= Force [A/fs]'
           write(LP,11)'auxiliary[6]= Ave_velocity [A/fs]'
           write(LP,11)'auxiliary[7]= Ave_acceleration [A/fs]'
           write(LP,11)'auxiliary[8]= DEP_velocity [A/fs]'
           write(LP,11)'auxiliary[9]= DEP_acceleration [A]'
!          Loop over all atom types
           DO K = 1,NTYPES
              IF( NOA(K).GT.0 ) THEN
                 WRITE(LP,*) XMASS(K)
                 WRITE(LP,*) ATYPE(K)
                 DO I=1,NP
                    IF ( KTYPEi(I).EQ.K) THEN
!                      Calculate fractional coordinates  
                       VSQ = 0.0d0
                       ASQ = 0.0d0         
                       DO D =1,3
                          R0F(D) = R0i(D,I)/CBOX(D) + .50d0
                          R1V(D) = R1i(D,I)/DELTA
                          VSQ = VSQ + R1V(D)**2
                          R2A(D) = R2i(D,I)/DELTSQ
                          ASQ = ASQ + R2A(D)**2
                       ENDDO 
!                      Adjust beam molecules
                       IF(R0F(DEPD).GT.1.d0.OR.R0F(DEPD).LT.0.d0) THEN
                         R0F(DEPD) = 0.0d0
                         R1V(1:3) = (/0.0d0,0.0d0,0.d0/)
                         R2A(1:3) = (/0.0d0,0.0d0,0.d0/)
                         ASQ = 0.0d0
                       ENDIF
!                      Caclulate force 
                       RNPSQ = 0.0d0
                       DO  J=1,3
                          RNPSQ = RNPSQ + RNP(J,I)**2
                       ENDDO
                       PRNRNP = SQRT(RNPSQ)
                       XX = VSQ*XMASS(KTYPEi(I))
                       TEMPKi = ENPR*XX/6.0d0/DELTSQ*ECONV 
!              WRITE(LP,13)R0F(1:3),I,ITRi(I),TEMPKi,MIMOLi(I),RCLIST(I)
!     &          ,PRNRNP,SQRT(VSQ),SQRT(ASQ),R1V(DEPD),R2A(DEPD)
              WRITE(LP,13)R0F(1:3),I,ITRi(I),TEMPKi,IMOLi(I),RCLIST(I)
     &          ,PRNRNP,SQRT(VSQ),SQRT(ASQ),R1V(DEPD),R2A(DEPD)
                    ENDIF
                 ENDDO
              ENDIF
           ENDDO
           CLOSE(LP)
        ENDIF
      ENDIF
!
  100 FORMAT(I7)
  300 FORMAT(3E20.11)
  350 FORMAT(I8,I5,3E20.11,I3,2I8)
  360 FORMAT(I8,3E20.11)
  110 FORMAT(I7,3E20.11,I8)
 1800 FORMAT(20A2)
10    format(a,1x,i6)
11    format(a)
12    format(a,1x,f26.6)
13    format(3(E12.4,1x),I8,I5,E12.4,1x,2I8,5(E12.4,1x))
14    format(3E20.11)

      return
      end

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Write out coord.d file for substrate after gas has been cleared
!
      SUBROUTINE WRITE_SUB
!
      USE MPIvars
      USE ANALYSIS
      USE STRUCTR
      USE POTS
      USE SPECIF
!
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
      INTEGER :: N,I
!
      IF ( mynode.EQ.0 ) THEN
        OPEN(UNIT=22,FILE='sub.d',STATUS='unknown')
        WRITE(22,1800) CTITLE
        WRITE(22,100) NP
        WRITE(22,300) TTIME,DELTA
        WRITE(22,300) (CUBE_r(N),N=1,3)
        WRITE(22,300) (CUBE(N),N=1,3)
        DO I=1,NP
           WRITE(22,350) I,KT2(KTYPEi(I)),
     &     (R0i(N,I),N=1,3),ITRi(i)
           WRITE(22,360) I,((R1i(N,I)/DELTA),N=1,3)
           WRITE(22,360) I,((R2i(N,I)/DELTSQ),N=1,3)
           WRITE(22,360) I,(R3i(N,I),N=1,3)
        ENDDO
        CLOSE(22)
      ENDIF
!
  100 FORMAT(I7)
  300 FORMAT(3E20.11)
  350 FORMAT(2I7,3E20.11,I3)
  360 FORMAT(I7,3E20.11)
 1800 FORMAT(20A2)
      RETURN
!
      END SUBROUTINE WRITE_SUB
!
!
      SUBROUTINE write_error

      USE dyn_array
      USE MPIvars
      USE STRUCTR
      USE MDSTP
      USE POTS
      USE SPECIF

      IMPLICIT none

      INCLUDE 'mpif.h'
     
C Local variables
      CHARACTER*20 filename
      INTEGER  :: nfile,NSHIFT,NDUNIT
     &,I,N,ierr

C Write current possitions 
      NSHIFT = 30


      if (nprocs+1.lt.10) then
          nfile=12  
          filename='error_node_0'
          write (filename(nfile+1:),'(i1)') mynode
          nfile = nfile+1
C      elseif (nprocs+1.lt.100) then
C          filename(nfile+1:) ='0'
C          nfile = nfile+1
C          write (filename(nfile+1:),'(i2)') in_beam
C          nfile = nfile+2
      else
          write (filename(nfile+1:),'(i2)') mynode
          nfile = nfile+2
      endif
      filename(nfile+1:) = '.d'
      NDUNIT = mynode+NSHIFT
      OPEN(mynode+NSHIFT,file=filename,STATUS='UNKNOWN')

      IF (mynode.EQ.0) THEN
        WRITE(NDUNIT,1800) CTITLE
        WRITE(NDUNIT,100) NP
        WRITE(NDUNIT,300) TTIME,DELTA
        WRITE(NDUNIT,300) (CUBE_r(N),N=1,3)
        WRITE(NDUNIT,300) (CUBE(N),N=1,3)
      ENDIF
        DO 10 I=1,NP_node
           WRITE(NDUNIT,350) NA(i),KT2(KTYPE_node(I)),
     &          (R0_node(N,I),N=1,3),itr_node(i)
           WRITE(NDUNIT,360) NA(i),((R1_node(N,I)/DELTA),N=1,3)
           WRITE(NDUNIT,360) NA(i),(R2_node(N,I)/DELTSQ,N=1,3)
           WRITE(NDUNIT,360) NA(i),(R3_node(N,I),N=1,3)
10      CONTINUE
C      CLOSE(NDUNIT)

c$$$      DO 321 I=1,NP_node !K=1,KEND
c$$$        jmin=nabors(i)
c$$$        jmax=nabors(i+1)-1
c$$$        DO 322 k=jmin,jmaxc$$$          J=list(k)
c$$$          if(lcheck(k).eq.0) go to 322
c$$$          IF (I.GE.J) GOTO 322
c$$$C           I=IVCT2B(K)
c$$$C           J=JVCT2B(K)
c$$$C          IF(NA(I).GE.NA(J)) GO TO 322
c$$$          DO 323 L=1,3
c$$$C            RNP_node(L,I)=RNP_node(L,I)+RPP(L,K)
c$$$C            RNP_node(L,J)=RNP_node(L,J)-RPP(L,K)
c$$$            RNP(L,NA(I))=RNP(L,NA(I)) + RPP(L,K)
c$$$            RNP(L,NA(J))=RNP(L,NA(J)) - RPP(L,K)
c$$$323       CONTINUE
c$$$322     CONTINUE
c$$$321   CONTINUE

       CALL DEALLOCATE_array
       CALL close_file
       CALL MPI_FINALIZE(ierr)
       STOP

  100 FORMAT(I7)
  200 FORMAT(4F12.6)
  300 FORMAT(3E20.11)
  350 FORMAT(2I7,3E20.11,I3)
  360 FORMAT(I7,5x,3E20.11)
 1200 FORMAT('NEIGHBOR LIST UPDATES: ',/)
 1300 FORMAT(10I4,/)
 1400 FORMAT(8X,'ENERGY(eV)',5X,'T',6X,'TSTEP(fs)  TIME(fs)',/)
 1500 FORMAT(I6,F14.8,F10.3,2F9.3)
 1800 FORMAT(20A2)
!
      END SUBROUTINE write_error
