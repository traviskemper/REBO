C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                              C
C  THESE PROGRAMS PERFORM MOLECULAR DYNAMICS SIMULATIONS WITH  C
C  BOND-ORDER POTENTIALS FOR HYDROCARBON AND OXYGEN,           C
C  AND LENNARD-JONES WITH INPUTTED PARAMETERS.                 C
C                                                              C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C   Units: mass = AMU's, length = Angstroms, time = fs         C
C          energy = eV's                                       C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                              C
C   Details about development and properties of potential      C
C   is described at :                                          C
C   1. CH                                                      C
C      D.W. Brenner*, O.A. Shenderova,                         C
C      J.A. Harrison, S.J. Stewart, B. Ni, S.B. Sinnott,       C
C      Journal of Physics: Condensed Matter 14,                C
C      783-802 (2002).                                         C
C   2. CHO                                                     C
C      Boris Ni, Ki-Ho Lee and Susan B. Sinnott                C
C      J. Phys.: Condens. Matter 16 (2004) 7261-7275           C
C                                                              C
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

      INCLUDE 'module_dynarray.f'
      INCLUDE 'module_MPI.f'
      INCLUDE 'module_TM.f'
      INCLUDE 'module_parameters.f'
      INCLUDE 'module_specif.f'
      INCLUDE 'module_STRUCTR.f'
      INCLUDE 'module_pot.f'
      INCLUDE 'module_tables.f'
      INCLUDE 'module_contm.f'
      INCLUDE 'module_MDSTP.f'
      INCLUDE 'module_analysis.f'
      INCLUDE 'module_beam.f'

      USE dyn_array
      USE MPIvars
      USE TM
      USE PARAMS
      USE SPECIF 
      USE STRUCTR
      USE POTS
      USE TABS
      USE CONTM
      USE MDSTP
      USE ANALYSIS
      USE BEAM
!
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
! MPI variables
!
      INTEGER  dest,tag,ierr,count
      INTEGER  source,status(MPI_STATUS_SIZE)
      INTEGER  start,J,I
      DOUBLE PRECISION :: startt,stopt,RNPSQ,PRNRNP
      REAL*8 :: PRDISTtmp,PRENtmp
!
! INITIAL MPI
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,mynode,ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Initialize simulation 
      CALL WRTHEAD    !Write head of output file
      CALL STDFLT
c open input/output files
       CALL open_file
       IF (mynode.EQ.0) WRITE(*,*)'Files open'
c initialize
      call setin
C      IF (mynode.EQ.0) WRITE(*,*)'setin in finished'
C read input data
C      startt=MPI_WTIME(ierr)
!      call read_input
       CALL READIN
!      WRITE(*,*) mynode,'Spacial decomp in main',NDP(:) 

C      IF (mynode.EQ.0) WRITE(*,*)'input read'
C      CALL write_errort=MPI_WTIME(ierr)
C      WRITE(*,*) 'node',mynode,'read_data',CALL write_errort-startt
C setup potential parameters
C      startt=MPI_WTIME(ierr)
      call setpp !param,mtablt,ljparam,ljcset
C      IF (mynode.EQ.0) WRITE(*,*)'Potential set'
C      CALL write_errort=MPI_WTIME(ierr)
C      WRITE(*,*) 'node',mynode,'setpp',CALL write_errort-startt
c setup predictor-corrector coefficients
C      startt=MPI_WTIME(ierr)
      call setpc
C      IF (mynode.EQ.0) WRITE(*,*)'Predictor corrector set'
C      CALL write_errort=MPI_WTIME(ierr)
!      WRITE(*,*) 'node',mynode,'setpc'!,CALL write_errort-startt
C read coord data
C      startt=MPI_WTIME(ierr)
      call read_coord
C      CALL write_errort=MPI_WTIME(ierr)
!      WRITE(*,*) 'node',mynode,'read_data' !,CALL write_errort-start
!     Set variables for post analysis
!      IF(RECAN) CALL checksub
      IF(MOLAN.OR.RECAN.OR.CLRGAS) CALL set_anly
!!      WRITE(*,*) 'node',mynode,'set_anly'
! Initialize cfg
      IF(PCFG) CALL SETCFG
C find neighbor nodes
      CALL nebor_node
C      startt=MPI_WTIME(ierr)
      CALL pass_atom_node
C      CALL write_errort=MPI_WTIME(ierr)
!      WRITE(*,*) 'node',mynode,'pass_atom_node'!,CALL write_errort-startt
C setup Langevin parameters
      call setgle
!      WRITE(*,*) 'node',mynode,'setgle'
C initialize random number generator
      call setran
!      WRITE(*,*) 'node',mynode,'setgle'
C find neighbor neighbor cells
      CALL cagut_nabor_cell
!      WRITE(*,*) 'node',mynode,'cagut_nabor_cell'
C      WRITE(*,*) 'Caguts nabor cell'
      CALL ljgut_nabor_cell 
!      WRITE(*,*) 'node',mynode,'ljgut_nabor_cell'
C      CALL write_errort=MPI_WTIME(ierr)
C      WRITE(*,*) 'node',mynode,'arrange_atom',CALL write_errort-startt
!      Write simulation options
       CALL WRITEOPT
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Begin calculation *
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      if(kflag.eq.6) then
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Energy minimization *
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
           kvc = 1
C           call minimize
           write(6,*) 'minimum energy= ',tote
      else
c**********************
C Dynamic Simulation  *
!**********************
!
            call setmd
!      WRITE(*,*) 'node',mynode,'setmd'
C Open optional output files
            DO I=1,NPtot_node
              IF(PRNPR.and.NA(I).EQ.PR1)THEN
                OPEN(UNIT=23,FILE='pair.out',STATUS='UNKNOWN')
                OPEN(UNIT=24,FILE='force.out',STATUS='UNKNOWN')
                OPEN(UNIT=25,FILE='pair.dat',STATUS='UNKNOWN')
              ENDIF
            ENDDO
!      WRITE(*,*) 'node',mynode,'optional open'
C Start the MD loop
            DO LSTEP=1,KVC
!      WRITE(*,*) 'node',mynode,'main loop'
c predictor
C      startt=MPI_WTIME(ierr)
              call cpred
C      CALL write_errort=MPI_WTIME(ierr)
!      WRITE(*,*) 'node',mynode,'cpred'!,CALL write_errort-startt
c calculate energy and forces
              CALL MODEL
C              WRITE(*,*) mynode,'LCHK=',LCHK
CCC Each node make it's neighbor list
              IF(LCHK.EQ.1) THEN
C      startt=MPI_WTIME(ierr)
                CALL pass_atom_node
C      CALL write_errort=MPI_WTIME(ierr)
!      WRITE(*,*) 'node',mynode,'pass_atom_node'!,CALL write_errort-startt
C      startt=MPI_WTIME(ierr)
                CALL colist
C      CALL write_errort=MPI_WTIME(ierr)
!      WRITE(*,*) 'node',mynode,'colist'!,CALL write_errort-startt
C      startt=MPI_WTIME(ierr)
                CALL LJVLIST
C      CALL write_errort=MPI_WTIME(ierr)
!      WRITE(*,*) 'node',mynode,'ljvlist'!,CALL write_errort-startt
              ENDIF
              DO I=1,NTYPES
                sss = sss + noa(I)
              ENDDO 
!              sss = noa(1)+noa(2)+noa(3)+noa(4)+noa(5)
C      startt=MPI_WTIME(ierr)
              if((ipot.eq.1).and.(sss.ne.0)) CALL CAGUTS
C      CALL write_errort=MPI_WTIME(ierr)
!      WRITE(*,*) 'node',mynode,'caguts'!,CALL write_errort-startt
C      startt=MPI_WTIME(ierr)
C Add two versions of pibond for CHF and CHO 
c  -travisk
            IF ( SYTP.EQ.1) THEN
                CALL PIBONDCHF
C              WRITE(*,*) mynode,'Potential type SYTP=',SYTP
            ELSEIF ( SYTP.EQ.2) THEN
                CALL PIBONDCHO
C              WRITE(*,*) mynode,'Potential type SYTP=',SYTP
            ELSEIF ( SYTP.EQ.3) THEN
                CALL PIBONDCHS
C              WRITE(*,*) mynode,'Potential type SYTP=',SYTP
            ELSE
              WRITE(*,*) mynode,'Potential type SYTP='
     &               ,SYTP,' is undefined'
              WRITE(*,*) noa(:)
              WRITE(*,*) SYTP
            ENDIF  

C      CALL write_errort=MPI_WTIME(ierr)
!      WRITE(*,*) 'node',mynode,'pibond'!,CALL write_errort-startt
C              if(noa(5)+noa(4).ne.0) call sili_germ
C              if((ipot.eq.2).and.((noa(1)+noa(2)+noa(3)).ne.0))
C     &        call tight_bind
C      startt=MPI_WTIME(ierr)
              IF(ILJ) call ljguts
C      CALL write_errort=MPI_WTIME(ierr)
!      WRITE(*,*) 'node',mynode,'ljguts'!,CALL write_errort-startt
C      startt=MPI_WTIME(ierr)
C              CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
C              CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
              IF( ADDFC.AND.mod(LSTEP,FSTEP).EQ.0) CALL ADD_FORCE
              Call pass_force_node
!
              IF( PRNLD) CALL load
!             Heat system
              IF ( HT ) CALL HEAT
!             Apply thermostats
              call thermos
!             Corrector
              call ccorr
!             move_tip
              IF (MVTIP.AND.mod(LSTEP,isteps).eq.0) call move_tip

C              WRITE(*,*) 'NMA=',NMA,'NP=',NP
C              IF (NMA.NE.NP) CALL MOVCTR
              CALL CHKNAB
!      WRITE(*,*) 'node',mynode,'CHKNAB'
c write out position file to be post converted to xmol format
C      startt=MPI_WTIME(ierr)
              IF (mod(LSTEP,NSTR).eq.0) THEN
                call reduce_coord
                IF(MOLAN.OR.RECAN) THEN
                  CALL molcheck
                  IF(MOLAN) CALL MASSPEC
                  IF(RECAN) CALL react
                  IF(CRSAN) CALL CRSLINK
                  IF(SMOLAN) CALL smolcheck
                  IF(SMOLAN) CALL surfmspec
                  IF(DEPPROF ) CALL DEPTHPROF
                ENDIF
                CALL write_coord
              ENDIF
!      WRITE(*,*) 'node',mynode,'write coord/xmol'
C      CALL write_errort=MPI_WTIME(ierr)
C      WRITE(*,*) 'node',mynode,'xmol',CALL write_errort-startt
C
C              WRITE(*,*) NMA,NP
              IF (KFLAG.EQ.5) CALL BERE
c
C       startt=MPI_WTIME(ierr)
              IF(mod(LSTEP,maxkb).eq.0) then
c generate and write data
                call write_data
              ENDIF
!      WRITE(*,*) 'node',mynode,'write 2'
C      CALL write_errort=MPI_WTIME(ierr)
C      WRITE(*,*) 'node',mynode,'write_data2',CALL write_errort-startt 
c volume scaling
c               call vscale
C             WRITE(*,*) 'Loop=',LSTEP,'finished'
C WRITE OUT POSITIONS FOR RESTART
C             call write_data3
           ENDDO
      ENDIF

      IF(mynode.EQ.0)WRITE(6,*)'Molecular dynamics simulation complete'
C       startt=MPI_WTIME(ierr)
C      CALL write_errort=MPI_WTIME(ierr)
C      WRITE(*,*) 'node',mynode,'write_data3',CALL write_errort-startt
C       CALL write_error
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Start post analysis
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(mynode.EQ.0)WRITE(6,*) 'Starting post analysis'
      IF( RECSUB ) CALL PRNSUBE !record substrate information for deposition
!
      IF ( CLRGAS )  call reduce_coord
      IF(mynode.EQ.0)  THEN
        IF ( CLRGAS.OR.MOLAN ) THEN
          WRITE(6,*) 'Starting post analysis'
          CALL molcheck
          IF( MOLAN ) THEN
            CALL MASSPEC
          ENDIF
          IF ( CLRGAS ) THEN
!         Move all gas phase molecules under substrate with v=0
            WRITE(6,*) ' Clearing gas molecules'
            call cleargas
            call write_coord
          ENDIF
        ENDIF
      ENDIF
      IF(mynode.EQ.0)WRITE(6,*) 'Finished post analysis'
!
      CALL DEALLOCATE_array
      CALL close_file
      CALL MPI_FINALIZE(ierr)
!
      STOP 
      END
C
C add included subroutines
C
      include 'subroutines.inc' 
