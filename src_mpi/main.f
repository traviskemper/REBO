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
      INCLUDE 'module_mpi.f'
      INCLUDE 'module_time.f'
      INCLUDE 'module_parameters.f'
      INCLUDE 'module_specif.f'
      INCLUDE 'module_structure.f'
      INCLUDE 'module_pot.f'
      INCLUDE 'module_tables.f'
      INCLUDE 'module_contm.f'
      INCLUDE 'module_md.f'

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
!     INITIAL MPI
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,mynode,ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Initialize simulation 
      CALL WRTHEAD    !Write head of output file
      CALL STDFLT
!     Open input/output files
      CALL open_file
      IF (mynode.EQ.0) WRITE(*,*)'Files open'
!     Initialize
      call setin
!
!     read input data
      CALL READIN
!
!     Setup potential parameters
      call setpp !param,mtablt,ljparam,ljcse
!
!     setup predictor-corrector coefficients
      call setpc
!
!     read coord data
      call read_coord
!
!     Initialize cfg
      IF(PCFG) CALL SETCFG
!
!     find neighbor nodes
      CALL nebor_node
      CALL pass_atom_node
!
!     setup Langevin parameters
      call setgle
!
!     initialize random number generator
      call setran
!
!     find neighbor neighbor cells
      CALL cagut_nabor_cell
      CALL ljgut_nabor_cell 
!
!     Write simulation options
      CALL WRITEOPT
!
!     Dynamic Simulation  
!
      CALL setmd
!     Start the MD loop
      DO LSTEP=1,KVC
!
!        Predictor
         call cpred
!
!        Calculate energy and forces
         CALL MODEL
!        Each node make it's neighbor list
         IF(LCHK.EQ.1) THEN
            CALL pass_atom_node
            CALL colist
            CALL LJVLIST
         ENDIF
!
!        Calculate pair forces
         CALL CAGUTS
!
!        Use proper pibond CHF/ CHO /CHS for bij term 
         IF ( SYTP.EQ.1) THEN
            CALL PIBONDCHF
         ELSEIF ( SYTP.EQ.2) THEN
            CALL PIBONDCHO
         ELSEIF ( SYTP.EQ.3) THEN
            CALL PIBONDCHS
         ELSE
            WRITE(*,*) mynode,'Potential type SYTP='
     &               ,SYTP,' is undefined'
            WRITE(*,*) noa(:)
            WRITE(*,*) SYTP
         ENDIF  
!         
!        Calculate lennard Jones forces 
         IF(ILJ) call ljguts
!
         Call pass_force_node
!
         IF( PRNLD) CALL load
!        Heat system
         IF ( HT ) CALL HEAT
!        Apply thermostats
         call thermos
!        Corrector
         call ccorr
!        move_tip
         IF (MVTIP.AND.mod(LSTEP,isteps).eq.0) call move_tip
!
         CALL CHKNAB
         IF (mod(LSTEP,NSTR).eq.0) THEN
            call reduce_coord
            CALL write_coord
         ENDIF
!
         IF (KFLAG.EQ.5) CALL BERE
!
         IF(mod(LSTEP,maxkb).eq.0) then
!        Generate and write data
            call write_data
         ENDIF
      ENDDO
!
      IF(mynode.EQ.0)WRITE(6,*)'Molecular dynamics simulation complete'
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      CALL DEALLOCATE_array
      CALL close_file
      CALL MPI_FINALIZE(ierr)
!
      STOP 
      END
!
!     add included subroutines
!
      include 'subroutines.inc' 
