C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Programer:
C     Travis Kemper
C     Department of Materials Science and Engineering
C     University of FLorida
C     traviskemper@ufl.edu
C
C     Version 1.0 9/25/09 T. W. Kemper
C     Version 5.2 10/06/09 T. W. Kemper
C          - Modify to track multible clusters
C     Version 9.7 11/09/09 T. W. Kemper
C          - Combine NN list and molecule list into analysis
C          - Print out molecular formula 
C            - subroutine: PRNEL
C          - Print out single reaction xmol.d file x.rec.xyz
C            - print out reaction with atom # of x.rec.xyz  
C     Version 9.8 11/17/09 T. W. Kemper
C          - Get rid of unused variables  
C          - Use dynamic array allocation
C     Version 11.4 12/09/09 T. W. Kemper
C          - add mass spectra and molecular depth profile
C          - prnel and comas made seperate subroutines
C     Version 11.6 12/29/09 T. W. Kemper
C          - add depth profile of atoms and molecules
C     Version 12 12/29/09 T. W. Kemper
C          - Update analysis to find backbone atoms
C     Version 13 01/18/10 T. W. Kemper
C          - Add depth profile of broken C-C 
C             and cross links
C     Version 14 01/18/10 T. W. Kemper
C          - Add max depth 
C     Version 15 01/18/10 T. W. Kemper
C          - copy of 14
C          - something  went wrong with v15 
C          - Add PS atom identifier
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Input files:
C       mol.xyz    !Molecule 
C       coord.d    !Finial structure file from simulation
C       beam.d     !Initial structure file
C     Output files:
C       rect.txt   !Reaction information  
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Read in initial and finial structures and molecular file
C Find beam atoms which have bonded to substrate
C Find nieghbors of bonded beam atoms 
C Compare to their orginal bonded states
C Compline list of atoms involved in reaction
C Find reaction in xmol.d file and check for transitions 

      INCLUDE '../module_rebo.f'
      INCLUDE 'module_reactv2.f'

      PROGRAM depanl

      USE rebo
      USE react

      IMPLICIT none

      INTEGER :: LMX
      PARAMETER ( LMX = 100 )
      INTEGER, DIMENSION(RMX) :: MRICNT,BMOL,GNTEMP,BCKBN
      INTEGER, DIMENSION(RMX,RMX) :: MRILST
      INTEGER :: A,E,K,O,M,B,C,P,N,I,Q,D,J,L,Qf,Di,Ef,Ei,Qi,T,IOS,
     & Io,Jo,Oo,IM,Nf,Bf,No,Bo,Of,Cf,Oi,Pf,Ii,Ci,Df,Ji,Po,Qo,Mi,H,
     & GN,BC,BK,Ff,Gf,Hf,Jf,Kf,Lf,Rf,Mf,Sf,Tf,Uf,Vf,Xf,
     & BBCNT,NBCNT,BCCNT,
     & CNT,
     & MOLCNTI,MOLCNTF,
     & MCNT,SCNT,CRSCNT,
     & RCNTS(AMX,2),PCNTS(AMX,2),
     & ELIST(LMX),ECNT,ELCNT(ELMX),
     & ATPS,Ityp
     & ,BMCNT,BKCNT 
     & ,ACNK,CNK(AMX)
     & ,RXC,RYC,RZC,RX,RY,RZ,NACS,IC,BKBND(12),BKID
      INTEGER, DIMENSION(AMX,NMX) :: NLISTI,NLISTF
      DOUBLE PRECISION ::  RCUTSQ,RHCUTSQ,HHCUTSQ,
     &  HBUF,DBUF,SMN(2,3),SMIN,SMAX,NBEAM,DEN,CNKMN,YCNK
      INTEGER , DIMENSION(AMX) :: BBLST,BCLST,NBLST,
     & MOLSTI,MPNTI,MOLSTF,MPNTF,
     & NCTI,NCTF,ITYPEM
      DOUBLE PRECISION, DIMENSION(5,AMX,3) :: RINT,RFIN
      DOUBLE PRECISION, DIMENSION(LMX,3) :: RM
      INTEGER :: NAM,NC,NAS,NPROD,NRECT,PNAM,PNC,NADEP
      CHARACTER(7) :: FIN,FOUT,MTITLE
      CHARACTER(2) :: ID
      CHARACTER(30) :: PREL,CDUM
      CHARACTER(16), DIMENSION(:), ALLOCATABLE :: PRAN
      CHARACTER(30), DIMENSION(:), ALLOCATABLE :: RECTANTS,PRODUCTS
      LOGICAL :: RECED(AMX),ASUBI(AMX),ASUBF(AMX),VERBOSE,XML
     &  ,BFBCK,RPREV,ERR1
      INTEGER :: ITYPEI(AMX),ALABI(AMX),ALABF(AMX),ADUM,BDUM
      DOUBLE PRECISION :: PETi,PETf,DPET,MAXDP(2),DPETOT
C Molecular analysis
      INTEGER :: MT,MTCNT,MNCNT(2,MMX),NM,NMB,NMU,DEP,MTPS
     & ,NABA,NAA,NABB,NAB,NT,CNCNT,HNCNT,BTPS,
     & DMN,DMX,DP,GPCNT(MMX),NN,NNTEMP(NMX),BANL(AMX,4)
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: MDCNT
      INTEGER, DIMENSION(:), ALLOCATABLE :: NMG,NAG,CRSD,NCNT
      DOUBLE PRECISION :: MMAS(MMX),DEV
     &  ,COM(3),MOLM,NAAR,NABR,NTR
     &  ,Hrl,Drl,DLST(4,MMX)
      CHARACTER(30), DIMENSION(MMX)::MMFRM
      CHARACTER(4), DIMENSION(2) :: PREF
      LOGICAL :: PBX

C Clean substrate
      INTEGER BCNT,CCNT
      INTEGER, DIMENSION(AMX) ::TYPSUB,THSUB,TYPC
      INTEGER, DIMENSION(5,AMX,3) ::RB,RC


C Read in booleans 
C          [true if parellel] [Read previous ] 
      READ(*,*)   PARL,            RPREV

C Set output
      VERBOSE=.FALSE.
      IF ( PARL ) THEN
       WRITE(*,*)'Analysizing parallel output beam.d out.d and xmol.d'
      ELSE
       WRITE(*,*)'Analysizing serial output beam.d coor.d and xmol.d'
      ENDIF 

C Should be set as parameters or data but not sure how to neatly do that
       ELTRN(6)  = 1
       ELTRN(1)  = 2
       ELTRN(8)  = 3
       ELTRN(9)  = 4
       ELTRN(16) = 5

C Set squares of cut offs
      RCUTSQ = RCUT*RCUT
      RHCUTSQ = RHCUT*RHCUT
      HHCUTSQ = HHCUT*HHCUT

C   Set buffer region
      DBUF = 7.0  !PMMA about 1/2 a coil
      HBUF = 20.  !Max hight to be counted as substrate

C Open needed files
      OPEN(UNIT=20,file='mol.xyz',status='old')
      OPEN(UNIT=22,FILE='output.d',STATUS='unknown')
      OPEN(UNIT=21,FILE='rect2.txt',STATUS='unknown')

C  Read in molecule structure
      WRITE(*,*) 'Reading in structures'
      READ(20,*) NAM,NC,DEN  !# atoms in mol.,# of cluster/beam
      READ(20,*) MTITLE
      DO I=1,NAM
         READ(20,*) ITYPEM(i),RM(i,1),RM(i,3),RM(i,3)
      ENDDO
      CLOSE(20)
      WRITE(*,*) 'Deposition of ',NC,' ',MTITLE,'at',DEN,'eV'
 
C Read in intial structure
      fin='beam.d'
      IF ( PARL ) THEN
        CALL readprlf(fin)
        BOX=BOX2   !Reset box for nborlst and molcheck
      ELSE
        CALL readcrdf(fin)
      ENDIF
      WRITE(*,*) 'Initial structure has',NA,' atoms'
      CALL nborlst      !Create intial neighbor list
      CALL anbond(ALABI)
      CALL molcheck(MOLCNTI,MPNTI,MOLSTI,ASUBI)
      NLISTI = NLIST
      NCTI   = NCT
      RINT   = R
      ITYPEI=ITYPE

C Read in final structure from deposition
      IF ( PARL ) THEN
        fin='out.d'
        CALL readprlf(fin)
        BOX=BOX2   !Reset box for nborlst and molcheck
      ELSE
        fin='coord.d'
        CALL readcrdf(fin)
      ENDIF
      WRITE(*,*) 'Finial structure has',NA,' atoms'
      CALL nborlst      !Create finial neighbor list
      CALL anbond(ALABF)
      CALL molcheck(MOLCNTF,MPNTF,MOLSTF,ASUBF)
      NLISTF = NLIST
      NCTF   = NCT
      RFIN   = R
      WRITE(*,*) '  Structures read in'

C Read outputfile
      DO L=1,10
      READ(22,*)
      ENDDO
      IF ( PARL ) THEN
      READ(22,*) A,PETi
      ELSE
      READ(22,*) A,B,T,PETi
      ENDIF 
      DO
        IF ( PARL ) THEN
          READ(22,*,IOSTAT=ios)  A,PETf
        ELSE
         READ(22,*,IOSTAT=ios) A,B,T,PETf
        ENDIF
        IF ( ios .gt. 0 ) THEN
          WRITE(*,*) 'Something is wrong with output.d'
          WRITE(*,*) 'at step ',A
          EXIT
        ELSE IF (ios .lt. 0 ) THEN
          EXIT
        ENDIF 
      ENDDO
      CLOSE(22) 
      DPET= (PETf-PETi)*100.0D0 
      WRITE(*,*) 'Change in Pot En:',DPET,' (meV)'
 
C Check to see if any atoms changed type
C   make sure parellel code is read incorrectly since the order 
C   gets mixed up
      DO I=1,NA
        IF ( ITYPEI(I).NE.ITYPE(I) ) THEN
          WRITE(*,*) 'Atom ',I,'changed type'
          WRITE(*,*) 'from a',ITYPEI(I), 'to a',ITYPE(I)
          STOP
        ELSE
        ENDIF      
      ENDDO

C Calculate numbers of atoms in substrate
      NAS = NA - NAM*NC

C Get max and min position of subtrate atoms and
      SMN(1,:) = (/1000.,1000.,1000./)
      SMN(2,:) = (/-1000.,-1000.,-1000./)
      ACNK  = 0
      CNKMN = 1000.0
      DO I=1,NA
        RECTRN(I) = I -1 !Initialize atom list to originial VMD values
        IF ( ASUBF(I) ) THEN
          DO N=1,3
            IF ( RFIN(1,I,N) .GT. SMN(2,N) ) SMN(2,N)=RFIN(1,I,N)
            IF ( RFIN(1,I,N) .LT. SMN(1,N) ) SMN(1,N)=RFIN(1,I,N)
          ENDDO
          IF ( RFIN(1,I,2).GT.HBUF ) THEN !record chunk
            ACNK  = ACNK + 1
	    CNK(ACNK) = I 
            IF ( RFIN(1,I,2).LE. CNKMN ) CNKMN = RFIN(1,I,2)
          ENDIF
        ENDIF
      ENDDO
C
C Warn if very large piece of substrate has broken off and print 
C  coordinates
      IF (SMN(2,2).GT.HBUF) THEN
        WRITE(*,*) 'Max of substrate is greater than ',HBUF
        WRITE(*,*) 'Likely a  large piece that has broken off'
        OPEN(UNIT=25,FILE='chunk.xyz',STATUS='unknown')
        WRITE(25,*) ACNK,' With min y:',CNKMN
        WRITE(25,*)'coordinates of large chunk which has broken off'
        DO C=1,ACNK
          I =  CNK(C) 
          YCNK  = RFIN(1,I,2)  - CNKMN
          WRITE(25,'(I7,3e16.8)') ITYPE(I),RFIN(1,I,1),YCNK,RFIN(1,I,3)
          ASUBF(I) = .FALSE.
        ENDDO 
        CLOSE(25)
      ENDIF

C Add buffer
      SMIN = SMN(1,2) - DBUF
      SMAX = -DBUF
C Initialize various counts
      BCNT  = 0   !Total number of atoms in the bulk
      CCNT  = 0   !Total number of atoms to be deleted 
      BBCNT = 0   !Number of beam in substrate atoms bonded
      NBCNT = 0   !Number of beam in substrate atoms nonbonded
      BCCNT = 0   !Number of beam atoms not in substrate

C Intilized atomic depth profile
      DEV = -5.0
      DMN = -4 !Hard set for PS with DEV=5.
      DMX = 20 !Hard set for PS
      MTPS = 500 !Max number of molecule types
      ATPS= 2    !number of types of atoms
      ALLOCATE(MDCNT(2,ATPS,DMN:DMX))
      ALLOCATE(CRSD(DMN:DMX))        !Crosslink depth profile
      ALLOCATE(NMG(MTPS))
      ALLOCATE(NAG(ATPS))

C Zero depth count
      DO A=1,ATPS
        DO T=1,2           !Bonded =1 Unbonded=2
          DO D=DMN,DMX
            MDCNT(T,A,D) = 0
            CRSD(D) = 0
          ENDDO
        ENDDO
        NAG(A) = 0
      ENDDO
      DO A=1,MTPS
        DO T=1,2        !Bonded =1 Unbonded=2 
          MNCNT(T,A) = 0 
        ENDDO
        NMG(A) = 0
      ENDDO
      MCNT = 0
      CRSCNT = 0        !Number of crosslinks

C Read in information from previous run
      IF ( RPREV ) THEN
       WRITE(*,*) 'Read in ../masspec.txt'
       OPEN(UNIT=27,FILE='../masspec.txt',STATUS='unknown')
       READ(27,*)
       READ(27,*) PNAM,PNC,NADEP  
       READ(27,*)
       READ(27,*)
       READ(27,*) CDUM,NAG(1)
       READ(27,*) CDUM,NAG(2)
       READ(27,*) CRSCNT
       READ(27,*) MTCNT 
       READ(27,*) CDUM,MAXDP(1),CDUM,MAXDP(2)
       READ(27,*) DPETOT
       WRITE(*,*) MTCNT,' molc types previously found'
       READ(27,*) 
       DO M=1,MTCNT
        READ(27,*) MT,MMFRM(MT),MMAS(MT),CDUM,NMB,NMU,NMG(MT),NM 
        MNCNT(1,MT) = NMB
        MNCNT(2,MT) = NMU 
        WRITE(*,*) M,' ',MMFRM(MT),NM
       ENDDO
       CLOSE(27)
       WRITE(*,*) 'Read in ../ST1.dat, ../ST2.dat and ../crsdep.dat'
       OPEN(UNIT=31,FILE='../ST1.dat',STATUS='unknown')
       OPEN(UNIT=32,FILE='../ST2.dat',STATUS='unknown')
       OPEN(UNIT=33,FILE='../crsdep.dat',STATUS='unknown')
       DO D=DMN,DMX
        READ(31,*) ADUM,BDUM,NABA,NAA
        READ(32,*) ADUM,BDUM,NABB,NAB
        READ(33,*) ADUM,BDUM,CRSD(D)
        MDCNT(1,1,D) = NABA
        MDCNT(1,2,D) = NABB
        MDCNT(2,1,D) = NAA - NABA
        MDCNT(2,2,D) = NAB - NABB
       ENDDO
       CLOSE(31)
       CLOSE(32)
       CLOSE(33)
      ELSE
       MTCNT = 0 
       PNAM  = 0
       PNC   = 0
       NADEP = 0 
       DPETOT = 0.0D0    !Increase in potential energy
       MAXDP(:) = 0.0D0  !Max depth of beam atoms
      ENDIF

      WRITE(*,*) 'Starting analysis'
C Find beam atoms bonded to substrate and atomic depth profile
      DO I=NAS+1,NA
        E = ITYPE(I)
        Ityp = ELTRN(E)
        Hrl = RFIN(1,I,2)
        IF ( E .EQ. 1 ) THEN
          IF ( Hrl.LT.MAXDP(2) ) MAXDP(2) = Hrl
        ELSEIF ( E.EQ. 6 ) THEN
          IF ( Hrl.LT.MAXDP(1) ) MAXDP(1) = Hrl
        ENDIF
        Drl = Hrl/DEV+0.5d0   !1/2
        DEP = NINT(Drl)
        IF  (ASUBF(I)) THEN
          IF ( DEP.LT.DMN .or. DEP.GT.DMX) THEN
            WRITE(*,*) 'Depth of beam atom',I,' at', Hrl,'
     &       bond to substrate is out of
     &      range of depth array. Change array size with DMN,DMX'
            WRITE(*,*)'Which is currently',DMN,DMX
            STOP
          ENDIF
          BBCNT = BBCNT + 1
          BBLST(BBCNT) = I 
          RECED(I) = .FALSE.  !Initialize reaction boolean
          MDCNT(1,Ityp,DEP) = MDCNT(1,Ityp,DEP) + 1
C        ELSEIF ( DEP.LE.DMX .and. Hrl.LE.SMN(2,2) ) THEN
        ELSEIF ( DEP.GE.DMN.and.DEP.LE.DMX) THEN
          IF ( DEP.LT.DMN.or.DEP.GT.DMX) THEN
            WRITE(*,*) 'Depth of beam atom',I 
            WRITE(*,*) 'At hight ',RFIN(1,I,2),Hrl
            WRITE(*,*) 'Which is less than',SMN(2,2)
            WRITE(*,*) 'With index',DEP,'Is out of range',DMN,DMX
            STOP
          ENDIF
          NBCNT = NBCNT + 1
          NBLST(NBCNT) = I  
          MDCNT(2,Ityp,DEP) = MDCNT(2,Ityp,DEP) + 1
        ELSE  !should be in gas
          BCCNT = BCCNT + 1
          BCLST(BCCNT) = I  
          NAG(Ityp) = NAG(Ityp) + 1
        ENDIF
      ENDDO
 
      IF ( VERBOSE ) THEN
        OPEN(unit=29,FILE='alst.txt',STATUS='unknown')
        WRITE(29,*) 'Bonded'
        WRITE(29,*) BBLST(:)
        WRITE(29,*) 'UnBonded'
        WRITE(29,*) NBLST(:)
        WRITE(29,*) 'Gas'
        WRITE(29,*) BCLST(:)
        CLOSE(29)
      ENDIF


C Print gas phase info
      OPEN(UNIT=27,FILE='../masspec.txt',STATUS='unknown')
      OPEN(UNIT=28,FILE='../masspec.dat',STATUS='unknown')
      OPEN(UNIT=31,FILE='../ST1.dat',STATUS='unknown')
      OPEN(UNIT=32,FILE='../ST2.dat',STATUS='unknown')
      NAA = NAG(1) 
      NAB = NAG(2) 
      N  = NAA  +NAB
      NT = NADEP+NAM*NC  !Total # of beam atoms
      WRITE(27,*)  '# of atoms/molc,# of molcs,# of atoms deposited'
      WRITE(27,*)  NAM,PNC+NC,NT
      WRITE(27,*) N,'atoms beam atoms found in gas' 
      WRITE(27,*) NT-N,'atoms in substrate'
      NTR = REAL(NT)
      NAAR = REAL(NAA)
      NABR = REAL(NAB)
      WRITE(27,*) 'C=',NAA,100*NAAR/NTR
      WRITE(27,*) 'H=',NAB,100*NABR/NTR

C  Print out atomic depth profile
        DO DP = DMN,DMX
          Drl=REAL(DP)- 0.5d0
          NABA = MDCNT(1,1,DP)
          NAA  = NABA  + MDCNT(2,1,DP) 
          NABB = MDCNT(1,2,DP)
          NAB  = NABB  + MDCNT(2,2,DP) 
          WRITE(31,*) DEN,Drl*DEV,NABA,NAA
          WRITE(32,*) DEN,Drl*DEV,NABB,NAB
        ENDDO
        WRITE(31,*) ''
        WRITE(32,*) ''
      DEALLOCATE (MDCNT)
      CLOSE(31)
      CLOSE(32)

C Convert the number of atoms in a beam to a float
      NBEAM =  FLOAT(NAM*NC)
C Intilized atomic depth profile
      BTPS= 3    !number of bond types
      ALLOCATE(MDCNT(4,BTPS,DMN:DMX))
      ALLOCATE(NCNT(4))

C Read in values from previous run
      IF (RPREV) THEN
        OPEN(UNIT=34,FILE='../bonddep.dat',STATUS='unknown')
        DO DP= DMN,DMX
          READ(34,*) ADUM,BDUM,TYPC(1:12)
          MDCNT(1,1,DP) = TYPC(1)
          MDCNT(1,2,DP) = TYPC(2)
          MDCNT(1,3,DP) = TYPC(3)
          MDCNT(2,1,DP) = TYPC(4)
          MDCNT(2,2,DP) = TYPC(5)
          MDCNT(2,3,DP) = TYPC(6)
          MDCNT(3,1,DP) = TYPC(7)
          MDCNT(3,2,DP) = TYPC(8)
          MDCNT(3,3,DP) = TYPC(9)
          MDCNT(4,1,DP) = TYPC(10)
          MDCNT(4,2,DP) = TYPC(11)
          MDCNT(4,3,DP) = TYPC(12)
        ENDDO  
        CLOSE(34)
      ELSE
        DO I=1,4
          DO J=1,3
            DO DP= DMN,DMX
              MDCNT(I,J,DP) = 0
            ENDDO
          ENDDO
       ENDDO
      ENDIF

C Intilized atomic depth profile
      BANL(:,:) = 0
      DO I =1,NAS   !Loop over all substrate atoms
C Check for broken bonds
        IF ( ITYPE(I).NE.1 ) THEN
          NCNT(:) = 0
          NCNT(1) = BANL(I,1)
          NCNT(3) = BANL(I,3)
          CNCNT = 0
          HNCNT = 0
          DO 2001 No=1,NCTI(I)
            Bo=NLISTI(I,No)  
            IF ( ITYPE(Bo).EQ.6) CNCNT = CNCNT + 1 
            IF ( ITYPE(Bo).EQ.1) HNCNT = HNCNT + 1 
            DO Nf=1,NCTF(I)
              Bf=NLISTF(I,Nf)         
              IF ( Bf.EQ.Bo ) GOTO 2001           
            ENDDO 
C  This should go through the xmol.d file and find when the bond broke
C    and use the possition at that time as the depth
C    but for intial testing
C   Set bond depth to be that of atom I 
            Hrl = RINT(1,I,2)  
            Drl = Hrl/DEV/+0.5d0
            DEP = NINT(Drl)
            IF ( ITYPE(Bo).EQ.6) THEN
              BANL(Bo,3) = BANL(Bo,3) - 1
              NCNT(3) = NCNT(3) + 1 
              DLST(3,NCNT(3)) = DEP
            ELSE
              NCNT(4) = NCNT(4) + 1 
              DLST(4,NCNT(4)) = DEP
            ENDIF
 2001     CONTINUE 
C Check for new neighbors [formed bonds]
          DO 2002 Nf=1,NCTF(I)
            Bf=NLISTF(I,Nf)
            DO No=1,NCTI(I)
              Bo=NLISTI(i,No)
              IF ( Bf.EQ.Bo ) GOTO 2002
            ENDDO
            Hrl = RINT(1,I,2)  
            Drl = Hrl/DEV/+0.5d0
            DEP = NINT(Drl)
            IF ( ITYPE(Bf).EQ.6) THEN
              BANL(Bf,1) = BANL(Bf,1) - 1
              NCNT(1) = NCNT(1) + 1 
              DLST(1,NCNT(1)) = DEP
            ELSE
              NCNT(2) = NCNT(2) + 1 
              DLST(2,NCNT(2)) = DEP
            ENDIF
 2002     CONTINUE 
C  B
C 1 - C-C formed, 2 - C-H formed, 3 - C-C broken, 4 - C-H broken
C Ityp
C  1 - aromatic C-C
C  2 - on back bone 
C  3 - Unkown 
            IF ( CNCNT.EQ.2.AND.HNCNT.EQ.1 ) THEN  !If arromatic bond
              Ityp = 1
            ELSEIF ( CNCNT.EQ.3.AND.HNCNT.EQ.0 )THEN
              Ityp = 1
            ELSEIF ( CNCNT.EQ.3.AND.HNCNT.EQ.1 )THEN
              Ityp = 2
            ELSEIF ( CNCNT.EQ.2.AND.HNCNT.EQ.1 )THEN
              Ityp = 2
            ELSE
              Ityp = 3
            ENDIF 
            DO B = 1,4
              DO C=1,NCNT(B)
                DEP = DLST(B,C) 
                IF ( DEP.GE.DMN.and.DEP.LE.DMX) THEN
                  MDCNT(B,Ityp,DEP) = MDCNT(B,Ityp,DEP) + 1
                ENDIF
              ENDDO
            ENDDO
        ENDIF 
      ENDDO

C Print out Bonding depth profile
      OPEN(UNIT=34,FILE='../bonddep.dat',STATUS='unknown')
      DO DP =DMN,DMX
        Drl = REAL(DP) - 0.5d0
        TYPC(1)  = MDCNT(1,1,DP) !3 C-C aromatic
        TYPC(2)  = MDCNT(1,2,DP) !4
        TYPC(3)  = MDCNT(1,3,DP) !5
        TYPC(4)  = MDCNT(2,1,DP) !6
        TYPC(5)  = MDCNT(2,2,DP) !7
        TYPC(6)  = MDCNT(2,3,DP) !8
        TYPC(7)  = MDCNT(3,1,DP) !9   C/C aromatic
        TYPC(8)  = MDCNT(3,2,DP) !10  C/C backbone
        TYPC(9)  = MDCNT(3,3,DP) !11  C/C unkown
        TYPC(10) = MDCNT(4,1,DP) !12
        TYPC(11) = MDCNT(4,2,DP) !13
        TYPC(12) = MDCNT(4,3,DP)
        WRITE(34,340) DEN,Drl*DEV,TYPC(1:12)
      ENDDO
      WRITE(34,*)'' 
      CLOSE(34)
      DEALLOCATE(MDCNT)

C Mass spec analysis on gas phase molecules
C Zero depth count
      ALLOCATE(MDCNT(2,MTPS,DMN:DMX))

      DO T=1,2        !Bonded =1 Unbonded=2 
        DO A=1,MTPS
          DO D=DMN,DMX
            MDCNT(T,A,D) = 0
          ENDDO
        ENDDO
      ENDDO

      IF ( RPREV ) THEN
       WRITE(*,*) 'Read in ../moldep.dat' 
       OPEN(UNIT=35,FILE='../moldep.dat',STATUS='unknown')
       DO M=1,MTCNT
        DO D=DMN,DMX
         READ(35,*) ADUM,BDUM,NMB,NM
         MDCNT(2,M,D) = NM - NMB
         MDCNT(1,M,D) = NMB
        ENDDO
        READ(35,*) 
       ENDDO
       CLOSE(35)
      ENDIF

C Loop over all molecules and output mass 
      WRITE(*,*) MOLCNTF,' molecules where formed'
      DO 1001 NMU=1,MOLCNTF
        Io=MPNTF(NMU)
        Jo=MPNTF(NMU+1)-1
        ECNT = 0
        DO Oo=Io,Jo
          IM=MOLSTF(Oo)
          ECNT = ECNT + 1
          ELIST(ECNT) = IM
        ENDDO
        ALLOCATE (PRAN(ECNT))
        CALL PRNEL(ECNT,ELIST,PREL,PRAN)
        DEALLOCATE (PRAN)
        CALL COMAS(ECNT,ELIST,ELCNT,COM,MOLM)
        Hrl = COM(2)
        Drl = Hrl/DEV+0.5d0   !1/2
        DEP = NINT(Drl)

C Preform molecular mass analysis on Unbonded and gas phase molecules
        DO MT=1,MTCNT  !Check molecular mass agianst recorded molecules
          Mi=INT(MOLM*1000)
          Mf=INT(MMAS(MT)*1000)
          IF ( Mi.EQ.Mf ) THEN
C            IF ( DEP.LE.DMX .and.Hrl.LE.SMN(2,2) ) THEN !Unbonded in substrate
            IF ( DEP.GE.DMN.and.DEP.LE.DMX) THEN
              IF ( DEP.LT.DMN.or.DEP.GT.DMX) THEN
                WRITE(*,*) 'Depth of atom',I,'of molc',MT
                WRITE(*,*) 'At hight ',RFIN(1,I,2),Hrl
                WRITE(*,*) 'With comass=',COM(2)
                WRITE(*,*) 'Which is less than',SMN(2,2)
                WRITE(*,*) 'With index',DEP,'Is out of range',DMN,DMX
                STOP
              ENDIF
              IF ( MNCNT(2,MT).EQ.0 ) THEN
                PREF(1)='sml_'
                PREF(2)='sbx_'
                PBX = .TRUE.
                CALL PRSTR(ECNT,ELIST,MT,PREF,PBX)
              ENDIF
              MNCNT(2,MT) = MNCNT(2,MT) + 1
              MDCNT(2,MT,DEP) = MDCNT(2,MT,DEP) + 1
            ELSE  !Gass phase
              IF ( NMG(MT).EQ.0 ) THEN
                PREF(1)='gml_'
                PREF(2)='gbx_'
                PBX = .FALSE.
                CALL PRSTR(ECNT,ELIST,MT,PREF,PBX)
              ENDIF
              NMG(MT) = NMG(MT) + 1
            ENDIF
            GOTO 1001
          ENDIF
        ENDDO
C If not previously encountered molecule start a new list
        MTCNT = MTCNT  + 1
        IF ( MTCNT.GE.MTPS) THEN
          WRITE(*,*) 'More than ',MTPS,' molc types found'
     &    ,'Increase MTPS'
          STOP
        ENDIF
        MMFRM(MTCNT) = PREL
        MMAS(MTCNT) = MOLM
C        IF ( DEP.LE.DMX .and. Hrl.LE.SMN(2,2) ) THEN !Unbonded in substrate
        IF ( DEP.GE.DMN.and.DEP.LE.DMX) THEN
          IF ( DEP.LT.DMN.or.DEP.GT.DMX) THEN
            WRITE(*,*) 'Depth of atom',I,'of molc',MTCNT
            WRITE(*,*) 'At hight ',RFIN(1,I,2),Hrl
            WRITE(*,*) 'With comass=',COM(2)
            WRITE(*,*) 'Which is less than',SMN(2,2)
            WRITE(*,*) 'With index',DEP,'Is out of range',DMN,DMX
            STOP
          ENDIF
          MNCNT(2,MTCNT) = MNCNT(2,MTCNT) + 1
          MDCNT(2,MTCNT,DEP) = MDCNT(2,MTCNT,DEP) + 1
          PREF(1)='sml_'
          PREF(2)='sbx_'
          PBX = .TRUE.
          CALL PRSTR(ECNT,ELIST,MTCNT,PREF,PBX)
        ELSE  !Gass phase
          NMG(MTCNT) = NMG(MTCNT) + 1
          PREF(1)='gml_'
          PREF(2)='gbx_'
          PBX = .FALSE.
          CALL PRSTR(ECNT,ELIST,MTCNT,PREF,PBX)
        ENDIF
 1001 CONTINUE
        

C Initialize Reacted list
      RCNTB = 0 

C Find which atoms the beam atoms bonded to

      DO RCNT=1,RMX ! Initialize lists
        MRICNT(RCNT)=0
        NNCNT(RCNT)=0
        BRKCNT(RCNT)=0 
      ENDDO
      RCNT = 0

C Print statistics
      WRITE(*,*) 'Printing output'
      WRITE(21,*) 'Deposition of :',MTITLE,'at',DEN
      WRITE(21,*) '' 
      WRITE(21,*) 'H-H cut off :',HHCUT
      WRITE(21,*) 'X-H cut off :',RHCUT
      WRITE(21,*) 'X-X cut off :',RCUT
      WRITE(21,*) "Min and Max of substrate =",SMN(1,2),SMN(2,2)
      WRITE(21,*) 'Buffer =',DBUF
      WRITE(21,*) 'Increase in Potential energy/atom',DPET,' (meV)'
      WRITE(21,4103) NAM*NC
      WRITE(21,4104) BCCNT,FLOAT(BCCNT)/NBEAM
      WRITE(21,4105) NBCNT,FLOAT(NBCNT)/NBEAM
      WRITE(21,4106) BBCNT,FLOAT(BBCNT)/NBEAM
      WRITE(21,*) 'List of xmol.d beam atoms in substrate bonded :'

      DO 1000 K=1,BBCNT
        I = BBLST(K)
        IF ( RECED(I) ) THEN
        ELSE
          BMCNT = 0
          BKCNT = 0 
          GN = 0
          NN = 0
          ERR1=.FALSE.
C Determine if atom I bonded to a substrate atom 
          DO 1101 Nf=1,NCTF(I) !Find new neighbors
            Bf=NLISTF(I,NF)
            GN = GN + 1
            GNTEMP(GN) = Bf  !General Neighbor
            IF ( ITYPE(Bf).EQ.1) GOTO 1101 !can't bond to H of poly
            DO No=1,NCTI(I)
              Bo=NLISTI(I,No)
              DO Po=1,NCTI(Bo) !make sure it wasnt orginaly a 2nd NN
                Ci=NLISTI(Bo,Po)
                IF (Bf.EQ.Ci.OR.Bf.EQ.Bo) GOTO 1101
              ENDDO
            ENDDO
            CALL POLYCHK(I,Bf,
     &        NLISTI,NLISTF,NCTI,NCTF,ASUBI,ASUBF,BFBCK)
            IF ( BFBCK ) THEN
              RECED(I) = .TRUE.
              BKCNT = BKCNT + 1
              BCKBN(BKCNT)=Bf  !Back bone atoms bonded
            ENDIF           
            NN = NN + 1
            NNTEMP(NN) = Bf  !New Neighbor
 1101     CONTINUE
          IF ( RECED(I) ) THEN
            RCNT = RCNT + 1
            RLIST(RCNT) = I
            BMCNT = BMCNT + 1
            BMOL(BMCNT) = I    !Bonded molecule
            WRITE(21,*) 'Bonded Mol #',RCNT
            WRITE(21,*) '  Has a reacted atom ',I
            WRITE(21,*) '  Which bonded to ',BKCNT,'substrate atoms'
            WRITE(21,*) ''
            WRITE(21,*) 'label delete Atoms all'
            WRITE(21,*) 'label delete Bonds all'
            WRITE(21,*) ''
            DO BC=1,BKCNT
              BK=BCKBN(BC)
              WRITE(21,*) 'label add Bonds 0/',I,' 0/',BK
            ENDDO
            DO N=1,NN
              Bf=NNTEMP(N)
              NNCNT(RCNT)= NNCNT(RCNT) + 1 
              NNLSTI(RCNT,NNCNT(RCNT))=I
              NNLSTJ(RCNT,NNCNT(RCNT))=Bf
            ENDDO
            DO 1102 Nf =1,NCTF(I)
             Bf=NLISTF(I,Nf)
             IF ( ITYPE(I).EQ.1.and.ITYPE(Bf).EQ.1 ) GOTO 1102
             DO BC=1,BKCNT
                BK=BCKBN(BC)
                IF ( Bf.EQ.BK ) GOTO 1102
             ENDDO
             RECED(Bf) = .TRUE.
             BMCNT = BMCNT + 1
             BMOL(BMCNT) = Bf 
             IF ( ITYPE(Bf).NE.1) THEN
              DO Mf=1,NCTF(Bf)
               Cf=NLISTF(Bf,Mf)
                IF ( Cf.NE.I ) THEN
                 CALL POLYCHK(Bf,Cf,
     &    NLISTI,NLISTF,NCTI,NCTF,ASUBI,ASUBF,BFBCK)
                 IF ( BFBCK ) THEN
                  BKCNT = BKCNT + 1
                  BCKBN(BKCNT)=Cf  !Back bone atoms bonded
                  WRITE(21,*) 'label add Bonds 0/',Bf-1,' 0/',Cf-1
                 ELSEIF( ITYPE(Cf).EQ.1 ) THEN
                  BMCNT = BMCNT + 1
                  BMOL(BMCNT) = Cf
                  RECED(Cf) = .TRUE.
                 ELSE  
                  BMCNT = BMCNT + 1
                  BMOL(BMCNT) = Cf
                  RECED(Cf) = .TRUE.
                  DO Of=1,NCTF(Cf)
                   Df=NLISTF(Cf,Of)
                   IF ( Df.NE.Bf.and.Df.NE.I) THEN
                    CALL POLYCHK(Cf,Df,
     &    NLISTI,NLISTF,NCTI,NCTF,ASUBI,ASUBF,BFBCK)
                    IF ( BFBCK ) THEN
                     BKCNT = BKCNT + 1
                     BCKBN(BKCNT)=Df  !Back bone atoms bonded
                    WRITE(21,*) 'label add Bonds 0/',Cf-1,' 0/',Df-1
                    ELSEIF ( ITYPE(Df).EQ.1) THEN
                     BMCNT = BMCNT + 1
                     BMOL(BMCNT) = Df
                     RECED(Df) = .TRUE.
                    ELSE
                     BMCNT = BMCNT + 1
                     BMOL(BMCNT) = Df
                     RECED(Df) = .TRUE.
                     DO Pf=1,NCTF(Df)
                      Ef=NLISTF(Df,Pf)
                      IF ( Ef.NE.Cf.and.Ef.NE.Bf.and.Ef.NE.I) 
     & THEN
                       CALL POLYCHK(Df,Ef,
     &    NLISTI,NLISTF,NCTI,NCTF,ASUBI,ASUBF,BFBCK)
                       IF ( BFBCK ) THEN
                        BKCNT = BKCNT + 1
                        BCKBN(BKCNT)=Ef  !Back bone atoms bonded
                    WRITE(21,*) 'label add Bonds 0/',Df-1,' 0/',Ef-1
                        ELSEIF( ITYPE(Ef).EQ.1) THEN
                         BMCNT = BMCNT + 1
                         BMOL(BMCNT) = Ef
                         RECED(Ef) = .TRUE.
                        ELSE
                         BMCNT = BMCNT + 1
                         BMOL(BMCNT) = Ef
                         RECED(Ef) = .TRUE.
                         DO Qf=1,NCTF(Ef)
                          Ff=NLISTF(EF,Qf)
                          IF ( Ff.NE.Df
     & .and.Ff.NE.Cf.and.Ff.NE.Bf.and.Ff.NE.I) THEN
                           CALL POLYCHK(Ef,Ff,
     &    NLISTI,NLISTF,NCTI,NCTF,ASUBI,ASUBF,BFBCK)
                           IF ( BFBCK ) THEN
                            BKCNT = BKCNT + 1
                            BCKBN(BKCNT)=Ff  !Back bone atoms bonded
                    WRITE(21,*) 'label add Bonds 0/',Ef-1,' 0/',Ff-1
                           ELSEIF( ITYPE(Ff).EQ.1 ) THEN
                            BMCNT = BMCNT + 1
                            BMOL(BMCNT) = Ff
                            RECED(Ff) = .TRUE.
                           ELSE
                            BMCNT = BMCNT + 1
                            BMOL(BMCNT) = Ff
                            RECED(Ff) = .TRUE.
                            DO Rf=1,NCTF(Ff)
                             Gf=NLISTF(Ff,Rf)
                             IF ( Gf.NE.Ef.and.Gf.NE.Df
     & .and.Gf.NE.Cf.and.Gf.NE.Bf.and.Gf.NE.I) THEN
                              CALL POLYCHK(Gf,Ff,
     &    NLISTI,NLISTF,NCTI,NCTF,ASUBI,ASUBF,BFBCK)
                              IF ( BFBCK ) THEN
                               BKCNT = BKCNT + 1
                               BCKBN(BKCNT)=Gf  !Back bone atoms bonded
                    WRITE(21,*) 'label add Bonds 0/',Ff-1,' 0/',Ff-1
                              ELSEIF( ITYPE(Gf).EQ.1 ) THEN
                               BMCNT = BMCNT + 1
                               BMOL(BMCNT) = Gf
                               RECED(Gf) = .TRUE.
                              ELSE
                               BMCNT = BMCNT + 1
                               BMOL(BMCNT) = Gf
                               RECED(Gf) = .TRUE.
                               DO Sf=1,NCTF(Gf)
                                Hf=NLISTF(Gf,Sf)
                                IF(Hf.NE.Ef.and.Hf.NE.Df.and.
     & Hf.NE.Ff.and.Hf.NE.Cf.and.Hf.NE.Bf.and.Hf.NE.I) THEN
                                 CALL POLYCHK(Hf,Gf,
     &    NLISTI,NLISTF,NCTI,NCTF,ASUBI,ASUBF,BFBCK)
                                 IF ( BFBCK ) THEN
                                  BKCNT = BKCNT + 1
                                  BCKBN(BKCNT)=Hf  !Back bone atoms bonded
                    WRITE(21,*) 'label add Bonds 0/',Gf-1,' 0/',Hf-1
                                 ELSEIF( ITYPE(Hf).EQ.1 ) THEN
                                  BMCNT = BMCNT + 1
                                  BMOL(BMCNT) = Hf
                                  RECED(Hf) = .TRUE.
                                 ELSE
                                  BMCNT = BMCNT + 1
                                  BMOL(BMCNT) = Hf
                                  RECED(Hf) = .TRUE.
                                  DO Tf=1,NCTF(Hf)
                                   II=NLISTF(Hf,Tf)
                                   IF(II.NE.Ef.and.II.NE.Df.and.
     & II.NE.Gf.and.II.NE.Ff.and.II.NE.Cf.and.II.NE.Bf.and.II.NE.I) THEN
                                    CALL POLYCHK(II,Hf,
     &    NLISTI,NLISTF,NCTI,NCTF,ASUBI,ASUBF,BFBCK)
                                    IF ( BFBCK ) THEN
                                     BKCNT = BKCNT + 1
                                     BCKBN(BKCNT)=II  !Back bone atoms bonded
                    WRITE(21,*) 'label add Bonds 0/',Hf-1,' 0/',II-1
                                    ELSEIF( ITYPE(II).EQ.1 ) THEN
                                     BMCNT = BMCNT + 1
                                     BMOL(BMCNT) = II
                                     RECED(Hf) = .TRUE.
                                    ELSE
                                     BMCNT = BMCNT + 1
                                     BMOL(BMCNT) = II
                                     RECED(II) = .TRUE.
                                     DO Uf=1,NCTF(II)
                                      Jf=NLISTF(II,Uf)
                                      IF(Jf.NE.Ef.and.Jf.NE.Df.and.
     & Jf.NE.Hf.and.Jf.NE.Gf.and.
     & Hf.NE.Ff.and.Hf.NE.Cf.and.Hf.NE.Bf.and.Hf.NE.I) THEN
                                      CALL POLYCHK(Jf,II,
     &    NLISTI,NLISTF,NCTI,NCTF,ASUBI,ASUBF,BFBCK)
                                      IF ( BFBCK ) THEN
                                       BKCNT = BKCNT + 1
                                       BCKBN(BKCNT)=Jf  !Back bone atoms bonded
                    WRITE(21,*) 'label add Bonds 0/',II-1,' 0/',Jf-1
                                      ELSEIF( ITYPE(Jf).EQ.1 ) THEN
                                       BMCNT = BMCNT + 1
                                       BMOL(BMCNT) = Jf
                                       RECED(Jf) = .TRUE.
                                      ELSE
                                       BMCNT = BMCNT + 1
                                       BMOL(BMCNT) = Jf
                                       RECED(Jf) = .TRUE.
                                       DO Vf=1,NCTF(Jf)
                                        Kf=NLISTF(Jf,Vf)
                                        IF(Kf.NE.Ef.and.Kf.NE.Df.and.
     & Kf.NE.Hf.and.Kf.NE.Gf.and.Kf.NE.II.and.
     & Kf.NE.Ff.and.Kf.NE.Cf.and.Kf.NE.Bf.and.Kf.NE.I) THEN
                                         CALL POLYCHK(Kf,Jf,
     &    NLISTI,NLISTF,NCTI,NCTF,ASUBI,ASUBF,BFBCK)
                                         IF ( BFBCK ) THEN
                                          BKCNT = BKCNT + 1
                                          BCKBN(BKCNT)=Kf  !Back bone atoms bonded
                    WRITE(21,*) 'label add Bonds 0/',Jf-1,' 0/',Kf-1
                                         ELSEIF(ITYPE(Kf).EQ.1)THEN
                                          BMCNT = BMCNT + 1
                                          BMOL(BMCNT) = Kf
                                          RECED(Kf) = .TRUE.
                                         ELSE
                                          BMCNT = BMCNT + 1
                                          BMOL(BMCNT) = Kf
                                          RECED(Kf) = .TRUE.
                                          DO Xf=1,NCTF(Kf)
                                           Lf=NLISTF(Kf,Xf)
C Check to see if any more non H nieghbors
       IF(Lf.EQ.1 ) THEN
       ELSEIF(Lf.EQ.Ef.or.Lf.EQ.Df.or.Lf.EQ.Jf.or.
     & Lf.EQ.Hf.or.Lf.EQ.Gf.or.Lf.EQ.II.or.
     & Lf.EQ.Ff.or.Lf.EQ.Cf.or.Lf.EQ.Bf.or.Lf.EQ.I) THEN
       ELSE                                      
         ERR1=.TRUE.
       ENDIF                              
                                          ENDDO
                                         ENDIF
                                        ENDIF
                                       ENDDO
                                      ENDIF
                                     ENDIF
                                    ENDDO
                                   ENDIF
                                  ENDIF
                                 ENDDO
                                ENDIF
                               ENDIF
                              ENDDO
                             ENDIF
                            ENDIF
                           ENDDO
                          ENDIF
                         ENDIF
                        ENDDO
                       ENDIF
                      ENDIF
                     ENDDO
                    ENDIF
                   ENDIF
                  ENDDO
                 ENDIF
                ENDIF
               ENDDO
              ENDIF           
 1102       CONTINUE 
            MCNT = MCNT +  1  !Total # of molecules found
       IF ( ERR1 ) THEN
        WRITE(*,*) 'Warning!! For reaction',RCNT
     &       ,'Molecule not complete'
        WRITE(21,*)  'Atom ',Ff-1,ITYPE(Ff),
     &        'has an unaccounted NN',Gf-1,ITYPE(Gf)
     &       ,'Check output structure for full molecule'
       ENDIF
            ALLOCATE (PRAN(BMCNT))
            CALL PRNEL(BMCNT,BMOL,PREL,PRAN)
            WRITE(21,*) PREL 
            WRITE(21,*) 'mol modselect 2 0 ',PRAN(:)
            DEALLOCATE (PRAN)
            PREF(1)='mol_'
            PREF(2)='rec_'
            PBX = .TRUE. 
            CALL PRSTR(BMCNT,BMOL,RCNT,PREF,PBX)
            WRITE(21,*) BKCNT,'backbone bonds formed'
            WRITE(21,*) ''
            CALL COMAS(BMCNT,BMOL,ELCNT,COM,MOLM)
            Hrl = COM(2)
            Drl = Hrl/DEV+0.5d0   !1/2
            DEP = NINT(Drl)
            IF ( DEP.LT.DMN .or. DEP.GT.DMX) THEN
              WRITE(*,*) 'Depth of moleclue bond to substrate out of
     &        range of depth array. Change array size with DMN,DMX'
              STOP
            ENDIF
C Preform molecular mass analysis on bonded molecules
            DO MT=1,MTCNT  !Check molecular mass agianst recorded molecules
              Mi=INT(MOLM*1000)
              Mf=INT(MMAS(MT)*1000)
              IF ( Mi.EQ.Mf ) THEN
                MNCNT(1,MT) = MNCNT(1,MT) + 1
                MDCNT(1,MT,DEP) = MDCNT(1,MT,DEP) + 1
                GOTO 1000
              ENDIF
            ENDDO
C If not previously encountered molecule start a new list
            MTCNT = MTCNT  + 1   !Total types of molecules found
            IF ( MTCNT.GE.MTPS) THEN
             WRITE(*,*) 'More than ',MTPS,' molc types found'
     &        ,'Increase MTPS'
             STOP
            ENDIF
            MMFRM(MTCNT) = PREL
            MMAS(MTCNT) = MOLM
            MNCNT(1,MTCNT) = MNCNT(1,MTCNT) + 1
            MDCNT(1,MTCNT,DEP) = MDCNT(1,MTCNT,DEP) + 1
C Record depth profile and total crosslinks
            IF ( BKCNT.GT.1 ) THEN
              CRSCNT = CRSCNT  + 1
              CRSD(DEP) = CRSD(DEP) + 1 
            ENDIF
          ENDIF
        ENDIF
 1000 CONTINUE 
      CLOSE(21)
      WRITE(*,*) '  Structures analyized'


C Print out molecule information
      WRITE(27,*) CRSCNT,'crosslinks formed'
      WRITE(27,*) MTCNT,'molecules found'
      WRITE(27,*) ' MaxDepth_C ',MAXDP(1),' MaxDepth_H ',MAXDP(2)
      DPETOT = DPETOT + DPET
      WRITE(27,*) DPETOT,' Change in potential energy (meV) '
      WRITE(27,*) 'Molc #,Formula,Mol Mass,#In sub bonded,
     &  #In sub unbonded,# in gas, Total'
C Loop of different molecule types
      OPEN(UNIT=35,FILE='../moldep.dat',STATUS='unknown')
      DO MT=1,MTCNT
        NMB =  MNCNT(1,MT)      !Number of molc of this type bonded 
        NMU =  MNCNT(2,MT)       !Number of molc of this type unbonded 
        NM =  NMB + NMU + NMG(MT)
        WRITE(27,2701)MT,MMFRM(MT),MMAS(MT),NMB,NMU,NMG(MT),NM
        WRITE(28,2701)MT,MMFRM(MT),MMAS(MT),NMB,NMU,NMG(MT),NM
C  Print out atomic depth profile
        DO DP = DMN,DMX
          Drl=REAL(DP)- 0.5d0
          NMB = MDCNT(1,MT,DP)
          NM = NMB + MDCNT(2,MT,DP)
          WRITE(35,*) MT,Drl*DEV,NMB,NM
        ENDDO
        WRITE(35,*) ''
      ENDDO
      CLOSE(27)
      CLOSE(35)

C Print cross link depth profile
      OPEN(UNIT=33,FILE='../crsdep.dat',STATUS='unknown')
      DO DP = DMN,DMX
        Drl=REAL(DP)- 0.5d0
        WRITE(33,*) DEN,Drl*DEV,CRSD(DP)
      ENDDO
      WRITE(33,*) ''
      CLOSE(33)

      DEALLOCATE (CRSD)
      DEALLOCATE (NMG)
      DEALLOCATE (MDCNT)
      WRITE(*,*) 'Finished'


 4103 FORMAT('Number of beam atoms                     :',I6)
 4104 FORMAT('Number of beam atoms not in substrate     :',I6,':',F7.3)
 4105 FORMAT('Number of beam atoms in substrate unbonded:',I6,':',F7.3)
 4106 FORMAT('Number of beam atoms in substrate bonded  :',I6,':',F7.3)
 4011 FORMAT(A2,I4,I10)
 4012 FORMAT(A2,I4,A2,I4,I10)
 2701 FORMAT(I6,' ',A10,F10.6,' AMU',4I6)
  401 FORMAT(I6,A30,2F8.3,I4)
  340 FORMAT(F6.1,F10.4,12I6)

      STOP

      ENDPROGRAM depanl

      INCLUDE '../sub.rw.f'
      INCLUDE 'sub.molcheckv9.2.f'
      INCLUDE '../sub.nborlstv2.f'
      INCLUDE 'sub.reactmolv4.f'
      INCLUDE 'sub.comas.f'
      INCLUDE 'sub.prnel.f'
      INCLUDE 'sub.anbond.f'
      INCLUDE 'sub.polychk.f'
      INCLUDE 'sub.prstrv2.f'

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

