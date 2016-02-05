!  module for beam builder mkbeam.f 
!
      MODULE beam
!
      INTEGER :: BAMX
      REAL*8 :: ECON,TBUF,BMH,THBUF,BMHT,MBMIN
      PARAMETER( BAMX = 10000 ) !Maximum number of beam atoms
      PARAMETER( ECON = 103.642695082851 ) ! eV/ ( AMU (A/fs)^2 )
      PARAMETER( MBMIN = 0.050d0 )         !min spacing for beam molecules in cfg files
      INTEGER :: NABM,BTYP(BAMX),STPS
      REAL*8 :: RB(BAMX,3),MMS,EB(3),VB(3)
      REAL*8, DIMENSION(3) :: DEPA,SZTHM
      INTEGER :: DEPD
! Molecule information
      INTEGER :: NPM,NC
      CHARACTER(5) :: MTITLE
      INTEGER, ALLOCATABLE :: ITYPEM(:)
      REAL*8,  ALLOCATABLE :: RM(:,:)
      REAL*8 :: MOLMASS
!
      INTEGER :: NAINM,NCLST,MOLMXR
!
      LOGICAL :: MKBM,CENTR,SHFT      
!
      END MODULE beam
