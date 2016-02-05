!
      MDRUN = .TRUE.
      CALC = SmallMolTest 
      SYSTEM = SmallMol 
      TEMP  = 300.0
      SEED = .38562
      LENNARDJONES = .TRUE. 
      THERMOSTATE = -1 ! 
! Parallel options
      SPACIALDECOMP = 2 2 1
! Integrator 
      INTEGRATOR = 1  ! 1- 3rd order PC,2 - Velocity Verlet
! Data output
      STEPS = 1000
      WOUT =  10
      WSTR =  10
      PRINTCFG = .TRUE.
      PRINTXMOL = .TRUE.
! Options 
      CENTER = .FALSE.   !Center structure at origin
      SHIFT = .FALSE.     ! Shift coordinates 
      PRINTPR = .FALSE.  !Print pair information between atom 1 and 2
      PRINTF  = .FALSE.            !Print out forces on each atom 
      PAIR1 = 1
      PAIR2 = 2
      WRITEINITIAL = .FALSE. 
! Thermostate options
      RETHERM = .FALSE.
      NTHERM  = 0
      CUSTOMTHERM = .FALSE.
      BOXTHERM = .FALSE.
      ! also needed for beam build
        RIGIDBOXMIN = 0.d0 0.0d0 0.0d0 
        THERMBOXMIN = 0.d0 20.0d0 0.0d0
!        RIGIDBOXMAX = 0.d0 0.0d0 0.0d0
!        THERMBOXMAX = 0.d0 0.0d0 0.0d0
! Heat 
      HEAT = .FALSE. !Heat structure
      TEMPMAX = 600.0d0
      DELTAT = 50 !K/ps
! Make beam
      MKBEAM = .FALSE.
      DEPDIM    = 3  !z
! Analysis
      MOLANALYSIS = .TRUE.  !Moleculare analysis 
      REACTANALYSIS = .TRUE.  !Reaction analysis
      DEPANALYSIS = .FALSE.  !Deposition analysis
      READREF = .FALSE. !Use referance strucuture as point of comparison
! Clear substrate
      CLEARGAS = .FALSE.

