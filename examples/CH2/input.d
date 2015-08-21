!
      MDRUN = .TRUE.         !Run MD
      TEMP  = 300.000        !System temperature (K)
      SEED = .123456         !Seed for random number generator
      LENNARDJONES = .TRUE. 
      THERMOSTATE = -1       ! KFLAG
! Integrator 
      INTEGRATOR = 1         ! 1- 3rd order PC,2 - Velocity Verlet
! Parallel options
      SPACIALDECOMP = 2 2 1
! Data output
      STEPS = 100              !Number of MD steps
      WOUT = 10               !Steps between data outputs
      WSTR = 10               !Steps between structure outputs
      PRINTCFG = .TRUE.      !Print cfg files for atomeye
      PRINTXMOL = .TRUE.     !Print xmol files for VMD

