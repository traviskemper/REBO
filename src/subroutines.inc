c source code for these subroutines is available in
c routine calling potentials/forces - application specific 
!
      INCLUDE 'open.f'
      INCLUDE 'times.f'
      INCLUDE 'readwrite.f'
      INCLUDE 'nborlist.f'
      INCLUDE 'cutfunc.f'
      INCLUDE 'splines.f'
      INCLUDE 'oxybcuint.f'
      INCLUDE 'pibondCH.f'
      INCLUDE 'pibondCHO.f'
      INCLUDE 'pibondCHS.f'
      INCLUDE 'pibondCHF.f'
      INCLUDE 'heat.f'
      INCLUDE 'default.f'
      INCLUDE 'settherm.f'
      include 'model.f' 
      include 'close.f'
c bond order potential subroutines 
      include 'caguts.f'
      include 'bcuint.f'
      include 'mtable.f'
      include 'param.f'
      include 'radic.f'
      include 'tor.f'
      include 'tcuint.f'
c Lennard-Jones routines 
      include 'ljcont.f'
      include 'ljcset.f'
      include 'ljguts.f'
      include 'ljparam.f'
c general MD routines 
      include 'setin.f' 
      include 'pred_corr.f'
      INCLUDE 'verlet.f'
      include 'thermostats.f'
      include 'random_numbers.f'
c generate data 
      include 'xmol.f'
C Data output
      INCLUDE 'cfg.f'
