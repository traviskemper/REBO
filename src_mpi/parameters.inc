C  NsendMAX: # of atoms send in one time
C  NTAB    : array size for potential table look-up
C  RTYPES  : maximum # of different types of atoms allowed for Lennard-Jones potential
C  natx    : number of atoms used in tight-binding matrices. MUST be set equal
C            to np before compiling.
C
      PARAMETER(NTAB=10000,NSendMAX=5000)
      PARAMETER(RTYPES=8)
      parameter(ikl=2)
      PARAMETER(VER=3.00D0)
c
c tight binding parameters 
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
C change natx to equal number of atoms
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C       parameter  (natx=1,nmax=4*natx,nnx=40)
C       parameter  (NElen=256)
C       parameter (Am=8.18555d0)
C       parameter (rctb=2.18d0)
C       parameter (ro=1.536329d0)C
C       parameter (xnc=6.5d0)
C       parameter (xmtb=3.30304d0)
C       parameter (xmc=8.655d0)
C       parameter (do=1.64d0)
C       parameter (dc=2.1052d0)
C       parameter (aa1=0.11487638d0)
C       parameter (aa2=0.10284431)
C       parameter (f0=-2.5909765118191,f1=0.5721151498619d0,
C     . f2=-1.7896349903996d-03,f3=2.3539221516757d-05,
C     . f4=-1.24251169551587d-07)
c       parameter (stb0=2.2504290109d-08, stb1=-1.4408640561d-06,
c     . stb2=2.1043303374d-05,stb3=6.6024390226d-05)
C        parameter (stb0=2.3397813383026395E-08,
C     .             stb1=-1.4932260865335458E-06,
C     .             stb2=-1.7754111704565561E-04,
C     .             stb3=-3.3923114578539729E-03)
C       parameter (ttb0=6.7392620074314d-03,ttb1=-8.1885359517898d-02,
C     . ttb2=0.1932365259144d0,ttb3=0.354287433238d0)


