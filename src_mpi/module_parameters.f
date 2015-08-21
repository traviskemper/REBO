C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C     Version 1.0 01/16/2010 T. W. Kemper                      C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C  Module of internal parameters 

      MODULE PARAMS

C Constants
      INTEGER :: DMS,NCC,MNDEPV
      parameter ( ncc=10000 ) !Max  of carbon NN
      PARAMETER ( DMS = 3) !# of vector dimensions
      REAL*8 :: PI,BOLZ,EPSI,AVO,ECONV,SIGMA,FCONV
      PARAMETER ( PI=3.141592654d0 
     &  ,BOLZ=1.380662d0
     &  ,SIGMA=1.0d0
     &  ,EPSI=11605.0D0  !Check -travisk
     &  ,AVO=6.02205d0
     &  ,ECONV=(1.0D0/(AVO*BOLZ*EPSI))*1.0D+07
     &  ,FCONV = 6.02205d-8 !AMU / ( MP A fs2)
     &  ) 

C     &  ,ECONV= 103.642695082851  ! eV/ ( AMU (A/fs)^2 )




C      INTEGER :: NPMAX,NSMAX,NLMAX,NPM1,NNMA,NMABIG

C      PARAMETER(NPMAX=80000,NSMAX=100000, NLMAX=190000)

C      PARAMETER(NPM1=NPMAX+1)

C      PARAMETER(NNMA=NPMAX*40)
      INTEGER :: ikl
      parameter(ikl=2)
C      PARAMETER(NMABIG=NNMA*9)
      PARAMETER(VER=4.10D0)

      INTEGER :: ktmax

C MPI parameters
      INTEGER :: NSendMAX   
      PARAMETER(NSendMAX = 5000 )


c
c tight binding parameters 
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
C change natx to equal number of atoms
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$
c$$$       parameter  (natx=1,nmax=4*natx,nnx=40)
c$$$       parameter  (NElen=256)
c$$$       parameter (Am=8.18555d0)
c$$$       parameter (rctb=2.18d0)
c$$$       parameter (ro=1.536329d0)
c$$$       parameter (xnc=6.5d0)
c$$$       parameter (xmtb=3.30304d0)
c$$$       parameter (xmc=8.655d0)
c$$$       parameter (do=1.64d0)
c$$$       parameter (dc=2.1052d0)
c$$$       parameter (aa1=0.11487638d0)
c$$$       parameter (aa2=0.10284431)
c$$$       parameter (f0=-2.5909765118191,f1=0.5721151498619d0,
c$$$     . f2=-1.7896349903996d-03,f3=2.3539221516757d-05,
c$$$     . f4=-1.24251169551587d-07)
c$$$c       parameter (stb0=2.2504290109d-08, stb1=-1.4408640561d-06,
c$$$c     . stb2=2.1043303374d-05,stb3=6.6024390226d-05)
c$$$        parameter (stb0=2.3397813383026395E-08,
c$$$     .             stb1=-1.4932260865335458E-06,
c$$$     .             stb2=-1.7754111704565561E-04,
c$$$     .             stb3=-3.3923114578539729E-03)
c$$$       parameter (ttb0=6.7392620074314d-03,ttb1=-8.1885359517898d-02,
c$$$     . ttb2=0.1932365259144d0,ttb3=0.354287433238d0)

C For spline calculations
      INTEGER isterr,istin,istout
      parameter (isterr = 0)
      parameter (istin  = 5)
      parameter (istout = 6)

c     numbers
      REAL*8 :: bigpos,bigneg
      parameter (bigpos = 9.d98)
      parameter (bigneg = -9.d98)


      END MODULE PARAMS
