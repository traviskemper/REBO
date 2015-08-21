C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C     Version 1.0 01/16/2010 T. W. Kemper                      C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C  Module of internal parameters 

      MODULE PARAMS

C Constants
      INTEGER :: DMS,NCC
      parameter ( ncc=10000 ) !Max  of carbon NN
      PARAMETER ( DMS = 3) !# of vector dimensions
      REAL*8 :: PI,BOLZ,EPSI,AVO,ECONV,SIGMA,Atonm,MPatoJ,BOLZeV
      PARAMETER ( PI=3.141592654d0 
     &  ,BOLZ=1.380662d0
     &  ,SIGMA=1.0d0
     &  ,EPSI=11604.50D0  !Check -travisk
     &  ,AVO=6.02205d0
     &  ,ECONV=(1.0D0/(AVO*BOLZ*EPSI))*1.0D+07
     &  ,Atonm = 0.1d0 !Angstroms to nm
     &  ,MPatoJ = 1.0d2      ! J m3 /MPa nm3
     &  ,BOLZeV = 8.6173430d-5
     &  ) 
C     &  ,ECONV= 103.642695082851  ! eV/ ( AMU (A/fs)^2 )
!
      INTEGER :: NPMAX,NSMAX,NLMAX,NPM1,NNMA,NMABIG
      PARAMETER(NPMAX=200000,NSMAX=200000, NLMAX=1000000)
C      PARAMETER(NTAB=10000)
      PARAMETER(NPM1=NPMAX+1)
!
      PARAMETER(NNMA=NPMAX*40)
C      parameter(ikl=2)
      PARAMETER(NMABIG=NNMA*9)
      PARAMETER(VER=4.10D0)
!
      INTEGER :: ktmax
!
!
c tight binding parameters 
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
C change natx to equal number of atoms
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

       parameter  (natx=1,nmax=4*natx,nnx=40)
       parameter  (NElen=256)
       parameter (Am=8.18555d0)
       parameter (rctb=2.18d0)
       parameter (ro=1.536329d0)
       parameter (xnc=6.5d0)
       parameter (xmtb=3.30304d0)
       parameter (xmc=8.655d0)
       parameter (do=1.64d0)
       parameter (dc=2.1052d0)
       parameter (aa1=0.11487638d0)
       parameter (aa2=0.10284431)
       parameter (f0=-2.5909765118191,f1=0.5721151498619d0,
     . f2=-1.7896349903996d-03,f3=2.3539221516757d-05,
     . f4=-1.24251169551587d-07)
c       parameter (stb0=2.2504290109d-08, stb1=-1.4408640561d-06,
c     . stb2=2.1043303374d-05,stb3=6.6024390226d-05)
        parameter (stb0=2.3397813383026395E-08,
     .             stb1=-1.4932260865335458E-06,
     .             stb2=-1.7754111704565561E-04,
     .             stb3=-3.3923114578539729E-03)
       parameter (ttb0=6.7392620074314d-03,ttb1=-8.1885359517898d-02,
     . ttb2=0.1932365259144d0,ttb3=0.354287433238d0)

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
