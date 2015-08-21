      SUBROUTINE PARAM
C
C
      USE POTS
      USE SPECIF
      USE PARAMS

      IMPLICIT none
!
      INTEGER :: I,J,K,L,M,N,I2D,IC,I3D,ITD,ii,nn,jj,kk
     &,KI,KJ,NTCH,NTO,NTS,NTH,NTC,xmin,xmax,ymin,ymax,zmin,zmax
      REAL*8 :: XHH,y(4),y1(4),y2(4),y12(4),dl(16)
     .     ,f(0:1,0:1,0:1), fx(0:1,0:1,0:1), fy(0:1,0:1,0:1)
     .     ,fz(0:1,0:1,0:1), fxy(0:1,0:1,0:1), fxz(0:1,0:1,0:1)
     .     ,fyz(0:1,0:1,0:1), fxyz(0:1,0:1,0:1),TRICOEF(0:3,0:3,0:3)
     &,dx1,dx2,dy1,dy2,dz1,dz2,dz3,dz4,xnlow,xnhg,ynlow,ynhg,znlow,znhg
     &,RCF,RFF,RHF
      REAL*8 CLM_value(4,16)
!
      IGC(1:16)  = 4
      IGC(17:18) = 3
      IGC(19:20) = 2
      IGC(21:25) = 1

      IGH(1:18)   = 3
      IGH(18:22)  = 2
      IGH(23:25)  = 1

      ATT=3.20D0
      XQM= 3.70D0

      PQ=PI/(XQM-ATT)

      do i=1,RTYPES
          do j=1,RTYPES
               AD(i,j)  = 0.0d0
               AXL(i,j) = 0.0d0
               BD(i,j)  = 0.0d0
               BXL(i,j) = 0.0d0
               CD(i,j)  = 0.0d0
               CXL(i,j) = 0.0d0
               DD(i,j)  = 0.0d0
               DXL(i,j) = 0.0d0
               ED(i,j)  = 0.0d0
               RB1(i,j) = 0.0d0
               RMAX(i,j) = 0.0d0
               PID(i,j) = 1.0d0
               CHI(i,j) = 1.0d0
               do k = 1,4
                 xdb(i,j,k) = 0.0d0
               enddo
          enddo
      enddo
      do i=1,NTYPES
          do j=1,NTYPES
               RB2(i,j) = 0.0d0
          enddo
      enddo
C*** IMPORTANT***************
C                           *
C TO INCLUDE DIHEDRAL TERMS *
C SET NDIHED=2, OTHERWISE   *
C SET NDIHED=10             *
C                           *
C****************************
C
      NDIHED=2
C CARBON
C
      AD(icarb,icarb)=12388.79197798375D0
      AXL(icarb,icarb)=4.720452312717397D0
      BD(icarb,icarb)=17.56740646508968D0
      BXL(icarb,icarb)=1.433213249951261D0
      CD(icarb,icarb)=30.71493208065162D0
      CXL(icarb,icarb)=1.382691250599169D0
      DD(icarb,icarb)=10953.54416216992D0
      DXL(icarb,icarb)=4.746539060659529D0
      ED(icarb,icarb)=0.3134602960832605d0
C Hydrogen
      AD(ihyd,ihyd)=29.6325931D0
      AXL(ihyd,ihyd)=1.715892169856421D0
      BD(ihyd,ihyd)=0.0D0
      BXL(ihyd,ihyd)=1.0D0
      CD(ihyd,ihyd)=0.0D0
      CXL(ihyd,ihyd)=1.0D0
      DD(ihyd,ihyd)=32.81735574722296D0
      DXL(ihyd,ihyd)=3.536298648376465D0
      ED(ihyd,ihyd)=0.3704714870452888d0
c CARBON-HYDROGEN
C
      AD(ihyd,icarb)=32.35518665873256
      AXL(ihyd,icarb)=1.434458059249837
      DD(ihyd,icarb)=149.9409872288120
      DXL(ihyd,icarb)= 4.102549828548784
      ED(ihyd,icarb)=0.3407757282257080
      BD(ihyd,icarb)=0.0D0
      BXL(ihyd,icarb)=1.0D0
      CD(ihyd,icarb)=0.0D0
      CXL(ihyd,icarb)=1.0D0
C Potential cut offs 
C  C-C
      RB1(icarb,icarb)=1.7d0
      RB2(icarb,icarb)=2.0d0

C H-H
      RB1(ihyd,ihyd)=1.10d0
      RB2(ihyd,ihyd)=1.70d0

C C-H
      RB1(ihyd,icarb)=1.3d0
      RB2(ihyd,icarb)=1.8d0

C H-C
      RB1(icarb,ihyd)=RB1(ihyd,icarb)
      RB2(icarb,ihyd)=RB2(ihyd,icarb)
!
! Angular function spline values
!   should be calculated with qucof x
      SPGC(:,1)= (/0.2817216000000E+00, 0.106291200E+01,0.2136736E+01,
     &  0.2533952000000E+01,0.15547360E+01,0.3863296000000E+00/)
      SPGC(:,2)=(/0.28172160E+00,0.1062912000000E+01,0.2136736000E+01,
     &  0.2533952000000E+01,0.1554736000000E+01,0.3863296000000E+00/)
      SPGC(:,3)=(/0.690066866E+00,0.5460691360000E+01,0.23013456800E+02,
     &  0.5491519344000E+02,0.6862037040000E+02,0.3470897779200E+02/)
      SPGC(:,4)=(/0.375449087E+0,.1407252749388E+01,0.2255103926323E+01,
     &  0.2028902219952E+01,0.1426981217906E+01,0.5063107994308E+00/)
      SPGC(:,5)=(/0.2718558E+0,0.4892727456293E+00,-0.4328199017473E+00,
     & -0.5616795197048E+00,0.1270874966906E+01,-0.3750409108350E-01/)
C
      SPGH(:,1)=(/270.467795364007301,1549.701314596994564
     &,3781.927258631323866,4582.337619544424228,2721.538161662818368,
     &630.658598136730774/)
      SPGH(:,2)=(/16.956325544514659,-21.059084522755980,
     &-102.394184748124742,-210.527926707779059,-229.759473570467513,
     &-94.968528666251945/)
      SPGH(:,3)=(/19.065031149937783,2.017732531534021,
     &-2.566444502991983,3.291353893907436,-2.653536801884563,
     &0.837650930130006/)

      XH(:,:,:) =0.0D0
      XH1(:,:,:)=0.0D0
      XH2(:,:,:)=0.0D0
!
! Zero bicubic spline coefficients
!
      DO 125 I=1,2
           DO 124 L=1,10
                DO 123 M=1,10
                     DO 122 J=1,16
                          CLM(I,L,M,J)=0.0
122                  CONTINUE
123             CONTINUE
124        CONTINUE
125   CONTINUE
C
C bicubic spline
C
      READ(15,*) I2D
C
C read integer values of bicubic spline

C
  119 CONTINUE
      READ(15,*) I,J,K,XHH
      IF(I.LE.0) GO TO 121
      XH(I,J,K)=XHH
      GO TO 119
  121 CONTINUE 
C
C
      XH1(2,3,1)=(XH(2,4,1)-XH(2,2,1))/2.0D0
      XH1(2,2,2)=(XH(2,3,2)-XH(2,1,2))/2.0D0
C
      XH2(2,2,2)=(XH(2,2,3)-XH(2,2,1))/2.0D0
      XH2(2,1,3)=(XH(2,1,4)-XH(2,1,2))/2.0D0
C
C Read bicubic spline coefficients
C
  126 CONTINUE
      READ(15,*,END=127) I,L,M
      READ(15,*) (CLM(I,L,M,J),J=1,16)
      GO TO 126
  127 CONTINUE
C
      IC=0
C
      DO 9 I=1,4
           DO 8 J=1,4
                     IC=IC+1
                     IN2(IC,1)=I-1
                     IN2(IC,2)=J-1
8               CONTINUE
9         CONTINUE
C
! Read in tricubic spline coefficients for CHF
C Zero tricubic spline coefficients 
C
      DO  I=1,3
        DO  L=1,10
          DO  M=1,10
            DO  N=1,10 
              XHF(I,L,M,N) = 0.0d0
              DO  J=1,64
                CLMF(I,L,M,N,J)=0.0d0
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      READ(19,*) I2D
219   CONTINUE
      READ(19,*) I,J,K,l,XHH
      IF(I.LE.0) GO TO 221
      IF( I .EQ. 3 ) I = iflor
      XHF(I,J,K,l)=XHH
      GO TO 219
221   CONTINUE

! begin debug
c$$$      DO  I=1,4
c$$$        DO  L=1,4
c$$$          DO  M=1,4
c$$$            DO  N=1,4
c$$$              XHH =XHF(I,L,M,N)
c$$$              IF(XHH.NE.0d0)WRITE(77,771) I,L,M,N,XHH
c$$$            ENDDO
c$$$          ENDDO
c$$$        ENDDO
c$$$      ENDDO
c$$$ 771  FORMAT('      XHF(',3(I6,','),I6,')=',F16.8,'d0')
! end debug

!
      DO nn=1,3*6*6*6
        READ(19,*) I,L,M,N
        DO j=1,16
             READ(19,*) CLM_value(1,j),CLM_value(2,j),
     &       CLM_value(3,j),CLM_value(4,j)
        ENDDO
        kk=0
        ! Change default atom numbers to correct values 
        IF( I .EQ. 3 ) I = iflor
        DO ii=1,16
          DO jj=1,4
            kk=kk+1
            CLMF(I,L,M,N,kk)=CLM_value(jj,ii)
          ENDDO
        ENDDO
      ENDDO
C
C Read tricubic spline coefficients
C
      IC=0
C
      DO 7 I=1,4
           DO 6 J=1,4
                DO 5 K=1,4
                     IC=IC+1
                     IN3(IC,1)=I-1
                     IN3(IC,2)=J-1
                     IN3(IC,3)=K-1
5               CONTINUE
6         CONTINUE
7     CONTINUE
C
      READ(14,*) I3D
  129 CONTINUE
      READ(14,*,END=130) L,M,N
      READ(14,*) (CLMN(1,L,M,N,I),I=1,64)
      GO TO 129
  130 CONTINUE
C
      READ(17,*) I3D
  229 CONTINUE


      READ(17,*,END=230) L,M,N
      READ(17,*) (CLMN(3,L,M,N,I),I=1,64)
      GO TO 229
  230 CONTINUE
C
      READ(18,*) I3D
  239 CONTINUE
      READ(18,*,END=240) L,M,N
      READ(18,*) (CLMN(2,L,M,N,I),I=1,64)
      GO TO 239
  240 CONTINUE
C
      DO 134 L=1,10
           DO 133 M=1,10
                DO 132 I=1,64
                     CLMN(1,L,M,10,I)=CLMN(1,L,M,9,I)
                     CLMN(2,L,M,10,I)=CLMN(2,L,M,9,I)
                     DO 131 N=6,10
                          CLMN(3,L,M,N,I)=CLMN(3,L,M,5,I)
131                  CONTINUE
132             CONTINUE
133        CONTINUE
134   CONTINUE
C
C Read tricubic spline coefficients for torsional potential
C
      READ(16,*) ITD
  135 CONTINUE
      READ(16,*,END=136) L,M,N
      READ(16,*) (TLMN(L,M,N,I),I=1,64)
      GO TO 135
  136 CONTINUE
C
      DO 140 L=1,10
           DO 139 M=1,10
                DO 138 N=4,10
                     DO 137 I=1,64
                          TLMN(L,M,N,I)=TLMN(L,M,3,I)
137                  CONTINUE
138             CONTINUE
139        CONTINUE
140   CONTINUE
C

      IF((ITD.NE.I2D).OR.(ITD.NE.I3D)) THEN
            WRITE(6,*) 'INCOMPATABLE POTENTIAL TYPES'
            WRITE(6,*) ITD,I2D,I3D
            include 'close.inc'
            STOP
      ENDIF
!
      IC=0
      DO I=1,4
        DO J=1,4
          IC=IC+1
          IN2(IC,1)=I-1
          IN2(IC,2)=J-1
        ENDDO
      ENDDO
! 
! Oxygen
      AD(ioxygen,ioxygen)=1105.00D0
      AXL(ioxygen,ioxygen)=1.325D0
      BD(ioxygen,ioxygen)=0.0D0
      BXL(ioxygen,ioxygen)=1.0D0
      CD(ioxygen,ioxygen)=0.0D0
      CXL(ioxygen,ioxygen)=1.0D0
      DD(ioxygen,ioxygen)=685.2555D0
      DXL(ioxygen,ioxygen)=1.173D0
      ED(ioxygen,ioxygen)=0.4065d0
!carbon-oxygen
      AD(ioxygen,icarb)=268.043
      AXL(ioxygen,icarb)=2.344
      DD(ioxygen,icarb)=81.0576
      DXL(ioxygen,icarb)=3.554
      ED(ioxygen,icarb)=9.13200
c
      BD(ioxygen,icarb)=0.0D0
      BXL(ioxygen,icarb)=1.0D0
      CD(ioxygen,icarb)=0.0D0
      CXL(ioxygen,icarb)=1.0D0
!hydogen-oxygen
      AD(ioxygen,ihyd)=884.5045
      AXL(ioxygen,ihyd)=1.704
      DD(ioxygen,ihyd)=717.14950
      DXL(ioxygen,ihyd)= 1.628
      ED(ioxygen,ihyd)=0.12400
c
      BD(ioxygen,ihyd)=0.0D0
      BXL(ioxygen,ihyd)=1.0D0
      CD(ioxygen,ihyd)=0.0D0
      CXL(ioxygen,ihyd)=1.0D0
C O-O
      RB1(ioxygen,ioxygen) = 1.55d0
      RB2(ioxygen,ioxygen) = 1.8d0
C C-O
      RB1(icarb,ioxygen) = 1.6d0
      RB2(icarb,ioxygen) = 1.9d0
C O-H
      RB1(ihyd,ioxygen) = 1.2d0
      RB2(ihyd,ioxygen) = 1.7d0
C Oxygen angular parameters
      a_o0 = -0.014d0
      a_o1 = 0.07d0
      a_o2 = -0.478d0
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   Set bicubic spline coefficients for P function
C   P (atomi,atomj,m+1,n+1) 
C   From table 4 : Boris Ni J. Phys. Condens. Matter 16 (2004) 7261-7275
C     Travis

C Initialize all spline values to 1 inorder to penalized overcoordination
      DO KI=1,RTYPES
        DO KJ=1,RTYPES
          DO NTCH=1,9
            DO NTO=1,9
              XHO(KI,KJ,NTCH,NTO) = 0.0d0
              XHO1(KI,KJ,NTCH,NTO) = 0.0d0
              XHO2(KI,KJ,NTCH,NTO) = 0.0d0
              XHO12(KI,KJ,NTCH,NTO) = 0.0d0
            ENDDO
          ENDDO
        ENDDO
      ENDDO
        
C   Set correct coordinations to zero

      XHO( 1,1,1,1) =   0.000000000000000D000
      XHO( 1,1,2,1 ) = 0.000000000000000D000
      XHO( 1,1,3,1) =  0.000000000000000D000
      XHO( 1,1,4,1) =  0.000000000000000D000
      XHO( 1,1,1,2) =  0.1390000000000D00     
      XHO( 1,1,2,2) =  0.3500000000000D00     
      XHO( 1,1,3,2) = -2.000000000000000D-2
      XHO( 1,1,1,3) =  0.2310D0
      XHO( 1,1,2,3) = -7.661999999999999E-005
      XHO( 1,1,1,4) = -6.2530D-3
      XHO( 1,2,1,1) =  0.0D0
      XHO( 1,2,2,1) =  0.0D0
      XHO( 1,2,3,1) =  0.0D0
      XHO( 1,2,4,1) =  0.00D0
      XHO( 1,2,1,2) =  0.3850D0
      XHO( 1,2,2,2) =  0.137D0
      XHO( 1,2,3,2) = -0.229D0
      XHO( 1,2,1,3) = -8.300000000000000D-002
      XHO( 1,2,2,3) = -0.238000000000000     
      XHO( 1,2,1,4) = -0.224000000000000     
      XHO( 1,2,4,2) =   1.00000000000000     
      XHO( 1,ioxygen,1,1) = -0.3900000D0
      XHO( 1,ioxygen,2,1) =  0.128000000000000     
      XHO( 1,ioxygen,3,1) =  6.600000000000000D-002
      XHO( 1,ioxygen,4,1) =  0.100000000000000     
      XHO( 1,ioxygen,1,2) =  0.109000000000000     
      XHO( 1,ioxygen,2,2) =  0.173000000000000     
      XHO( 1,ioxygen,3,2) =  0.135000000000000     
      XHO( 1,ioxygen,1,3) =  0.177000000000000     
      XHO( 1,ioxygen,2,3) =  9.689000000000000D-005
      XHO( 1,ioxygen,1,4) =  5.400000000000000D-002
      XHO( 1,ioxygen,4,2) =   12.5000D0
      XHO( ioxygen,1,1,1) = -0.288000D0
      XHO( ioxygen,1,2,1) =  0.607000D0
      XHO( ioxygen,1,3,1) =   19.0570D0
      XHO( ioxygen,1,4,1) =   19.0570D0
      XHO( ioxygen,1,1,2) =   1.02600D0
      XHO( ioxygen,1,2,2) =   19.0570D0
      XHO( ioxygen,1,3,2) =   19.0570D0
      XHO( ioxygen,1,1,3) =   19.0570D0
      XHO( ioxygen,1,2,3) =   19.0570D0
      XHO( ioxygen,1,1,4) =   19.0570D0
      XHO( ioxygen,2,1,1) = -2.200000000000000D-002
      XHO( ioxygen,2,2,1) = -1.100000000000000D-002
      XHO( ioxygen,2,3,1) =  7.500000000000000D-002
      XHO( ioxygen,2,4,1) =  8.200000000000000D-002
      XHO( ioxygen,2,1,2) = -6.560000000000000D-003
      XHO( ioxygen,2,2,2) =  7.500000000000000D-002
      XHO( ioxygen,2,3,2) =  8.200000000000000D-002
      XHO( ioxygen,2,1,3) =  7.500000000000000D-002
      XHO( ioxygen,2,2,3) =  8.200000000000000D-002
      XHO( ioxygen,2,1,4) =  8.200000000000000D-002
      XHO( 2,ioxygen,1,1) = -1.900000000000000D-002  !Possibe error should be (3,1,2,1)
      XHO( ioxygen,ioxygen,1,1) = -3.600000000000000D-002
      XHO( ioxygen,ioxygen,2,1) =  1.166000000000000D-003
      XHO( ioxygen,ioxygen,3,1) =  6.200000000000000D-002
      XHO( ioxygen,ioxygen,4,1) =  7.099999999999999D-002
      XHO( ioxygen,ioxygen,1,2) =  2.800000000000000D-002
      XHO( ioxygen,ioxygen,2,2) =  7.099999999999999D-002
      XHO( ioxygen,ioxygen,3,2) =  7.099999999999999D-002
      XHO( ioxygen,ioxygen,1,3) =  6.200000000000000D-002
      XHO( ioxygen,ioxygen,2,3) =  7.099999999999999D-002
      XHO( ioxygen,ioxygen,1,4) =  7.099999999999999D-002
C Overwrite previous spline values with refitted to dissociation energies
       XHO(icarb,icarb, 1 , 2 )= 2.6322758289717174E-001   !  
       XHO(icarb,icarb, 1 , 3 )= 1.8172173808666534E-001   !  
       XHO(icarb,icarb, 1 , 4 )= -3.3812822082183042E-001   !  
       XHO(icarb,icarb, 2 , 2 )= 6.8975879892882921E-003   !  
       XHO(icarb,icarb, 2 , 3 )= -1.4590703313216463E-001   !  
       XHO(icarb,icarb, 3 , 2 )= -1.8833795634238631E-001   !  
       XHO(icarb,ihyd, 1 , 2 )= 3.3634981543859399E-001   !  
       XHO(icarb,ihyd, 1 , 3 )= 6.3516806142203908E-001   !  
       XHO(icarb,ihyd, 1 , 4 )= -5.0686982615829446E-001   !  
       XHO(icarb,ihyd, 2 , 2 )= 1.9144348570035366E-001   !  
       XHO(icarb,ihyd, 2 , 3 )= -1.4434919479893391E-001   !  
       XHO(icarb,ihyd, 3 , 2 )= -3.5365613785081512E-001   !  
       XHO(icarb,ioxygen, 1 , 1 )= -4.9959316118488212E-001   !  
       XHO(icarb,ioxygen, 1 , 2 )= -1.3833764471605306E-001   !  
       XHO(icarb,ioxygen, 1 , 3 )= -1.1963189424243577E-001   !  
       XHO(icarb,ioxygen, 1 , 4 )= -1.9289035017487458E-001   !  
       XHO(icarb,ioxygen, 2 , 1 )= -2.0681692398078513E-001   !  
       XHO(icarb,ioxygen, 2 , 2 )= -2.4451392840821728E-001   !  
       XHO(icarb,ioxygen, 2 , 3 )= -7.5386271861458604E-002   !  
       XHO(icarb,ioxygen, 3 , 1 )= -1.8530763275453846E-001   !  
       XHO(icarb,ioxygen, 3 , 2 )= -2.0881444501255306E-001   !  
       XHO(icarb,ioxygen, 4 , 1 )= 5.8516357773566230E-001   !  
       XHO(ioxygen,icarb, 1 , 1 )= -1.0215691047467139E-001   !  
       XHO(ioxygen,icarb, 1 , 2 )= 5.7125255608847153E-002   !  
       XHO(ioxygen,icarb, 2 , 1 )= 3.7384717869475698E-001   !  
       XHO(ioxygen,ihyd, 1 , 1 )= -2.5469101731520657E-002   !  
       XHO(ioxygen,ihyd, 1 , 2 )= -2.2104792958594469E-002   !  
       XHO(ioxygen,ihyd, 2 , 1 )= -1.2921029783776575E-002   !  
       XHO(ioxygen,ioxygen, 1 , 1 )= -3.1882500229465302E-002   !  
       XHO(ioxygen,ioxygen, 1 , 2 )= 1.9665671695954869E-002   !  
       XHO(ioxygen,ioxygen, 2 , 1 )= -3.4852252855995043E-003   !  
C Calculate bicubic spline derivatives
C
      XHO1(1,1,2,2)=(XHO(1,1,3,2)-XHO(1,1,1,2))/2.0D0
      XHO1(1,1,2,1)=(XHO(1,1,3,1)-XHO(1,1,1,1))/2.0D0
      XHO1(1,1,3,1)=(XHO(1,1,4,1)-XHO(1,1,2,1))/2.0D0
C
      XHO2(1,1,2,2)=(XHO(1,1,2,3)-XHO(1,1,2,1))/2.0D0
      XHO2(1,1,1,2)=(XHO(1,1,1,3)-XHO(1,1,1,1))/2.0D0
      XHO2(1,1,1,3)=(XHO(1,1,1,4)-XHO(1,1,1,2))/2.0D0
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      XHO1(1,2,2,2)=(XHO(1,2,3,2)-XHO(1,2,1,2))/2.0D0
      XHO1(1,2,2,1)=(XHO(1,2,3,1)-XHO(1,2,1,1))/2.0D0
      XHO1(1,2,3,1)=(XHO(1,2,4,1)-XHO(1,2,2,1))/2.0D0
CCC            XHO1(1,2,3,2)=(XHO(1,2,4,2)-XHO(1,2,2,2))/2.0D0
C
      XHO2(1,2,2,2)=(XHO(1,2,2,3)-XHO(1,2,2,1))/2.0D0
      XHO2(1,2,1,2)=(XHO(1,2,1,3)-XHO(1,2,1,1))/2.0D0
      XHO2(1,2,1,3)=(XHO(1,2,1,4)-XHO(1,2,1,2))/2.0D0
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      XHO1(1,ioxygen,2,2)=(XHO(1,ioxygen,3,2)-XHO(1,ioxygen,1,2))/2.0D0
      XHO1(1,ioxygen,2,1)=(XHO(1,ioxygen,3,1)-XHO(1,ioxygen,1,1))/2.0D0
      XHO1(1,ioxygen,3,1)=(XHO(1,ioxygen,4,1)-XHO(1,ioxygen,2,1))/2.0D0
C
      XHO2(1,ioxygen,2,2)=(XHO(1,ioxygen,2,3)-XHO(1,ioxygen,2,1))/2.0D0
      XHO2(1,ioxygen,1,2)=(XHO(1,ioxygen,1,3)-XHO(1,ioxygen,1,1))/2.0D0
      XHO2(1,ioxygen,1,3)=(XHO(1,ioxygen,1,4)-XHO(1,ioxygen,1,2))/2.0D0

CC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
CCC          XHO1(3,1,1,1)=-15.0
CCC          XHO1(3,1,2,1)= 0.0
CCC          XHO1(3,1,1,2)= 0.0 

      XHO1(ioxygen,1,2,2)=(XHO(ioxygen,1,3,2)-XHO(ioxygen,1,1,2))/2.0D0
CCC      XHO1(ioxygen,1,2,1)=(XHO(ioxygen,1,3,1)-XHO(ioxygen,1,1,1))/2.0D0
      XHO1(ioxygen,1,3,1)=(XHO(ioxygen,1,4,1)-XHO(ioxygen,1,2,1))/2.0D0
C
ccc      XHO2(3,1,2,2)=(XHO(3,1,2,3)-XHO(3,1,2,1))/2.0D0
ccc      XHO2(3,1,1,2)=(XHO(3,1,1,3)-XHO(3,1,1,1))/2.0D0
ccc      XHO2(3,1,1,3)=(XHO(3,1,1,4)-XHO(3,1,1,2))/2.0D0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccC
c
      XHO1(ioxygen,2,2,2)=(XHO(ioxygen,2,3,2)-XHO(ioxygen,2,1,2))/2.0D0
      XHO1(ioxygen,2,2,1)=(XHO(ioxygen,2,3,1)-XHO(ioxygen,2,1,1))/2.0D0
      XHO1(ioxygen,2,3,1)=(XHO(ioxygen,2,4,1)-XHO(ioxygen,2,2,1))/2.0D0
C
      XHO2(ioxygen,2,2,2)=(XHO(ioxygen,2,2,3)-XHO(ioxygen,2,2,1))/2.0D0
      XHO2(ioxygen,2,1,2)=(XHO(ioxygen,2,1,3)-XHO(ioxygen,2,1,1))/2.0D0
      XHO2(ioxygen,2,1,3)=(XHO(ioxygen,2,1,4)-XHO(ioxygen,2,1,2))/2.0D0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccC
c
      XHO1(ioxygen,ioxygen,2,2)=(XHO(ioxygen,ioxygen,3,2)
     & -XHO(ioxygen,ioxygen,1,2))/2.0D0
      XHO1(ioxygen,ioxygen,2,1)=(XHO(ioxygen,ioxygen,3,1)
     & -XHO(ioxygen,ioxygen,1,1))/2.0D0
      XHO1(ioxygen,ioxygen,3,1)=(XHO(ioxygen,ioxygen,4,1)
     & -XHO(ioxygen,ioxygen,2,1))/2.0D0
C
      XHO2(ioxygen,ioxygen,2,2)=(XHO(ioxygen,ioxygen,2,3)
     & -XHO(ioxygen,ioxygen,2,1))/2.0D0
      XHO2(ioxygen,ioxygen,1,2)=(XHO(ioxygen,ioxygen,1,3)
     & -XHO(ioxygen,ioxygen,1,1))/2.0D0
      XHO2(ioxygen,ioxygen,1,3)=(XHO(ioxygen,ioxygen,1,4)
     & -XHO(ioxygen,ioxygen,1,2))/2.0D0

C  Calculate spline coefficients
C  Calculate spline coefficients
C   -travisk
C
      DO KI=1,RTYPES
        DO KJ=1,RTYPES
          DO NTCH=1,9
            DO NTO=1,9
                  xmin = NTCH
                  xmax = NTCH + 1
                  ymin = NTO
                  ymax = NTO + 1

                  y(1) = XHO(KI,KJ,NTCH,NTO)
                  y(2) = XHO(KI,KJ,NTCH+1,NTO)
                  y(3) = XHO(KI,KJ,NTCH+1,NTO+1)
                  y(4) = XHO(KI,KJ,NTCH,NTO+1)

                  y1(1) = XHO1(KI,KJ,NTCH,NTO)
                  y1(2) = XHO1(KI,KJ,NTCH+1,NTO)
                  y1(3) = XHO1(KI,KJ,NTCH+1,NTO+1)
                  y1(4) = XHO1(KI,KJ,NTCH,NTO+1)

                  y2(1) = XHO2(KI,KJ,NTCH,NTO)
                  y2(2) = XHO2(KI,KJ,NTCH+1,NTO)
                  y2(3) = XHO2(KI,KJ,NTCH+1,NTO+1)
                  y2(4) = XHO2(KI,KJ,NTCH,NTO+1)

                  y12(1) = 0.d0
                  y12(2) = 0.d0
                  y12(3) = 0.d0
                  y12(4) = 0.d0

                  call bcucof(xmin, xmax, ymin, ymax, y, y1, y2, y12,dl)
                  DO i =1,16  
                    CLMOX(KI,KJ,NTCH,NTO,i) = dl(i)
                  ENDDO
              ENDDO
            ENDDO   
         ENDDO   
      ENDDO    

C Set power factors



!
!   SULFUR-SULFUR
      ED(isulfur,isulfur)= 1.5818786797564960d-001 
      DD(isulfur,isulfur)= 4.4891090501082645d+003 
      DXL(isulfur,isulfur)= 1.7861952468560334d+000 
      AD(isulfur,isulfur)= 4.9060661549429469d+003 
      AXL(isulfur,isulfur)= 1.7771197129506962d+000 

      BD(isulfur,isulfur)=0.0D0
      BXL(isulfur,isulfur)=1.0D0
      CD(isulfur,isulfur)=0.0D0
      CXL(isulfur,isulfur)=1.0D0
C   SULFUR-CARBON
      ED(isulfur,icarb)= 9.5352969102182772d-001 
      DD(isulfur,icarb)= 7.5585813003331725d+002 
      DXL(isulfur,icarb)= 1.9034274245601037d+000 
      AD(isulfur,icarb)= 1.4237382337244990d+003 
      AXL(isulfur,icarb)= 1.9383112114776908d+000 

      BD(isulfur,icarb)=0.0D0
      BXL(isulfur,icarb)=1.0D0
      CD(isulfur,icarb)=0.0D0
      CXL(isulfur,icarb)=1.0D0
 
C   SULFUR-HYDROGEN
      ED(isulfur,ihyd)=0.08697809d0
      DD(isulfur,ihyd)=1033.7328d0
      DXL(isulfur,ihyd)=0.80874612d0
      AD(isulfur,ihyd)=1167.2251d0
      AXL(isulfur,ihyd)=0.84514851d0

      BD(isulfur,ihyd)=0.0D0
      BXL(isulfur,ihyd)=1.0D0
      CD(isulfur,ihyd)=0.0D0
      CXL(isulfur,ihyd)=1.0D0
C S-S
      RB1(isulfur,isulfur) = 2.3d0
      RB2(isulfur,isulfur) = 2.6d0
C C-S
      RB1(icarb,isulfur) = 2.1d0
      RB2(icarb,isulfur) = 2.4d0
C S-H
      RB1(ihyd,isulfur) = 1.5d0
      RB2(ihyd,isulfur) = 1.8d0 
C S8 - anglular
       SPGS(isulfur, 1 )= 1.8751407074297566E-002 
       SPGS(isulfur, 2 )= 6.8577808799679826E-002 
       SPGS(isulfur, 3 )= 8.0687840445684694E-002 
       SPGS(isulfur, 4 )= -7.6424998991209775E-002 
       SPGS(isulfur, 5 )= -5.4377598661551257E-002 
       SPGS(isulfur, 6 )= 2.5375402533424827E-002 
       SPGS(isulfur, 7 )= 1.4768843195291413E-002 
       SPGS(ihyd, 1 )= 2.9472201695595075E-003 
       SPGS(ihyd, 2 )= 1.5955601114131474E-003 
       SPGS(ihyd, 3 )= 1.8274954548820197E-002 
       SPGS(ihyd, 4 )= 2.7447754723426216E-002 
       SPGS(ihyd, 5 )= 1.1384795826984137E-002 
       SPGS(ihyd, 6 )= -1.9797112425817422E-002 
       SPGS(ihyd, 7 )= -1.1964385623099929E-002 
!       SPGS(icarb, 1 )= 6.9834457599274716E-002 
!       SPGS(icarb, 2 )= 7.3311352125000867E-002 
!       SPGS(icarb, 3 )= 2.7473838391582317E-001 
!       SPGS(icarb, 4 )= 2.1558231533511055E-001 
!       SPGS(icarb, 5 )= -5.2438060922228649E-002 
!       SPGS(icarb, 6 )= -1.8804941429676142E-001 
!       SPGS(icarb, 7 )= -2.8428333140184517E-002 
! weaken X-S-C angle for thiirane  
       SPGS(icarb, 1 )= 6.9657382304200835E-002 
       SPGS(icarb, 2 )= 5.7953774115340036E-002 
       SPGS(icarb, 3 )= 2.1053917695323954E-001 
       SPGS(icarb, 4 )= 1.3771562314358143E-001 
       SPGS(icarb, 5 )= -6.6320937972030744E-002 
       SPGS(icarb, 6 )= -1.3076853315410142E-001 
       SPGS(icarb, 7 )= -2.0138279792248997E-003 
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   Set bicubic spline coefficients for P functiona for sulfur
C   P (atomi,atomj,m+1,n+1) 
C    -travisk 
C Initialize all spline values to 0 
      DO KI=1,RTYPES
        DO KJ=1,RTYPES
          DO NTC=0,9
           DO NTH=0,9
            DO NTS=0,9
              XHS(KI,KJ,NTC,NTH,NTS) = 0.0d0
              XHS1(KI,KJ,NTC,NTH,NTS)= 0.0d0
              XHS2(KI,KJ,NTC,NTH,NTS)= 0.0d0
              XHS3(KI,KJ,NTC,NTH,NTS)= 0.0d0
              XHS12(KI,KJ,NTC,NTH,NTS)= 0.0d0
              XHS13(KI,KJ,NTC,NTH,NTS)= 0.0d0
              XHS23(KI,KJ,NTC,NTH,NTS)= 0.0d0
              XHS123(KI,KJ,NTC,NTH,NTS)= 0.0d0
              DO i =1,64
                CLMSX(i,NTC,NTH,NTS,KI,KJ) = 0.0d0
              ENDDO
             ENDDO    
            ENDDO   
          ENDDO    
        ENDDO   
      ENDDO    
! Original REBO II values for non sulfur containing compounds 
      XHS(icarb,icarb,0,1,0) = 0.0000000000000000D+00
      XHS(icarb,icarb,1,1,0) = 3.0266974734813094D-03
      XHS(icarb,icarb,0,2,0) = 7.8607002547458685D-03
      XHS(icarb,icarb,0,3,0) = 1.6125364564267377D-02
      XHS(icarb,icarb,2,1,0) = 3.1795308307312573D-03
      XHS(icarb,icarb,1,2,0) = 6.3262482411195846D-03
      XHS(icarb,ihyd,0,0,0)   =  0.0000000000000000D+00
      XHS(icarb,ihyd,0,1,0)   = 0.2093367328250380D0  
      XHS(icarb,ihyd,0,2,0)   =-6.4449615432525347D-02
      XHS(icarb,ihyd,0,3,0)   = -0.3039275463461620D0     
      XHS(icarb,ihyd,1,0,0)   = 1.0000000000000000D-02
      XHS(icarb,ihyd,2,0,0)   = -0.1220421462782555D0    
      XHS(icarb,ihyd,1,1,0)   =  -0.1251234006287090D0    
      XHS(icarb,ihyd,1,2,0)   = -0.2989052457826418D0    
      XHS(icarb,ihyd,3,0,0)   = -0.3075847050665519D0  
      XHS(icarb,ihyd,2,1,0)   = -0.3005291724067579D0      
! From fitv19 fit
       XHS( 1 , 1 , 0 , 0 ,  1 )= -1.2789986828116587E-002   !  
       XHS( 1 , 1 , 0 , 0 ,  2 )= 1.2103636547684291E-001   !  
       XHS( 1 , 1 , 0 , 0 ,  3 )= -1.6793556091300475E-002   !  
       XHS( 1 , 1 , 0 , 1 ,  1 )= 6.7123449081973291E-002   !  
       XHS( 1 , 1 , 0 , 1 ,  2 )= -7.7046126624899861E-002   !  
       XHS( 1 , 1 , 0 , 2 ,  1 )= 1.2601368234600176E-002   !  
       XHS( 1 , 1 , 1 , 0 ,  1 )= 9.5020583372489489E-002   !  
       XHS( 1 , 1 , 1 , 0 ,  2 )= -1.6501737349813836E-001   !  
       XHS( 1 , 1 , 1 , 1 ,  1 )= -1.2817255312646481E-001   !  
       XHS( 1 , 1 , 2 , 0 ,  1 )= -9.2699822134416687E-002   !  
       XHS( 1 , 2 , 0 , 0 ,  1 )= 9.3977124833385850E-001   !  
       XHS( 1 , 2 , 0 , 0 ,  2 )= 4.4642061591068544E-001   !  
       XHS( 1 , 2 , 0 , 0 ,  3 )= -2.2901058789241233E-001   !  
       XHS( 1 , 2 , 0 , 1 ,  1 )= 2.8986125274688296E-001   !  
       XHS( 1 , 2 , 0 , 1 ,  2 )= -4.4193864932197552E-001   !  
       XHS( 1 , 2 , 0 , 2 ,  1 )= -2.8520731512159192E-001   !  
       XHS( 1 , 2 , 1 , 0 ,  1 )= 4.2793940654159995E-001   !  
       XHS( 1 , 2 , 1 , 0 ,  2 )= -3.0323920425182616E-001   !  
       XHS( 1 , 2 , 1 , 1 ,  1 )= -3.8640549472107072E-001   !  
       XHS( 1 , 2 , 2 , 0 ,  1 )= -2.1558431558861113E-001   !  
       XHS(icarb,isulfur, 0 , 0 ,  0 )= 3.4118902134013784E-001   !  
       XHS(icarb,isulfur, 0 , 0 ,  1 )= 4.8693356702919893E-001   !  
       XHS(icarb,isulfur, 0 , 0 ,  2 )= 4.7384376191096816E-001   !  
       XHS(icarb,isulfur, 0 , 0 ,  3 )= 3.4763706058518457E-001   !  
       XHS(icarb,isulfur, 0 , 1 ,  0 )= 3.7965939239444002E-001   !  
       XHS(icarb,isulfur, 0 , 1 ,  1 )= 4.0897454706802949E-001   !  
       XHS(icarb,isulfur, 0 , 1 ,  2 )= 3.4796693492123965E-001   !  
       XHS(icarb,isulfur, 0 , 2 ,  0 )= 3.3888722745620037E-001   !  
       XHS(icarb,isulfur, 0 , 2 ,  1 )= 4.2934469243193241E-001   !  
       XHS(icarb,isulfur, 0 , 3 ,  0 )= 4.4053197067010474E-001   !  
       XHS(icarb,isulfur, 1 , 0 ,  0 )= 5.1997334061515188E-001   !  
       XHS(icarb,isulfur, 1 , 0 ,  1 )= 4.2880332834528923E-001   !  
       XHS(icarb,isulfur, 1 , 0 ,  2 )= 3.3939261090969586E-001   !  
       XHS(icarb,isulfur, 1 , 1 ,  0 )= 3.1424231405956898E-001   !  
       XHS(icarb,isulfur, 1 , 1 ,  1 )= 3.3798314862732465E-001   !  
       XHS(icarb,isulfur, 1 , 2 ,  0 )= 4.0990475940460358E-001   !  
       XHS(icarb,isulfur, 2 , 0 ,  0 )= 3.1751085086129271E-001   !  
       XHS(icarb,isulfur, 2 , 0 ,  1 )= 3.8480064082934806E-001   !  
       XHS(icarb,isulfur, 2 , 1 ,  0 )= 4.0000088357551822E-001   !  
       XHS(icarb,isulfur, 3 , 0 ,  0 )= 3.9660912277409577E-001   !  
       XHS(isulfur,icarb, 0 , 0 ,  0 )= -2.3148630194998890E-001   !  
       XHS(isulfur,icarb, 0 , 0 ,  1 )= -2.8438022092407877E-001   !  
       XHS(isulfur,icarb, 0 , 1 ,  0 )= -2.1384692439737088E-001   !  
       XHS(isulfur,icarb, 1 , 0 ,  0 )= -2.7653778493425218E-001   !  
       XHS(isulfur,ihyd, 0 , 0 ,  0 )= -1.8317967218679648E-003   !  
       XHS(isulfur,ihyd, 0 , 0 ,  1 )= -5.1966074866549921E-003   !  
       XHS(isulfur,ihyd, 0 , 1 ,  0 )= -6.2573236143447022E-003   !  
       XHS(isulfur,ihyd, 1 , 0 ,  0 )= -1.0703421761948984E-002   !  
       XHS(isulfur,isulfur, 0 , 0 ,  0 )= 6.9221312114058900E-004   !  
       XHS(isulfur,isulfur, 0 , 0 ,  1 )= 1.4050360544063662E-002   !  
       XHS(isulfur,isulfur, 0 , 1 ,  0 )= 1.1147118182644430E-002   !  
       XHS(isulfur,isulfur, 1 , 0 ,  0 )= 6.9983136135269583E-003   !  

!
! set S+ 3 coordinated structures to bond energy -0.038778eV (300K)
! Not used 
!       XHS(isulfur,icarb, 0 , 0 ,  2 )= -7.0361220999010543E-002   !  
!       XHS(isulfur,icarb, 0 , 1 ,  1 )= -6.9866760028864494E-002   !  
!       XHS(isulfur,icarb, 0 , 2 ,  0 )= -6.3648804942289416E-002   !  
!       XHS(isulfur,icarb, 1 , 0 ,  1 )= -2.0999664129477391E-001   !  
!       XHS(isulfur,icarb, 1 , 1 ,  0 )= -2.0873573983732285E-001   !  
!       XHS(isulfur,icarb, 2 , 0 ,  0 )= -2.0788754682742161E-001   !  
!       XHS(isulfur,ihyd, 0 , 0 ,  2 )= 3.4043137367705922E-002   !  
!       XHS(isulfur,ihyd, 0 , 1 ,  1 )= 3.4083174051413367E-002   !  
!       XHS(isulfur,ihyd, 0 , 2 ,  0 )= 3.4440123771939390E-002   !  
!       XHS(isulfur,ihyd, 1 , 0 ,  1 )= 3.3959443512503684E-002   !  
!       XHS(isulfur,ihyd, 1 , 1 ,  0 )= 3.4197253138897210E-002   !  
!       XHS(isulfur,ihyd, 2 , 0 ,  0 )= 1.7583895861427168E-002   !  
! Overcoordination penalty for S-S 
       XHS(isulfur,isulfur, 0 , 0 ,  2 )= 1.0887650521456638E-001   !  
       XHS(isulfur,isulfur, 0 , 1 ,  1 )= 1.0278991439624474E-001   !  
       XHS(isulfur,isulfur, 0 , 2 ,  0 )= 1.0536816677057703E-001   !  
       XHS(isulfur,isulfur, 1 , 0 ,  1 )= 1.0625804625879121E-001   !  
       XHS(isulfur,isulfur, 1 , 1 ,  0 )= 1.0507661825674175E-001   !  
       XHS(isulfur,isulfur, 2 , 0 ,  0 )= 4.9400255424925099E-002   !  
!
C Calculate Spline Derivatives
      DO KI=1,RTYPES
        DO KJ=1,RTYPES
          DO NTC=1,8
            DO NTH=0,9
              DO NTS=0,9
                XHS1(KI,KJ,NTC,NTH,NTS)=(XHS(KI,KJ,NTC+1,NTH,NTS)
     &           -XHS(KI,KJ,NTC-1,NTH,NTS))/2.0D0  
              ENDDO   
            ENDDO   
          ENDDO    
        ENDDO   
      ENDDO    
      DO KI=1,RTYPES
        DO KJ=1,RTYPES
          DO NTC=0,9
            DO NTH=1,8
              DO NTS=0,9
                XHS2(KI,KJ,NTC,NTH,NTS)=(XHS(KI,KJ,NTC,NTH+1,NTS)
     &           -XHS(KI,KJ,NTC,NTH-1,NTS))/2.0D0  
              ENDDO   
            ENDDO   
          ENDDO    
        ENDDO   
      ENDDO    
      DO KI=1,RTYPES
        DO KJ=1,RTYPES
          DO NTC=0,9
            DO NTH=0,9
              DO NTS=1,8
                XHS3(KI,KJ,NTC,NTH,NTS)=(XHS(KI,KJ,NTC,NTH,NTS+1)
     &           -XHS(KI,KJ,NTC,NTH,NTS-1))/2.0D0  
              ENDDO   
            ENDDO   
          ENDDO    
        ENDDO   
      ENDDO    
      DO KI=1,RTYPES
        DO KJ=1,RTYPES
          DO NTC=1,8
            DO NTH=1,8
              DO NTS=1,8
                dx1 = XHS1(KI,KJ,NTC,NTH+1,NTS)
                dx2 = XHS1(KI,KJ,NTC,NTH-1,NTS)
                XHS12(KI,KJ,NTC,NTH,NTS) = (dx1 -  dx2)/2.0d0
                dx1 = XHS1(KI,KJ,NTC,NTH,NTS+1)
                dx2 = XHS1(KI,KJ,NTC,NTH,NTS-1)
                XHS13(KI,KJ,NTC,NTH,NTS) = (dx1 -  dx2)/2.0d0
                dy1 = XHS2(KI,KJ,NTC,NTH,NTS+1)
                dy2 = XHS2(KI,KJ,NTC,NTH,NTS-1)
                XHS23(KI,KJ,NTC,NTH,NTS) = (dy1 -  dy2)/2.0d0
                dz1 = XHS3(KI,KJ,NTC+1,NTH+1,NTS)
                dz2 = XHS3(KI,KJ,NTC-1,NTH+1,NTS)
                dy1 = (dz1 - dz2)/2.0d0
                dz3 = XHS3(KI,KJ,NTC+1,NTH-1,NTS)
                dz4 = XHS3(KI,KJ,NTC-1,NTH-1,NTS)
                dy2 = (dz3 - dz4)/2.0d0                
                XHS123(KI,KJ,NTC,NTH,NTS) = (dy1 -  dy2)/2.0d0
              ENDDO
            ENDDO   
          ENDDO    
        ENDDO   
      ENDDO    
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
C  Calculate spline coefficients
C  -travisk
C
      DO KI=1,RTYPES
        DO KJ=1,RTYPES
          DO NTC=0,4
           DO NTH=0,4
            DO NTS=0,4
                  xnlow = FLOAT(NTC)
                  xnhg = FLOAT(NTC + 1)
                  ynlow = FLOAT(NTH)
                  ynhg = FLOAT(NTH + 1)
                  znlow = FLOAT(NTS)
                  znhg = FLOAT(NTS + 1)

                  xmin = (NTC)
                  xmax = (NTC + 1)
                  ymin = (NTH)
                  ymax = (NTH + 1)
                  zmin = (NTS)
                  zmax = (NTS + 1)


                  f(0,0,0) = XHS(KI,KJ,xmin,ymin,zmin)
                  f(1,0,0) = XHS(KI,KJ,xmax,ymin,zmin)
                  f(0,1,0) = XHS(KI,KJ,xmin,ymax,zmin)
                  f(0,0,1) = XHS(KI,KJ,xmin,ymin,zmax)
                  f(1,1,0) = XHS(KI,KJ,xmax,ymax,zmin)
                  f(0,1,1) = XHS(KI,KJ,xmin,ymax,zmax)
                  f(1,1,1) = XHS(KI,KJ,xmax,ymax,zmax)

                  fx(0,0,0) = XHS1(KI,KJ,xmin,ymin,zmin)
                  fx(1,0,0) = XHS1(KI,KJ,xmax,ymin,zmin)
                  fx(0,1,0) = XHS1(KI,KJ,xmin,ymax,zmin)
                  fx(0,0,1) = XHS1(KI,KJ,xmin,ymin,zmax)
                  fx(1,1,0) = XHS1(KI,KJ,xmax,ymax,zmin)
                  fx(0,1,1) = XHS1(KI,KJ,xmin,ymax,zmax)
                  fx(1,1,1) = XHS1(KI,KJ,xmax,ymax,zmax)

                  fy(0,0,0) = XHS2(KI,KJ,xmin,ymin,zmin)
                  fy(1,0,0) = XHS2(KI,KJ,xmax,ymin,zmin)
                  fy(0,1,0) = XHS2(KI,KJ,xmin,ymax,zmin)
                  fy(0,0,1) = XHS2(KI,KJ,xmin,ymin,zmax)
                  fy(1,1,0) = XHS2(KI,KJ,xmax,ymax,zmin)
                  fy(0,1,1) = XHS2(KI,KJ,xmin,ymax,zmax)
                  fy(1,1,1) = XHS2(KI,KJ,xmax,ymax,zmax)

                  fz(0,0,0) = XHS3(KI,KJ,xmin,ymin,zmin)
                  fz(1,0,0) = XHS3(KI,KJ,xmax,ymin,zmin)
                  fz(0,1,0) = XHS3(KI,KJ,xmin,ymax,zmin)
                  fz(0,0,1) = XHS3(KI,KJ,xmin,ymin,zmax)
                  fz(1,1,0) = XHS3(KI,KJ,xmax,ymax,zmin)
                  fz(0,1,1) = XHS3(KI,KJ,xmin,ymax,zmax)
                  fz(1,1,1) = XHS3(KI,KJ,xmax,ymax,zmax)

                  fxy(0,0,0) = XHS12(KI,KJ,xmin,ymin,zmin)
                  fxy(1,0,0) = XHS12(KI,KJ,xmax,ymin,zmin)
                  fxy(0,1,0) = XHS12(KI,KJ,xmin,ymax,zmin)
                  fxy(0,0,1) = XHS12(KI,KJ,xmin,ymin,zmax)
                  fxy(1,1,0) = XHS12(KI,KJ,xmax,ymax,zmin)
                  fxy(0,1,1) = XHS12(KI,KJ,xmin,ymax,zmax)
                  fxy(1,1,1) = XHS12(KI,KJ,xmax,ymax,zmax)

                  fxz(0,0,0) = XHS13(KI,KJ,xmin,ymin,zmin)
                  fxz(1,0,0) = XHS13(KI,KJ,xmax,ymin,zmin)
                  fxz(0,1,0) = XHS13(KI,KJ,xmin,ymax,zmin)
                  fxz(0,0,1) = XHS13(KI,KJ,xmin,ymin,zmax)
                  fxz(1,1,0) = XHS13(KI,KJ,xmax,ymax,zmin)
                  fxz(0,1,1) = XHS13(KI,KJ,xmin,ymax,zmax)
                  fxz(1,1,1) = XHS13(KI,KJ,xmax,ymax,zmax)

                  fyz(0,0,0) = XHS23(KI,KJ,xmin,ymin,zmin)
                  fyz(1,0,0) = XHS23(KI,KJ,xmax,ymin,zmin)
                  fyz(0,1,0) = XHS23(KI,KJ,xmin,ymax,zmin)
                  fyz(0,0,1) = XHS23(KI,KJ,xmin,ymin,zmax)
                  fyz(1,1,0) = XHS23(KI,KJ,xmax,ymax,zmin)
                  fyz(0,1,1) = XHS23(KI,KJ,xmin,ymax,zmax)
                  fyz(1,1,1) = XHS23(KI,KJ,xmax,ymax,zmax)

                  fxyz(0,0,0) = XHS123(KI,KJ,xmin,ymin,zmin)
                  fxyz(1,0,0) = XHS123(KI,KJ,xmax,ymin,zmin)
                  fxyz(0,1,0) = XHS123(KI,KJ,xmin,ymax,zmin)
                  fxyz(0,0,1) = XHS123(KI,KJ,xmin,ymin,zmax)
                  fxyz(1,1,0) = XHS123(KI,KJ,xmax,ymax,zmin)
                  fxyz(0,1,1) = XHS123(KI,KJ,xmin,ymax,zmax)
                  fxyz(1,1,1) = XHS123(KI,KJ,xmax,ymax,zmax)
!
                  call tcucofair(xnlow,xnhg,ynlow,ynhg,znlow,znhg,
     &                 f, fx, fy, fz, fxy, fxz, fyz, fxyz,TRICOEF )
                  ii = 0
                  DO k = 0,3
                    DO j = 0,3
                      DO i= 0,3
                        ii = ii + 1
                        CLMSX(ii,NTC,NTH,NTS,KI,KJ) =TRICOEF(i,j,k)
                      ENDDO
                    ENDDO
                  ENDDO
! begin debug
!           WRITE(53,*) 'Coefficients for',KI,KJ,NTC,NTH,NTS
!           WRITE(53,*) CLMSX(KI,KJ,NTC,NTH,NTS,1:64)
! end debug
                ENDDO
              ENDDO
            ENDDO
         ENDDO
      ENDDO
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! FLOURINE-CARBON
!
       AD(iflor,icarb)=219.7799D0
       AXL(iflor,icarb)=2.1763D0
       DD(iflor,icarb)=909.2022D0
       DXL(iflor,icarb)= 3.7128D0
       ED(iflor,icarb)=0.0D0
	
      BD(iflor,icarb)=0.0D0
      BXL(iflor,icarb)=1.0D0
      CD(iflor,icarb)=0.0D0
      CXL(iflor,icarb)=1.0D0
C
      AD(icarb,iflor)=AD(iflor,icarb)
      AXL(icarb,iflor)=AXL(iflor,icarb)
      BD(icarb,iflor)=BD(iflor,icarb)
      BXL(icarb,iflor)=BXL(iflor,icarb)
      CD(icarb,iflor)=CD(iflor,icarb)
      CXL(icarb,iflor)=CXL(iflor,icarb)
      DD(icarb,iflor)=DD(iflor,icarb)
      DXL(icarb,iflor)=DXL(iflor,icarb)
      ED(icarb,iflor)=ED(iflor,icarb)
!
      RB1(iflor,icarb)=1.7d0
      RB2(iflor,icarb)=2.0d0
!
      RMAX(iflor,icarb)=RB2(iflor,icarb)
      PID(iflor,icarb)=PI/(RB2(iflor,icarb)-RB1(iflor,icarb))      
!
      RB1(icarb,iflor)=RB1(iflor,icarb)
      RB2(icarb,iflor)=RB2(iflor,icarb)
!
      RMAX(icarb,iflor)=RB2(icarb,iflor)
      PID(icarb,iflor)=PI/(RB2(icarb,iflor)-RB1(icarb,iflor))
!
! 	Fluorine
!
      AD(iflor,iflor)=146.8149D0
      AXL(iflor,iflor)=2.8568D0
      BD(iflor,iflor)=0.0D0
      BXL(iflor,iflor)=1.0D0
      CD(iflor,iflor)=0.0D0
      CXL(iflor,iflor)=1.0D0
      DD(iflor,iflor)=16451.97D0
      DXL(iflor,iflor)=6.8149D0
      ED(iflor,iflor)=0.0d0
!
      RB1(iflor,iflor)=1.7d0
      RB2(iflor,iflor)=2.0d0
      RMAX(iflor,iflor)=RB2(iflor,iflor)
      PID(iflor,iflor)=PI/(RB2(iflor,iflor)-RB1(iflor,iflor))
!
! 	Hydrogen-Fluorine
!
      AD(iflor,ihyd)=571.1737214533D0
      AXL(iflor,ihyd)=3.0919704916D0
      BD(iflor,ihyd)=0.0D0
      BXL(iflor,ihyd)=1.0D0
      CD(iflor,ihyd)=0.0D0
      CXL(iflor,ihyd)=1.0D0
      DD(iflor,ihyd)=887.051318618D0
      DXL(iflor,ihyd)=3.7789072764D0
      ED(iflor,ihyd)=0.0d0
!
      RB1(iflor,ihyd)=1.3d0
      RB2(iflor,ihyd)=1.8d0
      RMAX(iflor,ihyd)=RB2(iflor,ihyd)
      PID(iflor,ihyd)=PI/(RB2(iflor,ihyd)-RB1(iflor,ihyd))
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!      
      AD(ihyd,iflor)=AD(iflor,ihyd)
      AXL(ihyd,iflor)=AXL(iflor,ihyd)
      BD(ihyd,iflor)=BD(iflor,ihyd)
      BXL(ihyd,iflor)=BXL(iflor,ihyd)
      CD(ihyd,iflor)=CD(iflor,ihyd)
      CXL(ihyd,iflor)=CXL(iflor,ihyd)
      DD(ihyd,iflor)=DD(iflor,ihyd)
      DXL(ihyd,iflor)=DXL(iflor,ihyd)
      ED(ihyd,iflor)=ED(iflor,ihyd)

      RB1(ihyd,iflor)=RB1(iflor,ihyd)
      RB2(ihyd,iflor)=RB2(iflor,ihyd)
  
      RMAX(ihyd,iflor)=RB2(ihyd,iflor)
      PID(ihyd,iflor)=PI/(RB2(ihyd,iflor)-RB1(ihyd,iflor))
!
      AD(icarb,ihyd)=AD(ihyd,icarb)
      AXL(icarb,ihyd)=AXL(ihyd,icarb)
      BD(icarb,ihyd)=BD(ihyd,icarb)
      BXL(icarb,ihyd)=BXL(ihyd,icarb)
      CD(icarb,ihyd)=CD(ihyd,icarb)
      CXL(icarb,ihyd)=CXL(ihyd,icarb)
      DD(icarb,ihyd)=DD(ihyd,icarb)
      DXL(icarb,ihyd)=DXL(ihyd,icarb)
      ED(icarb,ihyd)=ED(ihyd,icarb)
!
      AD(icarb,ioxygen)=AD(ioxygen,icarb)
      AXL(icarb,ioxygen)=AXL(ioxygen,icarb)
      BD(icarb,ioxygen)=BD(ioxygen,icarb)
      BXL(icarb,ioxygen)=BXL(ioxygen,icarb)
      CD(icarb,ioxygen)=CD(ioxygen,icarb)
      CXL(icarb,ioxygen)=CXL(ioxygen,icarb)
      DD(icarb,ioxygen)=DD(ioxygen,icarb)
      DXL(icarb,ioxygen)=DXL(ioxygen,icarb)
      ED(icarb,ioxygen)=ED(ioxygen,icarb)
!
      AD(ihyd,ioxygen)=AD(ioxygen,ihyd)
      AXL(ihyd,ioxygen)=AXL(ioxygen,ihyd)
      BD(ihyd,ioxygen)=BD(ioxygen,ihyd)
      BXL(ihyd,ioxygen)=BXL(ioxygen,ihyd)
      CD(ihyd,ioxygen)=CD(ioxygen,ihyd)
      CXL(ihyd,ioxygen)=CXL(ioxygen,ihyd)

      DD(ihyd,ioxygen)=DD(ioxygen,ihyd)
      DXL(ihyd,ioxygen)=DXL(ioxygen,ihyd)
      ED(ihyd,ioxygen)=ED(ioxygen,ihyd)
!
!dong, these are oxygen related stuff from Susan
      ROH=0.96d0
!
!    Pair potential values for Sulfur
!
      AD(icarb,isulfur)=AD(isulfur,icarb)
      AXL(icarb,isulfur)=AXL(isulfur,icarb)
      BD(icarb,isulfur)=BD(isulfur,icarb)
      BXL(icarb,isulfur)=BXL(isulfur,icarb)
      CD(icarb,isulfur)=CD(isulfur,icarb)
      CXL(icarb,isulfur)=CXL(isulfur,icarb)
      DD(icarb,isulfur)=DD(isulfur,icarb)
      DXL(icarb,isulfur)=DXL(isulfur,icarb)
      ED(icarb,isulfur)=ED(isulfur,icarb)
!
      AD(ihyd,isulfur)=AD(isulfur,ihyd)
      AXL(ihyd,isulfur)=AXL(isulfur,ihyd)
      BD(ihyd,isulfur)=BD(isulfur,ihyd)
      BXL(ihyd,isulfur)=BXL(isulfur,ihyd)
      CD(ihyd,isulfur)=CD(isulfur,ihyd)
      CXL(ihyd,isulfur)=CXL(isulfur,ihyd)
!
      DD(ihyd,isulfur)=DD(isulfur,ihyd)
      DXL(ihyd,isulfur)=DXL(isulfur,ihyd)
      ED(ihyd,isulfur)=ED(isulfur,ihyd)
!
! PID = frac{/pi}{(r_{cut}^{outer} - r_{cut}^{inner})}
!
      RMAX(icarb,icarb)=RB2(icarb,icarb)
      PID(icarb,icarb)=PI/(RB2(icarb,icarb)-RB1(icarb,icarb))

      RMAX(icarb,ihyd)=RB2(icarb,ihyd)
      PID(icarb,ihyd)=PI/(RB2(icarb,ihyd)-RB1(icarb,ihyd))

      RMAX(ihyd,icarb)=RB2(ihyd,icarb)
      PID(ihyd,icarb)=PI/(RB2(ihyd,icarb)-RB1(ihyd,icarb))
      PIDT=PI/0.30D0

      RMAX(ihyd,ihyd)=RB2(ihyd,ihyd)
      PID(ihyd,ihyd)=PI/(RB2(ihyd,ihyd)-RB1(ihyd,ihyd))
      RMAX(ioxygen,ioxygen)=RB2(ioxygen,ioxygen)
      PID(ioxygen,ioxygen)=PI/
     & (RB2(ioxygen,ioxygen)-RB1(ioxygen,ioxygen))

      RMAX(icarb,ioxygen)=RB2(icarb,ioxygen)
      PID(icarb,ioxygen)=PI/(RB2(icarb,ioxygen)-RB1(icarb,ioxygen))
      RB1(ioxygen,icarb)=RB1(icarb,ioxygen)
      RB2(ioxygen,icarb)=RB2(icarb,ioxygen)
      RMAX(ioxygen,icarb)=RB2(ioxygen,icarb)
      PID(ioxygen,icarb)=PI/(RB2(ioxygen,icarb)-RB1(ioxygen,icarb))


      RMAX(ihyd,ioxygen)=RB2(ihyd,ioxygen)
      PID(ihyd,ioxygen)=PI/(RB2(ihyd,ioxygen)-RB1(ihyd,ioxygen))
      RB1(ioxygen,ihyd)=RB1(ihyd,ioxygen)
      RB2(ioxygen,ihyd)=RB2(ihyd,ioxygen)
      RMAX(ioxygen,ihyd)=RB2(ioxygen,ihyd)
      PID(ioxygen,ihyd)=PI/(RB2(ioxygen,ihyd)-RB1(ioxygen,ihyd))

      RMAX(isulfur,isulfur)=RB2(isulfur,isulfur)
      PID(isulfur,isulfur)=PI/
     & (RB2(isulfur,isulfur)-RB1(isulfur,isulfur))

      RMAX(icarb,isulfur)=RB2(icarb,isulfur)
      PID(icarb,isulfur)=PI/(RB2(icarb,isulfur)-RB1(icarb,isulfur))
      RB1(isulfur,icarb)=RB1(icarb,isulfur)
      RB2(isulfur,icarb)=RB2(icarb,isulfur)
      RMAX(isulfur,icarb)=RB2(isulfur,icarb)
      PID(isulfur,icarb)=PI/(RB2(isulfur,icarb)-RB1(isulfur,icarb))

      RMAX(ihyd,isulfur)=RB2(ihyd,isulfur)
      PID(ihyd,isulfur)=PI/(RB2(ihyd,isulfur)-RB1(ihyd,isulfur))
      RB1(isulfur,ihyd)=RB1(ihyd,isulfur)
      RB2(isulfur,ihyd)=RB2(ihyd,isulfur)
      RMAX(isulfur,ihyd)=RB2(isulfur,ihyd)
      PID(isulfur,ihyd)=PI/(RB2(isulfur,ihyd)-RB1(isulfur,ihyd))


      WRITE(13,*) 'Checking cut-off values:'
      DO I = 1, RTYPES
        DO J = 1,RTYPES
         IF (RB1(I,J).GT.0.0d0) THEN
           RB1(J,I) = RB1(I,J)
         ELSEIF  (RB1(I,J).GT.0.0d0) THEN
           RB1(J,I) = RB1(I,J)
         ELSE
           WRITE(13,*) 'RB1 not set for',I,J
         ENDIF
         IF (RB2(I,J).GT.0.0d0) THEN
           RB2(J,I) = RB2(I,J)
         ELSEIF  (RB2(I,J).GT.0.0d0) THEN
           RB2(J,I) = RB2(I,J)
         ELSE
           WRITE(13,*) 'RB2 not set for',I,J
         ENDIF
         RMAX(I,J) = RB2(I,J)
         PID(I,J) = PI/(RB2(I,J) - RB1(I,J))
        ENDDO
      ENDDO

      DO 12 I=1,RTYPES
           DO 11 J=1,RTYPES
                DO 10 K=1,RTYPES
                     XDB(I,J,K)=0.0d0
                     REG(I,J,K)=1.0d0
10              CONTINUE
11         CONTINUE
12    CONTINUE 
!
      XXDB=4.0D0
      RHH=0.7415886997d0
      RCH=1.09d0
      RCF=1.2718d0
      RFF=1.4119d0
      RHF=0.9378d0
      XDB(2,2,2)=4.0D0
C
      XDB(2,1,2)=4.0D0
      XDB(2,2,1)=4.0D0
C
      XDB(2,1,1)=4.0D0
C
      REG(2,1,2)=EXP(XDB(2,1,2)*(RHH-RCH))
      REG(2,2,1)=EXP(XDB(2,2,1)*(RCH-RHH))
C
      XDB(iflor,iflor,iflor)=4.0D0
!
      XDB(iflor,1,iflor)=4.0D0
      XDB(iflor,iflor,1)=4.0D0
!
      XDB(iflor,1,1)=4.0D0
      XDB(iflor,2,2)=4.0D0
      XDB(iflor,2,1)=4.0D0
      XDB(iflor,1,2)=4.0D0
      XDB(iflor,2,iflor)=4.0D0
      XDB(iflor,iflor,2)=4.0D0
!
      XDB(1,iflor,1)=0.0D0
      XDB(1,iflor,iflor)=0.0d0
C
! Quick fix for F3 not sure of implications
!      REG(iflor,iflor,iflor)=800.0d0
!     end fix
      REG(iflor,1,iflor)=EXP(XDB(iflor,1,iflor)*(RFF-RCF))
      REG(iflor,iflor,1)=EXP(XDB(iflor,iflor,1)*(RCF-RFF))

      REG(iflor,1,2)=EXP(XDB(iflor,1,2)*(RHF-RCF))
      REG(iflor,2,1)=EXP(XDB(iflor,2,1)*(RCF-RHF))

      DO 61 I=1,RTYPES
           DO 60 J=1,RTYPES
                RLIST(I,J)=(RMAX(I,J)+RLL)**2
                RMAX(I,J)=RMAX(I,J)**2
   60      CONTINUE
   61 CONTINUE

      RETURN
110   FORMAT(4I5)
120   FORMAT(4E20.6)
  
      END SUBROUTINE PARAM
