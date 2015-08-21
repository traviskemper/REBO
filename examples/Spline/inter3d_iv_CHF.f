C
C FILE INTER3D_iv_CHF.F
c
c Create input file coef3d
c
      IMPLICIT REAL*8(A-H,O-Z)
      integer ix,jj,kk
C
      dimension IN(64,3),II(8,3),A(64,64),B(64),C(64)
      DIMENSION YYC(4,4,4),YPC(10,10,10,3),
     & YYH(4,4,4),YPH(10,10,10,3),YYF(4,4,4),YPF(10,10,10,3)
C
      DATA II/0,1,0,1,0,1,0,1,0,0,1,1,0,0,1,1,0,0,0,0,1,1,1,1/
C
      DATA YYC/64*0.0D0/
      DATA YPC/3000*0.0D0/
      DATA YYH/64*0.0D0/
      DATA YPH/3000*0.0D0/
      DATA YYF/64*0.0D0/
      DATA YPF/3000*0.0D0/
      data IN/192*0/
C
      REAL RWKSP(8406)
      CALL IWKIN(8406)
C
C     Y(NH,NC,NF)
C
C

      write(1,*) 4

C From CH
c done
      yyh(1,1,1)=0.0000000000d0
C From CH2
c done
      yyh(2,1,1)=0.2093367328250380
C From CH3
c done
      yyh(3,1,1)=-6.4449615432525347E-02
C From CH4
c done
      yyh(4,1,1)=-0.3039275463461620
C From C2H2
c done
      yyh(1,2,1)=0.01D0
      yyc(2,1,1)=0.0D0
C From (CH3)HC=CH(CH3)
c done
      yyc(2,2,1)=3.0266974734813094E-03
      yyh(1,3,1)=-0.1220421462782555
C From C2H4
c done
       yyc(3,1,1)=7.8607002547458685E-03
       yyh(2,2,1)=-0.1251234006287090
C From C2H6
c done
      yyc(4,1,1)=1.6125364564267377E-02
      yyh(3,2,1)= -0.2989052457826418
C From i-C4H10
c done
      yyc(2,3,1)=3.1795308307312573E-03
      yyh(1,4,1)=-0.3075847050665519
C From c-C6H12
c done
      yyc(3,2,1)=6.3262482411195846E-03
      yyh(2,3,1)=-0.3005291724067579

c CF
      YYF(1,1,1)=0.0D0

c CF2
      YYF(1,1,2)=0.19024451572976
C CHF
      YYF(2,1,1)=0.19024451572976
      YYH(1,1,2)=0.2093367328250380
C CF3
      YYF(1,1,3)=0.3485710783254
C CHF2
      yyf(2,1,2)=0.338998160620986
      YYH(1,1,3)=-0.05232782973885675
C CH2F
      YYF(3,1,1)=0.3299679089558
      YYH(2,1,2)=-0.05911750735142984
C CF4
      YYF(1,1,4)=-0.3039275463461620
C CHF3
      YYF(2,1,3)=-0.312969839786
      YYH(1,1,4)=-0.282462006562285
C CH2F2
      YYF(3,1,2)=-0.314008475695567
      YYH(2,1,3)=-0.295893104887462
C CH3F
      YYF(4,1,1)=-0.305552728291
      YYH(3,1,2)=-0.303367434289228
C C2F2
      YYF(1,2,1)=0.35D0
      YYC(1,1,2)=0.0D0
C (CH3)FC=CF(CH3)
      YYF(1,3,1)=0.40310716301456
      YYC(1,2,2)=0.003405089871339
C C2F4
      YYF(1,2,2)=0.1329152919986
      YYC(1,1,3)=0.00880500239758
c C2H2F2
      YYF(2,2,1)=0.133000000495
      YYH(1,2,2)=-0.125673698097359
      YYC(2,1,2)=0.0034677627824
C C2F6
      YYF(1,2,3)=-0.29218488638
      YYC(1,1,4)=0.0019058417195
C C2H2F4
      YYF(2,2,2)=-0.290995698488
      YYH(1,2,3)=-0.283830004269
      YYC(2,1,3)=-0.0085455706363
C C2H4F2
      YYF(3,2,1)=-0.274974489660
      YYH(2,2,2)=-0.297831534909713
      YYC(3,1,2)=-0.00577350464522
C i-C4H9F
      yyf(1,4,1)=0.06844317489
      yyc(1,3,2)=0.00153020602237
c c-C6F12
      yyf(1,3,2)=-0.20259034048
      yyc(1,2,3)=-0.000299517913389
c c-C6H6F6
      yyf(2,3,1)=-0.1829045094487
      yyh(1,3,2)=-0.312299492750992
      yyc(2,2,2)=-0.000677461804849
 
C
C Use centered derivatives for slope
C
      YPH(3,1,1,1)=(YYH(4,1,1)-YYH(2,1,1))/2.0D0
      YPH(2,2,1,1)=(YYH(3,2,1)-YYH(1,2,1))/2.0D0
C
      YPH(2,2,1,2)=(YYH(2,3,1)-YYH(2,1,1))/2.0D0
      YPH(1,3,1,2)=(YYH(1,4,1)-YYH(1,2,1))/2.0D0
cccc
      YPF(1,1,3,3)=(YYF(1,1,4)-YYF(1,1,2))/2.0D0
      YPF(1,2,2,3)=(YYF(1,2,3)-YYF(1,2,1))/2.0D0
C
      YPF(1,2,2,2)=(YYF(1,3,2)-YYF(1,1,2))/2.0D0
      YPF(1,3,1,2)=(YYF(1,4,1)-YYF(1,2,1))/2.0D0

      YPH(1,1,3,3)=(YYH(1,1,4)-YYH(1,1,2))/2.0D0
      YPF(3,1,1,1)=(YYF(4,1,1)-YYF(2,1,1))/2.0D0

      YPH(1,2,2,3)=(YYH(1,2,3)-YYH(1,2,1))/2.0D0
      YPF(2,2,1,1)=(YYF(3,2,1)-YYF(1,2,1))/2.0D0

      YPH(1,2,2,2)=(YYH(1,3,2)-YYH(1,1,2))/2.0D0
      YPF(2,2,1,2)=(YYF(2,3,1)-YYF(2,1,1))/2.0D0

c 

      do 2 ix=1,4
        do 3 jj=1,4
          do 444 kk=1,4
            if ((ix+jj+kk).gt.6) goto 444
            write(1,*) 2,ix,jj,kk,yyh(ix,jj,kk)
            write(1,*) 1,ix,jj,kk,yyc(ix,jj,kk)
            write(1,*) 3,ix,jj,kk,yyf(ix,jj,kk)
444       continue
3       continue
2     continue

      write(1,*) -1,-1,-1,-1,0.0
C
c***
C
C
      IC=0
C
      DO 6 I=1,4
           DO 5 J=1,4
                DO 4 K=1,4
                     IC=IC+1
                     IN(IC,1)=I-1
                     IN(IC,2)=J-1
                     IN(IC,3)=K-1
4               CONTINUE
5          CONTINUE
6     CONTINUE
C
c     Carbon-Carbon

      DO 30 L=1,6
C
           DO 20 M=1,6
C
                DO 10 N=1,6
C
                     IC=-7
C
                     DO 9 I=1,8
C
                          IC=IC+8
C
                          LI=L+II(I,1)
                          MI=M+II(I,2)
                          NI=N+II(I,3)
                          B(IC)=YYC(LI,MI,NI)
                          B(IC+1)=YPC(LI,MI,NI,1)
                          B(IC+2)=YPC(LI,MI,NI,2)
                          B(IC+3)=YPC(LI,MI,NI,3)
                          B(IC+4)=0.0D0
                          B(IC+5)=0.0D0
                          B(IC+6)=0.0D0
                          B(IC+7)=0.0D0
                          T=LI
                          U=MI
                          V=NI
                          DO 8 J=1,64
                               Y=(T**IN(J,1))*
     &                           (U**IN(J,2))*(V**IN(J,3))
                               A(IC,J)=Y
                               A(IC+1,J)=IN(J,1)*Y/T
                               A(IC+2,J)=IN(J,2)*Y/U
                               A(IC+3,J)=IN(J,3)*Y/V
                               A(IC+4,J)=Y*IN(J,1)/T*IN(J,2)/U
                               A(IC+5,J)=Y*IN(J,1)/T*IN(J,3)/V
                               A(IC+6,J)=Y*IN(J,2)/U*IN(J,3)/V
                               A(IC+7,J)=Y*IN(J,1)/T*IN(J,2)/U*IN(J,3)/V
8                         CONTINUE
9                    CONTINUE
c
C
C    DSLARG(N,A,LDA,B,IPATH,X) solves a set of linear equations
C
C     N: Number of equations (input)
C
C     A: NxN matrix containing the coefficients of the linear system (input)
C
C   LDA: Leading dimension of A (input)
C
C     B: Vector of length N containing the RHS of the linear system (input)
C
C IPATH: Path indicator. (input)
C                        IPATH = 1 A*X = B solved
C                        IPATH = 2 trans(A)*X = B solved
C
C     X: Vector of length N containing solution to the linear system (output).
C
                     CALL DLSARG(64,A,64,B,1,C)
C
                      WRITE(1,110) 1,L,M,N
                      WRITE(1,120) (C(I),I=1,64)

C                
10              CONTINUE
20         CONTINUE
30    CONTINUE

c     Carbon-Hydrogen

      DO 31 L=1,6
C
           DO 21 M=1,6
C
                DO 11 N=1,6
C
                     IC=-7
C
                     DO 19 I=1,8
C
                          IC=IC+8
C
                          LI=L+II(I,1)
                          MI=M+II(I,2)
                          NI=N+II(I,3)
                          B(IC)=YYH(LI,MI,NI)
                          B(IC+1)=YPH(LI,MI,NI,1)
                          B(IC+2)=YPH(LI,MI,NI,2)
                          B(IC+3)=YPH(LI,MI,NI,3)
                          B(IC+4)=0.0D0
                          B(IC+5)=0.0D0
                          B(IC+6)=0.0D0
                          B(IC+7)=0.0D0
                          T=LI
                          U=MI
                          V=NI
                          DO 18 J=1,64
                               Y=(T**IN(J,1))*
     &                           (U**IN(J,2))*(V**IN(J,3))
                               A(IC,J)=Y
                               A(IC+1,J)=IN(J,1)*Y/T
                               A(IC+2,J)=IN(J,2)*Y/U
                               A(IC+3,J)=IN(J,3)*Y/V
                               A(IC+4,J)=Y*IN(J,1)/T*IN(J,2)/U
                               A(IC+5,J)=Y*IN(J,1)/T*IN(J,3)/V
                               A(IC+6,J)=Y*IN(J,2)/U*IN(J,3)/V
                               A(IC+7,J)=Y*IN(J,1)/T*IN(J,2)/U*IN(J,3)/V
18                         CONTINUE
19                    CONTINUE
c
C
C    DSLARG(N,A,LDA,B,IPATH,X) solves a set of linear equations
C
C     N: Number of equations (input)
C
C     A: NxN matrix containing the coefficients of the linear system (input)
C
C   LDA: Leading dimension of A (input)
C
C     B: Vector of length N containing the RHS of the linear system (input)
C
C IPATH: Path indicator. (input)
C                        IPATH = 1 A*X = B solved
C                        IPATH = 2 trans(A)*X = B solved
C
C     X: Vector of length N containing solution to the linear system (output).
C
                     CALL DLSARG(64,A,64,B,1,C)
C
                      WRITE(1,110) 2,L,M,N
                      WRITE(1,120) (C(I),I=1,64)

C                
11              CONTINUE
21         CONTINUE
31    CONTINUE

c     Carbon-fluorine

      DO 32 L=1,6
C
           DO 22 M=1,6
C
                DO 12 N=1,6
C
                     IC=-7
C
                     DO 29 I=1,8
C
                          IC=IC+8
C
                          LI=L+II(I,1)
                          MI=M+II(I,2)
                          NI=N+II(I,3)
                          B(IC)=YYF(LI,MI,NI)
                          B(IC+1)=YPF(LI,MI,NI,1)
                          B(IC+2)=YPF(LI,MI,NI,2)
                          B(IC+3)=YPF(LI,MI,NI,3)
                          B(IC+4)=0.0D0
                          B(IC+5)=0.0D0
                          B(IC+6)=0.0D0
                          B(IC+7)=0.0D0
                          T=LI
                          U=MI
                          V=NI
                          DO 28 J=1,64
                               Y=(T**IN(J,1))*
     &                           (U**IN(J,2))*(V**IN(J,3))
                               A(IC,J)=Y
                               A(IC+1,J)=IN(J,1)*Y/T
                               A(IC+2,J)=IN(J,2)*Y/U
                               A(IC+3,J)=IN(J,3)*Y/V
                               A(IC+4,J)=Y*IN(J,1)/T*IN(J,2)/U
                               A(IC+5,J)=Y*IN(J,1)/T*IN(J,3)/V
                               A(IC+6,J)=Y*IN(J,2)/U*IN(J,3)/V
                               A(IC+7,J)=Y*IN(J,1)/T*IN(J,2)/U*IN(J,3)/V
28                         CONTINUE
29                    CONTINUE
c
C
C    DSLARG(N,A,LDA,B,IPATH,X) solves a set of linear equations
C
C     N: Number of equations (input)
C
C     A: NxN matrix containing the coefficients of the linear system (input)
C
C   LDA: Leading dimension of A (input)
C
C     B: Vector of length N containing the RHS of the linear system (input)
C
C IPATH: Path indicator. (input)
C                        IPATH = 1 A*X = B solved
C                        IPATH = 2 trans(A)*X = B solved
C
C     X: Vector of length N containing solution to the linear system (output).
C
                     CALL DLSARG(64,A,64,B,1,C)
C
                      WRITE(1,110) 3,L,M,N
                      WRITE(1,120) (C(I),I=1,64)

C                
12              CONTINUE
22         CONTINUE
32    CONTINUE

100   FORMAT(16F5.2)
110   FORMAT(4I5)
120   FORMAT(4E20.11)
      STOP
      END
