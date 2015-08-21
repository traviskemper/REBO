      SUBROUTINE PIBONDCH
!
      USE MDSTP
      USE STRUCTR
      USE PARAMS
      USE POTS
      USE SPECIF
!
      IMPLICIT none
!
      INTEGER :: I,J,K,L,M,N,JBEGIN,JEND,JN,MM,KJ,KIKJ,NK,IG,LBEGIN,LEND
     &  ,KN,KK,KI,JJ,IG1,NL,LN,KL,NH,NC,NBEGIN,MBEGIN,MN,NEND,NN,MEND
      REAL*8 :: GANGLE,DGDTHET,DTEMP,GANGLE1,DGDTHET1
     &,RAD,ATOR,DATORI,DATORJ,SIJ,RSQIJ
     &,FC,DFC,XX,PX,EXX,DCTDIJ,DCTDIK,GS,XTEMP,GFX
     &,XSIJ,SSUMK,CONK,SDALIK,QI,SS
     &,ALI,DALI,DALDIK,RSQ3,RSQ2,RR,COSTH
     &,XSJI,SSUML,CONL,QJ,SDALJL
     &,ALJ,DALJ,DALDJL,DCTDJK,DCTDJI,DCTDJL,EXNIJ,DIJ,BIJ
     &,DBDZJ,VATT,DRADI,DRADJ,DRDC,CONJUG,XNT1,XNT2
     &,BTOT,DBTORI,DBTORJ,DBTORC,BTOR,DCTDIL,EXNJI,DJI,BJI,DBZI
     &,DATORC,SINK2,RCK,FCK,DFCK,SINL2,RCL,FCL,DFCL
     &,DT1DIK,DT1DJK,DT1DJL,DT1DIJ,CRKX,CRLX,CRKY,CRLY
     &,DT1DIL,CRKZ,CRLZ,CW,BT,AA,AAA1,AT2,RP1,RP2,RP3,RP4,RP5
     &,REP,VDBDI,VDBDJ,VDRDC,VDRDI,VDRDJ,RP,DWR,DBDZI,DDR

      REAL*8 :: XK(250,3),XL(250,3),XSIK(250)
     &,XSJK(250),XSIL(250),XSJL(250),CJ(3),CK(3),CL(3),RK(3),RL(3)
     &,DT2DIK(3),DT2DJL(3),DT2DIJ(3)
     &,DEXNI(2),DEXNJ(2),XNI(RTYPES),XNJ(RTYPES)
     &,CFUNI(ncc),CFUNJ(ncc),DCFUNI(ncc),DCFUNJ(ncc)
     &,COSK(250),COSL(250),SINK(250),SINL(250)
     &,DCTJK(250),DCTIJ(250),DCTIK(250),DCTIL(250),DCTJI(250)
     &,DCTJL(250)

C
C
C Sum over bonds between atoms I and J
C
      DO 40 I=1,NP
           JBEGIN=NABORS(I)
           JEND=NABORS(I+1)-1
           IF(JBEGIN.GT.JEND) GO TO 40
           KI=KTYPE(I)
C
           DO 30 J=JBEGIN,JEND
                IF(LCHECK(J).ne.1) GO TO 30
                JN=LIST(J)
                IF(I.GE.JN) GO TO 30
                DO 401 MM=1,3
                CJ(MM)=COR(J,MM)
401             CONTINUE
                SIJ=RCOR(J)
                RSQIJ=SIJ*SIJ
                KJ=KTYPE(JN)
                KIKJ=KI+KJ
C
C I side of bond
C
                NK=0
                XSIJ=0.0d0
                SSUMK=0.0d0
                CONK=0.0D0
                XNI(1)=XHC(I,1)
                XNI(2)=XHC(I,2)
                XNI(KJ)=XNI(KJ)-WW(J)
                QI=XNI(1)+XNI(2)-2.0D0
                SDALIK=0.0D0
                IF(JBEGIN.EQ.JEND) GO TO 21
C
                DO 20 K=JBEGIN,JEND
                     ALI=0.0D0
                     DALI=0.0D0
                     DALDIK=0.0D0
                     IF(K.EQ.J) GO TO 20
                     IF(LCHECK(K).ne.1) GO TO 20
                     KN=LIST(K)
                     KK=KTYPE(KN)
                     NK=NK+1
                     S3=RCOR(K)
                     RSQ3=S3*S3
                     RSQ2=0.0D0
                     DO 402 MM=1,3
                          XK(NK,MM)=COR(K,MM)-CJ(MM)
                          RSQ2=RSQ2+XK(NK,MM)*XK(NK,MM)
402                  CONTINUE
                     SS=2.0d0*SIJ*S3
                     RR=RSQIJ-RSQ3
                     COSTH=(RSQIJ+RSQ3-RSQ2)/SS
                     IF(COSTH.GT.1.0D0) COSTH=1.0D0
                     IF(COSTH.LT.-1.0D0) COSTH=-1.0D0
                     COSK(NK)=COSTH
                     SINK(NK)=SQRT(1.0D0-COSTH*COSTH)
                     IF(ACOS(COSTH).GT.PI) SINK(NK)=-SINK(NK)
                     IG=IGC(INT(-COSTH*12.0D0)+13)
                     IF(KI.EQ.1) THEN
                          IF(IG.NE.4) THEN
                               GANGLE=SPGC(1,IG)+SPGC(2,IG)*COSTH
                               DGDTHET=SPGC(2,IG)
                               DO 45 JJ=3,6
                                    GANGLE=GANGLE+SPGC(JJ,IG)
     &                                           *(COSTH**(JJ-1))
                                    DGDTHET=DGDTHET+SPGC(JJ,IG)*
     &                                     (JJ-1)*(COSTH**(JJ-2))
45                             CONTINUE
                          ELSE
                               ALI=0.0D0
                               DALI=0.0D0
                               IF(QI.LT.XQM) THEN
C####
c                                    QQE=QI-2.0D0
c                                    QQM=QQE/(QI-XQM)
c                                    ALI=ATT*EXP(0.05D0*QQE*QQM)
c                                    DALI=ALI*0.05D0*QQM*(2.0D0-QQM)
C####
                                    ALI=1.0D0
                                    IF(QI.GT.ATT) THEN
                                         DTEMP=PQ*(QI-ATT)
                                         ALI=(1.0D0+COS(DTEMP))/2.0D0
                                         DALI=-PQ/2.0D0*SIN(DTEMP)
                                    ENDIF
C####
                               ENDIF
                               GANGLE=SPGC(1,IG)+SPGC(2,IG)*COSTH
                               DGDTHET=SPGC(2,IG)
                               IG1=IG+1
                               GANGLE1=SPGC(1,IG1)+SPGC(2,IG1)*COSTH
                               DGDTHET1=SPGC(2,IG1)
                               DO 545 JJ=3,6
                                    GANGLE=GANGLE+SPGC(JJ,IG)
     &                                           *(COSTH**(JJ-1))
                                    DGDTHET=DGDTHET+SPGC(JJ,IG)*
     &                                     (JJ-1)*(COSTH**(JJ-2))
C
                                    GANGLE1=GANGLE1+SPGC(JJ,IG1)
     &                                           *(COSTH**(JJ-1))
                                    DGDTHET1=DGDTHET1+SPGC(JJ,IG1)*
     &                                     (JJ-1)*(COSTH**(JJ-2))
545                            CONTINUE
                               DALDIK=DALI*(GANGLE1-GANGLE)
                               GANGLE=GANGLE+ALI*(GANGLE1-GANGLE)
                               DGDTHET=DGDTHET+ALI*(DGDTHET1-DGDTHET)
                          ENDIF
                     ELSE
                     IG=IGH(INT(-COSTH*12.0D0)+13)
                     GANGLE=SPGH(1,IG)+SPGH(2,IG)*COSTH
                     DGDTHET=SPGH(2,IG)
                     DO 46 JJ=3,6
                          GANGLE=GANGLE+SPGH(JJ,IG)*(COSTH**(JJ-1))
                          DGDTHET=DGDTHET+SPGH(JJ,IG)*
     &                             (JJ-1)*(COSTH**(JJ-2))
46                   CONTINUE
                     ENDIF
                     FC=WW(K)
                     DFC=DWW(K)
                     CFUNI(NK)=0.0d0
                     DCFUNI(NK)=0.0d0
                     IF(KK.EQ.1) THEN
                          XX=XHC(KN,1)+XHC(KN,2)-FC-2.0d0
                          IF(XX.LT.3.0d0) THEN
                               IF(XX.LE.2.0d0) THEN
                                    CFUNI(NK)=1.0d0
                               ELSE
                                    PX=PI*(XX-2.0d0)
                                    CFUNI(NK)=(1.0+COS(PX))/2.0d0
                                    DCFUNI(NK)=-FC*SIN(PX)*PI/2.0d0
                               ENDIF
                          ENDIF
                     ENDIF
                     CONK=CONK+FC*CFUNI(NK)
C
                     IF(XDB(KI,KJ,KK).NE.0.0D0) THEN
                          EXX=REG(KI,KJ,KK)
     &                            *EXP(XDB(KI,KJ,KK)*(SIJ-S3))
                     ELSE
                          EXX=1.0
                     ENDIF
C
                     DCTDJK=-2.0d0/SS
                     DCTDIJ=(RR+RSQ2)/(SS*RSQIJ)
                     DCTDIK=(-RR+RSQ2)/(SS*RSQ3)
                     DCTJK(NK)=DCTDJK
                     DCTIJ(NK)=DCTDIJ
                     DCTIK(NK)=DCTDIK
                     GS=GANGLE*EXX
                     SSUMK=SSUMK + FC*GS
                     XTEMP=FC*EXX*DGDTHET
                     GFX=GS*FC*XDB(KI,KJ,KK)
                     XSIJ=XSIJ + XTEMP*DCTDIJ+GFX/SIJ
                     XSIK(NK)=(GS*DFC-GFX)/S3+XTEMP*DCTDIK
                     SDALIK=SDALIK+EXX*FC*DALDIK
                     XSJK(NK) = XTEMP*DCTDJK
20              CONTINUE
21              CONTINUE
C
C J side of bond
C
                NL=0
                XSJI=0.0d0
                SSUML=0.0d0
                CONL=0.0d0
                LBEGIN=NABORS(JN)
                LEND=NABORS(JN+1)-1
                XNJ(1)=XHC(JN,1)
                XNJ(2)=XHC(JN,2)
                XNJ(KI)=XNJ(KI)-WW(J)
                QJ=XNJ(1)+XNJ(2)-2.0D0
                SDALJL=0.0D0
                IF(LBEGIN.EQ.LEND) GO TO 11
C
                DO 10 L=LBEGIN,LEND
                     ALJ=0.0D0
                     DALJ=0.0D0
                     DALDJL=0.0D0
                     LN=LIST(L)
                     IF(LN.EQ.I) GO TO 10
                     IF(LCHECK(L).ne.1) GO TO 10
                     KL=KTYPE(LN)
                     NL=NL+1
                     S3=RCOR(L)
                     RSQ3=S3*S3
                     RSQ2=0.0D0
                     DO 403 MM=1,3
                          XL(NL,MM)=COR(L,MM)+CJ(MM)
                          RSQ2=RSQ2+XL(NL,MM)*XL(NL,MM)
403                  CONTINUE
                     SS=2.0d0*SIJ*S3
                     RR=RSQIJ-RSQ3
                     COSTH=(RSQIJ+RSQ3-RSQ2)/SS
                     IF(COSTH.GT.1.0D0) COSTH=1.0D0
                     IF(COSTH.LT.-1.0D0) COSTH=-1.0D0
                     COSL(NL)=COSTH
                     SINL(NL)=SQRT(1.0D0-COSTH*COSTH)
                     IF(ACOS(COSTH).GT.PI) SINL(NL)=-SINL(NL)
                     IF(KJ.EQ.1) THEN
                          IG=IGC(INT(-COSTH*12.0D0)+13)
                          IF(IG.NE.4) THEN
                               GANGLE=SPGC(1,IG)+SPGC(2,IG)*COSTH
                               DGDTHET=SPGC(2,IG)
                               DO 47 JJ=3,6
                                    GANGLE=GANGLE+SPGC(JJ,IG)
     &                                       *(COSTH**(JJ-1))
                                    DGDTHET=DGDTHET+SPGC(JJ,IG)
     &                                      *(JJ-1)*(COSTH**(JJ-2))
47                             CONTINUE
                          ELSE
                               ALJ=0.0D0
                               DALJ=0.0D0
                               IF(QJ.LT.XQM) THEN
C###
c                                    QQE=QJ-2.0D0
c                                    QQM=QQE/(QJ-XQM)
c                                    ALJ=ATT*EXP(0.05D0*QQE*QQM)
c                                    DALJ=ALJ*0.05D0*QQM*(2.0D0-QQM)
                                    ALJ=1.0D0
                                    IF(QJ.GT.ATT) THEN
                                         DTEMP=PQ*(QJ-ATT)
                                         ALJ=(1.0D0+COS(DTEMP))/2.0D0
                                         DALJ=-PQ/2.0D0*SIN(DTEMP)
                                    ENDIF
C####
                               ENDIF
                               GANGLE=SPGC(1,IG)+SPGC(2,IG)*COSTH
                               DGDTHET=SPGC(2,IG)
                               IG1=IG+1
                               GANGLE1=SPGC(1,IG1)+SPGC(2,IG1)*COSTH
                               DGDTHET1=SPGC(2,IG1)
                               DO 546 JJ=3,6
                                    GANGLE=GANGLE+SPGC(JJ,IG)
     &                                           *(COSTH**(JJ-1))
                                    DGDTHET=DGDTHET+SPGC(JJ,IG)*
     &                                     (JJ-1)*(COSTH**(JJ-2))
C
                                    GANGLE1=GANGLE1+SPGC(JJ,IG1)
     &                                           *(COSTH**(JJ-1))
                                    DGDTHET1=DGDTHET1+SPGC(JJ,IG1)*
     &                                     (JJ-1)*(COSTH**(JJ-2))
546                            CONTINUE
                               DALDJL=DALJ*(GANGLE1-GANGLE)
                               GANGLE=GANGLE+ALJ*(GANGLE1-GANGLE)
                               DGDTHET=DGDTHET+ALJ*(DGDTHET1-DGDTHET)
                          ENDIF
                     ELSE
                     IG=IGH(INT(-COSTH*12.0D0)+13)
                     GANGLE=SPGH(1,IG)+SPGH(2,IG)*COSTH
                     DGDTHET=SPGH(2,IG)
                     DO 48 JJ=3,6
                          GANGLE=GANGLE+SPGH(JJ,IG)*(COSTH**(JJ-1))
                          DGDTHET=DGDTHET+SPGH(JJ,IG)*(JJ-1)
     &                             *(COSTH**(JJ-2))
48                   CONTINUE
                     ENDIF
                     FC=WW(L)
                     DFC=DWW(L)
                     CFUNJ(NL)=0.0d0
                     DCFUNJ(NL)=0.0d0
c                     IF((KL+KIKJ).EQ.3) THEN
                     IF(KL.EQ.1) THEN
                          XX=XHC(LN,1)+XHC(LN,2)-FC-2.0d0
                          IF(XX.LT.3.0d0) THEN
                               IF(XX.LE.2.0d0) THEN
                                    CFUNJ(NL)=1.0d0
                               ELSE
                                    PX=PI*(XX-2.0d0)
                                    CFUNJ(NL)=(1.0+COS(PX))/2.0d0
                                    DCFUNJ(NL)=-FC*SIN(PX)*PI/2.0d0
                               ENDIF
                          ENDIF
                     ENDIF
                     CONL=CONL+FC*CFUNJ(NL)
C
                     IF(XDB(KJ,KI,KL).NE.0.0d0) THEN
                          EXX=REG(KJ,KI,KL)
     &                            *EXP(XDB(KJ,KI,KL)*(SIJ-S3))
                     ELSE
                          EXX=1.0d0
                     ENDIF
C
                     DCTDIL=-2.0d0/SS
                     DCTDJI=(RR+RSQ2)/(SS*RSQIJ)
                     DCTDJL=(-RR+RSQ2)/(SS*RSQ3)
                     DCTIL(NL)=DCTDIL
                     DCTJI(NL)=DCTDJI
                     DCTJL(NL)=DCTDJL
                     GS=GANGLE*EXX
                     SSUML=SSUML + FC*GS
                     XTEMP=FC*EXX*DGDTHET
                     GFX=GS*FC*XDB(KJ,KI,KL)
                     XSJI=XSJI + XTEMP*DCTDJI + GFX/SIJ
                     XSJL(NL)=(GS*DFC-GFX)/S3+XTEMP*DCTDJL
                     SDALJL=SDALJL+EXX*FC*DALDJL
                     XSIL(NL) = XTEMP*DCTDIL
10              CONTINUE
11              CONTINUE
C
                EXNIJ=0.0d0
                DEXNI(1)=0.0d0
                DEXNI(2)=0.0d0
                IF(KI.EQ.1) THEN
                      NH=INT(XNI(2)+1.0D-12)
                      NC=INT(XNI(1)+1.0D-12)
                      IF((ABS(FLOAT(NH)-XNI(2)).GT.1.0d-8).OR.
     &                   (ABS(FLOAT(NC)-XNI(1)).GT.1.0d-8)) THEN
                           CALL BCUINT(KI,KJ,XNI(2),XNI(1),NH,NC,
     &                                   EXNIJ,DEXNI(2),DEXNI(1))
                      ELSE
                           EXNIJ=XH(KJ,NH,NC)
                           DEXNI(2)=XH1(KJ,NH,NC)
                           DEXNI(1)=XH2(KJ,NH,NC)
                      ENDIF
                ENDIF
C
                EXNJI=0.0d0
                DEXNJ(1)=0.0d0
                DEXNJ(2)=0.0d0
                IF(KJ.EQ.1) THEN
                      NH=INT(XNJ(2)+1.0D-12)
                      NC=INT(XNJ(1)+1.0D-12)
                      IF((ABS(FLOAT(NH)-XNJ(2)).GT.1.0d-8).OR.
     &                   (ABS(FLOAT(NC)-XNJ(1)).GT.1.0d-8)) THEN
                           CALL BCUINT(KJ,KI,XNJ(2),XNJ(1),NH,NC,
     &                                   EXNJI,DEXNJ(2),DEXNJ(1))
                      ELSE
                           EXNJI=XH(KI,NH,NC)
                           DEXNJ(2)=XH1(KI,NH,NC)
                           DEXNJ(1)=XH2(KI,NH,NC)
                      ENDIF
                ENDIF
C
                DIJ=(1.0d0+EXNIJ+SSUMK)
                BIJ=DIJ**(-0.50D0)
                DJI=(1.0d0+EXNJI+SSUML)
                BJI=DJI**(-0.50D0)
                DBDZI=-0.50D0*BIJ/DIJ
                DBDZJ=-0.50D0*BJI/DJI
                VATT=EXX1(J)
                DRADI=0.0d0
                DRADJ=0.0d0
                DRDC=0.0d0
                CONJUG = 1.0D0 + (CONK**2) + (CONL**2)
                XNT1=XNI(1)+XNI(2)-1.0d0
                XNT2=XNJ(1)+XNJ(2)-1.0d0
                CALL RADIC(KI,KJ,XNT1,XNT2,CONJUG,RAD,DRADI,DRADJ,DRDC)
                BTOT=(BJI+BIJ+RAD)
C
C Dihedral terms
C
                IF(KIKJ.NE.NDIHED) GO TO 231
C
                     DBTORI=0.0d0
                     DBTORJ=0.0d0
                     DBTORC=0.0d0
                     BTOR=0.0d0
                     CALL TOR(XNT1,XNT2,CONJUG,ATOR,DATORI,DATORJ
     &                         ,DATORC)
C
                     IF(ABS(ATOR).LE.1.0d-08) GO TO 231
C
                     IF(JBEGIN.EQ.JEND) GO TO 230
                     IF(LBEGIN.EQ.LEND) GO TO 230
                     NK=0
                     DO 220 K=JBEGIN,JEND
                          IF(K.EQ.J) GO TO 220
                          IF(LCHECK(K).ne.1) GO TO 220
                          NK=NK+1
                          IF(ABS(SINK(NK)).LT.1.0D-01) GO TO 220
                          SINK2=SINK(NK)*SINK(NK)
                          KN=LIST(K)
                          DO 404 MM=1,3
                               CK(MM)=COR(K,MM)
404                       CONTINUE
                          RCK=RCOR(K)
C
                          IF(KTYPE(KN).EQ.2) THEN
                          FCK=1.0D0
                          DFCK=0.0D0
                          IF(RCK.GE.1.60D0) GO TO 220
                          IF(RCK.GE.1.30D0) THEN
                               DTEMP=PIDT*(RCK-1.30D0)
                               FCK=(1.0d0+COS(DTEMP))/2.0d0
                               DFCK=-PIDT/2.0d0*SIN(DTEMP)
                          ENDIF
                          ELSE
                          FCK=WW(K)
                          DFCK=DWW(K)
                          ENDIF
                          NL=0
                               DO 210 L=LBEGIN,LEND
                                    LN=LIST(L)
                                    IF(LN.EQ.I) GO TO 210
                                    IF(LCHECK(L).ne.1) GO TO 210
                                    NL=NL+1
                                    IF(ABS(SINL(NL)).LT.1.0D-01)
     &                                 GO TO 210
                                    SINL2=SINL(NL)*SINL(NL)
                                    DO 405 MM=1,3
                                         CL(MM)=COR(L,MM)
405                                 CONTINUE
                                    RCL=RCOR(L)
C
                                    IF(KTYPE(LN).EQ.2) THEN
                                    FCL=1.0D0
                                    DFCL=0.0D0
                                    IF(RCL.GE.1.60D0) GO TO 210
                                    IF(RCL.GE.1.30D0) THEN
                                         DTEMP=PIDT*(RCL-1.30D0)
                                         FCL=(1.0d0+COS(DTEMP))/2.0d0
                                         DFCL=-PIDT/2.0d0*SIN(DTEMP)
                                     ENDIF
                                     ELSE
                                     FCL=WW(L)
                                     DFCL=DWW(L)
                                     ENDIF
C
                                    T1=RCK*RCL*SIJ*SIJ
     &                                 *SINK(NK)*SINL(NL)
C
                                    DT1DIK=1.0/RCK/RCK
     &                              -DCTIK(NK)/SINK2*COSK(NK)
C
                                    DT1DJK=-DCTJK(NK)/SINK2
     &                                      *COSK(NK)
C
                                    DT1DJL=1.0/RCL/RCL
     &                              -DCTJL(NL)/SINL2*COSL(NL)
C
                                    DT1DIL=-DCTIL(NL)/SINL2
     &                                      *COSL(NL)
C
                                    DT1DIJ=2.0/SIJ/SIJ
     &                             -DCTIJ(NK)/SINK2*COSK(NK)
     &                             -DCTJI(NL)/SINL2*COSL(NL)
C
                                    CRKX=CK(2)*CJ(3)-CJ(2)*CK(3)
                                    CRLX=CJ(2)*CL(3)-CL(2)*CJ(3)
                                    CRKY=CK(3)*CJ(1)-CJ(3)*CK(1)
                                    CRLY=CJ(3)*CL(1)-CL(3)*CJ(1)
                                    CRKZ=CK(1)*CJ(2)-CJ(1)*CK(2)
                                    CRLZ=CJ(1)*CL(2)-CL(1)*CJ(2)
C
                                    T2=CRKX*CRLX+CRKY*CRLY
     &                                 +CRKZ*CRLZ
C
                                    CW=T2/T1
                                    BT=(1.0d0-CW*CW)
                                    BTOR=BTOR+BT*FCK*FCL
C   
                                    DT2DIK(1)=-CJ(3)*CRLY+CJ(2)*CRLZ
                                    DT2DIK(2)=-CJ(1)*CRLZ+CJ(3)*CRLX
                                    DT2DIK(3)=-CJ(2)*CRLX+CJ(1)*CRLY
C
                                    DT2DJL(1)=-CJ(2)*CRKZ+CJ(3)*CRKY
                                    DT2DJL(2)=-CJ(3)*CRKX+CJ(1)*CRKZ
                                    DT2DJL(3)=-CJ(1)*CRKY+CJ(2)*CRKX
C
                                    DT2DIJ(1)=CK(3)*CRLY-CL(3)*CRKY
     &                                     -CK(2)*CRLZ+CL(2)*CRKZ
C
                                    DT2DIJ(2)=CK(1)*CRLZ-CL(1)*CRKZ
     &                                     -CK(3)*CRLX+CL(3)*CRKX

                                    DT2DIJ(3)=CK(2)*CRLX-CL(2)*CRKX
     &                                     -CK(1)*CRLY+CL(1)*CRKY
                                    AA=-VATT*2.0d0*CW/T1*ATOR
     &                                  *FCL*FCK
                                    AAA1=VATT*BT*ATOR
                                    AT2=AA*T2
C
                                    RP1=-DT1DIJ*AT2
                                    RP2=-DT1DIK*AT2
     &                                 +AAA1*FCL*DFCK/RCK
                                    RP3=-DT1DJL*AT2
     &                                 +AAA1*FCK*DFCL/RCL
                                    RP4=-DT1DJK*AT2
                                    RP5=-DT1DIL*AT2
C
                                    DO 406 MM=1,3
                                         REP=RP1*CJ(MM) + AA*DT2DIJ(MM)
                                         RNP(I,MM)=RNP(I,MM) + REP
                                         RNP(JN,MM)=RNP(JN,MM) - REP
C
                                         REP=RP2*CK(MM) + AA*DT2DIK(MM)
                                         RNP(I,MM)=RNP(I,MM) + REP
                                         RNP(KN,MM)=RNP(KN,MM) - REP
C
                                         REP=RP3*CL(MM) + AA*DT2DJL(MM)
                                         RNP(JN,MM)=RNP(JN,MM) + REP
                                         RNP(LN,MM)=RNP(LN,MM) - REP
C
                                         REP=RP4*XK(NK,MM)
                                         RNP(JN,MM)=RNP(JN,MM) + REP
                                         RNP(KN,MM)=RNP(KN,MM) - REP
C
                                         REP=RP5*XL(NL,MM)
                                         RNP(I,MM)=RNP(I,MM) + REP
                                         RNP(LN,MM)=RNP(LN,MM) - REP
406                                 CONTINUE
210                            CONTINUE
220                       CONTINUE
230                  CONTINUE
C
                BTOT=BTOT+BTOR*ATOR
                DRADI=DRADI+DATORI*BTOR
                DRADJ=DRADJ+DATORJ*BTOR
                DRDC=DRDC+DATORC*BTOR
C
231             CONTINUE
C
C END DIHEDRAL FORCES
C
c********
                 TOTE=TOTE-BTOT*VATT

                eatom(i) = eatom(i) - btot*vatt/2.0d0
                eatom(jn)= eatom(jn) - btot*vatt/2.0d0
CC Print pair interaction

                VDBDI=VATT*DBDZI
                VDBDJ=VATT*DBDZJ
                VDRDC=VATT*DRDC
                VDRDI=VATT*DRADI
                VDRDJ=VATT*DRADJ
C
                RP= VDBDI*XSIJ + VDBDJ*XSJI + BTOT*DEXX1(J)
                DO 407 MM=1,3
                     REP=RP*CJ(MM)
                     RNP(I,MM)=RNP(I,MM)+REP
                     RNP(JN,MM)=RNP(JN,MM)-REP
407             CONTINUE
C
C Add many-body forces
C
C I side of bond
C
                IF(JBEGIN.EQ.JEND) GO TO 23
                NK=0
C
                DO 22 K=JBEGIN,JEND
                     IF(K.EQ.J) GO TO 22
                     IF(LCHECK(K).ne.1) GO TO 22
                     KN=LIST(K)
                     KK=KTYPE(KN)
                     DWR=DWW(K)/RCOR(K)
                     NK=NK+1
C
C First Neighbors
C
                     RP1=VDBDI*(XSIK(NK)+DWR*DEXNI(KK))
     &                  +DWR*(VDRDI+VDRDC*CFUNI(NK))
     &                  +VDBDI*DWR*SDALIK
                     RP2=VDBDI*XSJK(NK)
                     DO 408 MM=1,3
                          REP=RP1*COR(K,MM)
                          RNP(I,MM)=RNP(I,MM)+REP
                          RNP(KN,MM)=RNP(KN,MM)-REP
C
C Angular Forces
C
                          REP=RP2*XK(NK,MM)
                          RNP(JN,MM)=RNP(JN,MM)+REP
                          RNP(KN,MM)=RNP(KN,MM)-REP
408                  CONTINUE
C
C Second Neighbors via RADIC
C
                     DDR=VDRDC*DCFUNI(NK)*2.0D0*CONK
                     IF(DDR.EQ.0.0) GO TO 22
                     MBEGIN=NABORS(KN)
                     MEND=NABORS(KN+1)-1
                     IF(MBEGIN.GE.MEND) GO TO 22
C
                     DO 17 M=MBEGIN,MEND
                          IF(LCHECK(M).ne.1) GO TO 17
                          MN=LIST(M)
                          IF(MN.EQ.KN) GO TO 17
                          RP=DDR*DWW(M)/RCOR(M)
                          DO 409 MM=1,3
                               REP=RP*COR(M,MM)
                               RNP(KN,MM)=RNP(KN,MM)+REP
                               RNP(MN,MM)=RNP(MN,MM)-REP
409                       CONTINUE
17                   CONTINUE
22              CONTINUE
23              CONTINUE
C
C J side of bond
C
                IF(LBEGIN.EQ.LEND) GO TO 30
                NL=0
C
                DO 12 L=LBEGIN,LEND
                     LN=LIST(L)
                     IF(LN.EQ.I) GO TO 12
                     IF(LCHECK(L).ne.1) GO TO 12
                     KL=KTYPE(LN)
                     DWR=DWW(L)/RCOR(L)
                     NL=NL+1
C
C First Neighbors
C
                     RP1=VDBDJ*(XSJL(NL)+DWR*DEXNJ(KL))
     &                  +DWR*(VDRDJ+VDRDC*CFUNJ(NL))
     &                  +VDBDJ*DWR*SDALJL
                     RP2=VDBDJ*XSIL(NL)
                     DO 410 MM=1,3
                          REP=RP1*COR(L,MM)
                          RNP(JN,MM)=RNP(JN,MM)+REP
                          RNP(LN,MM)=RNP(LN,MM)-REP
C
C Angular Forces
C
                          REP=RP2*XL(NL,MM)
                          RNP(I,MM)=RNP(I,MM)+REP
                          RNP(LN,MM)=RNP(LN,MM)-REP
410                  CONTINUE
C
C Second Neighbors via RADIC
C
                     DDR=VDRDC*DCFUNJ(NL)*2.0D0*CONL
                     IF(DDR.EQ.0.0) GO TO 12
                     NBEGIN=NABORS(LN)
                     NEND=NABORS(LN+1)-1
                     IF(NBEGIN.GE.NEND) GO TO 12
C
                     DO 18 N=NBEGIN,NEND
                          IF(LCHECK(N).ne.1) GO TO 18
                          NN=LIST(N)
                          IF(NN.EQ.LN) GO TO 18
                          RP=DDR*DWW(N)/RCOR(N)
                          DO 411 MM=1,3
                          REP=RP*COR(N,MM)
                          RNP(LN,MM)=RNP(LN,MM)+REP
                          RNP(NN,MM)=RNP(NN,MM)-REP
411                  CONTINUE
18                   CONTINUE
12              CONTINUE
30         CONTINUE
40    CONTINUE
      RETURN
      END SUBROUTINE PIBONDCH

