C Parallel pibond for CHO interactions 

      SUBROUTINE PIBONDCHO 

      USE MPIvars
      USE dyn_array      
      USE MDSTP
      USE STRUCTR
      USE PARAMS
      USE POTS
      USE SPECIF

      IMPLICIT none

      INTEGER :: I,J,K,L,M,N,JBEGIN,JEND,JN,MM,KJ,KIKJ,NK
     &,KN,KK,KI,JJ,IG1,NL,LN,KL,NH,NC,NBEGIN,MBEGIN,MN,NEND,NN,MEND
     &,NOX,NTCH
     &,NF,IG
      REAL*8 :: GANGLE,DGDTHET,DTEMP,GANGLE1,DGDTHET1
     &,RAD,ATOR,DATORI,DATORJ,SIJ,RSQIJ
     &,FC,DFC,XX,PX,EXX,DCTDIJ,DCTDIK,GS,XTEMP,GFX
     &,XSIJ,SSUMK,CONK,SDALIK,QI,SS
     &,ALI,DALI,DALDIK,RSQ3,RSQ2,RR,COSTH
     &,XSJI,SSUML,CONL,LBEGIN,LEND,QJ,SDALJL
     &,ALJ,DALJ,DALDJL,DCTDJK,DCTDJI,DCTDJL,EXNIJ,DIJ,BIJ
     &,DBDZJ,VATT,DRADI,DRADJ,DRDC,CONJUG,XNT1,XNT2
     &,BTOT,DBTORI,DBTORJ,DBTORC,BTOR,DCTDIL,EXNJI,DJI,BJI,DBZI
     &,DATORC,SINK2,RCK,FCK,DFCK,SINL2,RCL,FCL,DFCL
     &,DT1DIK,DT1DJK,DT1DJL,DT1DIJ,CRKX,CRLX,CRKY,CRLY
     &,DT1DIL,CRKZ,CRLZ,CW,BT,AA,AAA1,AT2,RP1,RP2,RP3,RP4,RP5
     &,REP,VDBDI,VDBDJ,VDRDC,VDRDI,VDRDJ,RP,DWR,DBDZI,DDR
     &,EXNI2J,EXNJ2I
     &,pno_min,pno_max,delt_po,pie0
C    & ,pij_min
     &,fijmid,dfijmid,pijtmp,fjimid,dfjimid

      REAL*8 :: XK(250,3),XL(250,3),XSIK(250)
     &,XSJK(250),XSIL(250),XSJL(250),CJ(3),CK(3),CL(3),RK(3),RL(3)
     &,DT2DIK(3),DT2DJL(3),DT2DIJ(3)
     &,DEXNI(RTYPES),DEXNJ(RTYPES),XNI(RTYPES),XNJ(RTYPES)
     &,CFUNI(ncc),CFUNJ(ncc),DCFUNI(ncc),DCFUNJ(ncc)
     &,COSK(250),COSL(250),SINK(250),SINL(250)
     &,DCTJK(250),DCTIJ(250),DCTIK(250),DCTIL(250),DCTJI(250)
     &,DCTJL(250)
     &,DEXNI2(RTYPES),DEXNJ2(RTYPES)
     &,dmax,dmin
C
C Find number hydrogens and carbons connected to each atom
C
      DO 500 I=1,NPtot_node
           JBEGIN=NABORS(I)
           JEND=NABORS(I+1)-1
C Initialize number of neighbor types  to 1
C probably should be zero with proper array allocations
C -travisk
           DO NN = 1,RTYPES
             XHC(I,NN)=1.0d0
           ENDDO
           IF(JBEGIN.GT.JEND) GO TO 500
C
           DO 490 J=JBEGIN,JEND
                IF(LCHECK(J).ne.1) GO TO 490
                JN=LIST(J)
                XHC(I,KTYPE_node(JN))=XHC(I,KTYPE_node(JN))+WW(J)
490        CONTINUE
C                IF (mynode.EQ.0) THEN
C              WRITE(297,122) i,NA(i),XHC(i,1),XHC(i,2),XHC(i,3)
C                ENDIF
C                IF (mynode.EQ.1) THEN
C              WRITE(298,122) i,NA(i),XHC(i,1),XHC(i,2),XHC(i,3)
C                ENDIF
C122             FORMAT (2(1x,I4),3(F12.8))
C
500   CONTINUE
C
C Sum over bonds between atoms I and J
C
      DO 40 I=1,NP_node
           JBEGIN=NABORS(I)
           JEND=NABORS(I+1)-1
           IF(JBEGIN.GT.JEND) GO TO 40
           KI=KTYPE_node(I)
C
           DO 30 J=JBEGIN,JEND
                IF(LCHECK(J).ne.1) GO TO 30
                JN=LIST(J)
                IF(NA(I).GE.NA(JN)) GO TO 30
                DO 401 MM=1,3
                CJ(MM)=COR(J,MM)
401             CONTINUE
                SIJ=RCOR(J)
                RSQIJ=SIJ*SIJ
                KJ=KTYPE_node(JN)
                KIKJ=KI+KJ

C
C I side of bond
C
                NK=0
                XSIJ=0.0d0
                SSUMK=0.0d0
                CONK=0.0D0
C -travisk
                XNI(:)=XHC(I,:)
                XNI(KJ)=XNI(KJ)-WW(J)
                QI =-1.0d0* FLOAT(RTYPES)
                DO NN = 1,RTYPES
                  QI = QI + XNI(NN)
                ENDDO
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
                     KK=KTYPE_node(KN)
c$$$                     WRITE(*,*) JN,I,KN
c$$$                     WRITE(*,*) 'TYPE'
c$$$                     WRITE(*,*) KJ,KI,KK
c$$$                     WRITE(*,*) '  '
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
                     GANGLE=0.0d0
                     DGDTHET=0.0d0
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
                     ELSEIF (KI.EQ.IHYD.AND.KJ.NE.IOXYGEN) THEN
                       IG=IGH(INT(-COSTH*12.0D0)+13)
                       GANGLE=SPGH(1,IG)+SPGH(2,IG)*COSTH
                       DGDTHET=SPGH(2,IG)
                       DO JJ=3,6
                          GANGLE=GANGLE+SPGH(JJ,IG)*(COSTH**(JJ-1))
                          DGDTHET=DGDTHET+SPGH(JJ,IG)*
     &                             (JJ-1)*(COSTH**(JJ-2))
                       ENDDO
C Oxygen angular function
C   -travisk
                     ELSEIF (KI.eq.ioxygen) THEN
                       GANGLE =  a_o0 + a_o1 * (a_o2 - costh)**2
                       DGDTHET =  2.d0 * a_o1 * (costh - a_o2)

                     ENDIF
                     FC=WW(K)
                     DFC=DWW(K)
                     CFUNI(NK)=0.0d0
                     DCFUNI(NK)=0.0d0
                     IF(KK.EQ.ICARB) THEN 
C.AND.KI.NE.IOXYGEN.AND.
C     &                  KJ.NE.IOXYGEN ) THEN
                          XX = - FC -1.0d0*FLOAT(RTYPES)
                          DO NN = 1,RTYPES
                            XX=XX + XHC(KN,NN)
                          ENDDO
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
            IF (PRNPR.and.NA(I).EQ.PR1.AND.NA(JN).EQ.PR2) THEN
              WRITE(23,*)  'CONK',CONK,FC,CFUNI(NK)
            ENDIF
C FROM EQU.8 
C E^/LAMBDA_IJK
                     IF(XDB(KI,KJ,KK).NE.0.0D0) THEN
                          EXX=REG(KI,KJ,KK)
     &                            *EXP(XDB(KI,KJ,KK)*(SIJ-S3))
                     ELSE
                          EXX=1.0d0
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
C -travisk
                XNJ(:)=XHC(JN,:)
                XNJ(KI)=XNJ(KI)-WW(J)
                QJ = -1.0d0*FLOAT(RTYPES)
                DO NN = 1,RTYPES
                  QJ=QJ + XNJ(NN)
                ENDDO
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
                     KL=KTYPE_node(LN)
                     NL=NL+1
c$$$                     WRITE(*,*) I,JN,LN
c$$$                     WRITE(*,*) 'TYPE'
c$$$                     WRITE(*,*) KI,KJ,KL
c$$$                     WRITE(*,*) '
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
                     GANGLE=0.0d0
                     DGDTHET=0.0d0
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
                                    ALJ=1.0D0
                                    IF(QJ.GT.ATT) THEN
                                         DTEMP=PQ*(QJ-ATT)
                                         ALJ=(1.0D0+COS(DTEMP))/2.0D0
                                         DALJ=-PQ/2.0D0*SIN(DTEMP)
                                    ENDIF
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
                     ELSEIF (KJ.eq.ihyd.AND.KI.NE.IOXYGEN) THEN
                       IG=IGH(INT(-COSTH*12.0D0)+13)
                       GANGLE=SPGH(1,IG)+SPGH(2,IG)*COSTH
                       DGDTHET=SPGH(2,IG)
                       DO JJ=3,6
                          GANGLE=GANGLE+SPGH(JJ,IG)*(COSTH**(JJ-1))
                          DGDTHET=DGDTHET+SPGH(JJ,IG)*(JJ-1)
     &                             *(COSTH**(JJ-2))
                       ENDDO
C Oxygen angular function
C   -travisk
                     ELSEIF (KJ.eq.ioxygen) THEN
                       GANGLE =  a_o0 + a_o1 * (a_o2 - costh)**2
                       DGDTHET =  2.d0 * a_o1 * (costh - a_o2)
                     ENDIF
                     FC=WW(L)
                     DFC=DWW(L)
                     CFUNJ(NL)=0.0d0
                     DCFUNJ(NL)=0.0d0
c                     IF((KL+KIKJ).EQ.3) THEN
                     IF(KL.EQ.icarb) THEN 
C !.and.KI.NE.ioxygen.and.
C     &                  KJ.NE.ioxygen ) THEN
                          XX = - FC -1.0d0*FLOAT(RTYPES)
            IF (PRNPR.and.NA(I).EQ.PR1.AND.NA(JN).EQ.PR2) THEN
               WRITE(23,*)  'XX0',XX,FC,FLOAT(RTYPES)
            ENDIF
                          DO NN = 1,RTYPES
                            XX=XX + XHC(LN,NN)
                          ENDDO
             IF (PRNPR.and.NA(I).EQ.PR1.AND.NA(JN).EQ.PR2) THEN
               WRITE(23,*)  'XXsum',XX
            ENDIF
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
            IF (PRNPR.and.NA(I).EQ.PR1.AND.NA(JN).EQ.PR2) THEN
              WRITE(23,*)  'CONL',CONL,FC,CFUNJ(NL),XX
               WRITE(23,*)  'XHC',XHC(NL,:)
             ENDIF
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
C Spline evaluation

C  Pij(I,NC,NH,NO)
C  -travisk
                EXNIJ=0.0d0
                DEXNI(:)=0.0d0
                EXNI2J = 0.0d0
                DEXNI2(:) = 0.0d0
C Carbon - Carbon spline PCC(NCO,NH)
C add oxygen as another carbon??
C same as AIREBO_CHO for now
C  -travis

C 
      pno_min = 1.01d0
      pno_max = 1.21d0
      delt_po = pno_max - pno_min
      pie0 = dacos(-1.d0)
C      pij_min =-0.39d0

                NH=INT(XNI(ihyd)+1.0D-12)
                NC=INT(XNI(icarb)+1.0D-12)
                NF=INT(XNI(iflor)+1.0D-12)
                IF (NF.NE.1) THEN
                  WRITE(*,*) 'NF not equal 1, something wrong'
                  CALL write_error
                ENDIF
                NOX=INT(XNI(ioxygen)+1.0D-12)
                NTCH =INT(XNI(ihyd)+XNI(icarb) +1.0D-12)
                IF(KI.NE.ioxygen.and.KJ.NE.ioxygen) THEN !1
                  IF(XNI(ioxygen).LE.pno_min ) THEN
                      IF( (ABS(FLOAT(NH)-XNI(ihyd)).GT.1.0d-8).OR.
     &                    (ABS( FLOAT(NC+NOX)-XNI(icarb)-XNI(ioxygen))
     &                    .GT.1.0d-8)) THEN 
                           CALL TCUINT(I,JN,KI,KJ,XNI(ihyd),XNI(icarb)
     &                     +XNI(ioxygen)-1.0d0,XNI(iflor),NH,NC+NOX-1
     &                  ,NF,EXNIJ,DEXNI(ihyd),DEXNI(icarb),DEXNI(iflor))
                      ELSE
                           EXNIJ=XH(KJ,NH,NC,NF)
                           DEXNI(ihyd)=XH1(KJ,NH,NC,NF)
                           DEXNI(icarb)=XH2(KJ,NH,NC,NF)
                           DEXNI(iflor)=XH3(KJ,NH,NC,NF)
                      ENDIF
                           DEXNI(ioxygen)=DEXNI(icarb)
                  ELSEIF ( XNI(ioxygen).GE.pno_max) THEN 
                    IF((ABS(FLOAT(NTCH)-XNI(icarb)-XNI(ihyd)).GT.1.0d-8)
     &             .or.(ABS(FLOAT(NOX) - XNI(ioxygen)).GT.1.0d-8)) THEN
                       CALL OXYBCUINT (KI,KJ,XNI(icarb)+XNI(ihyd)-1.0d0,
     &                               XNI(ioxygen),NTCH-1,NOX,EXNIJ
     &                              ,DEXNI(icarb),DEXNI(ioxygen) )
                        DEXNI(ihyd) = DEXNI(icarb)
c$$$                      IF ( EXNIJ.LT.pij_min) THEN
c$$$                         EXNIJ = pij_min
c$$$                      ENDIF
                    ELSE
                      EXNIJ  = XHO(KI,KJ,NTCH-1,NOX)
                      DEXNI(icarb) = XHO1(KI,KJ,NTCH-1,NOX)
                      DEXNI(ioxygen) = XHO2(KI,KJ,NTCH-1,NOX)
                      DEXNI(ihyd) = DEXNI(icarb)
                    ENDIF
                  ELSE !switching region 
                      IF( (ABS(FLOAT(NH)-XNI(ihyd)).GT.1.0d-8).OR.
     &                    (ABS( FLOAT(NC+NOX)-XNI(icarb)-XNI(ioxygen))
     &                    .GT.1.0d-8)) THEN 
                           CALL TCUINT(I,JN,KI,KJ,XNI(ihyd),XNI(icarb)
     &                     +XNI(ioxygen)-1.0d0,XNI(iflor),NH,NC+NOX-1
     &                  ,NF,EXNIJ,DEXNI(ihyd),DEXNI(icarb),DEXNI(iflor))
                      ELSE 
                           EXNIJ=XH(KJ,NH,NC,NF)
                           DEXNI(ihyd)=XH1(KJ,NH,NC,NF)
                           DEXNI(icarb)=XH2(KJ,NH,NC,NF)
                           DEXNI(iflor)=XH3(KJ,NH,NC,NF)
                      ENDIF
                           DEXNI(ioxygen)=DEXNI(icarb)
                    IF((ABS(FLOAT(NTCH)-XNI(icarb)-XNI(ihyd)).GT.1.0d-8)
     &             .or.(ABS(FLOAT(NOX) - XNI(ioxygen)).GT.1.0d-8)) THEN
                       CALL OXYBCUINT (KI,KJ,XNI(icarb)+XNI(ihyd)-1.0d0,
     &                               XNI(ioxygen),NTCH-1,NOX,EXNI2J
     &                              ,DEXNI2(icarb),DEXNI2(ioxygen) )
                        DEXNI2(ihyd) = DEXNI2(icarb)
c$$$                      IF ( EXNIJ.LT.pij_min) THEN
c$$$                         EXNI2J = pij_min
c$$$                      ENDIF
                    ELSE
                      EXNI2J  = XHO(KI,KJ,NTCH-1,NOX)
                      DEXNI2(icarb) = XHO1(KI,KJ,NTCH-1,NOX)
                      DEXNI2(ioxygen) = XHO2(KI,KJ,NTCH-1,NOX)
                      DEXNI2(ihyd) = DEXNI2(icarb)
                    ENDIF
C combine 2 spines
          fijmid = 0.5d0*(1.d0+cos(pie0*(XNI(ioxygen)-pno_min)
     &              /delt_po))
          dfijmid = - pie0 * 0.5d0 *sin(pie0*(XNI(ioxygen)-pno_min)
     &              /delt_po)/delt_po
          pijtmp = EXNIJ*fijmid + EXNI2J*(1-fijmid)
          DEXNI(icarb)=DEXNI(icarb)*fijmid+DEXNI2(icarb)*(1-fijmid)
          DEXNI(ihyd)=DEXNI(ihyd)*fijmid+DEXNI2(ihyd)*(1-fijmid)
          DEXNI(ioxygen) = DEXNI(ioxygen)*fijmid + 
     &                      DEXNI2(ioxygen)*(1-fijmid) +
     &                      (EXNIJ-EXNI2J)*dfijmid
          EXNIJ = pijtmp                                        
                  ENDIF
                ELSEIF ( KI.EQ.ioxygen.OR.KJ.EQ.ioxygen) THEN  
                  IF( (ABS(FLOAT(NTCH)-XNI(icarb)-XNI(ihyd)).GT.1.0d-8)
     &              .or.(ABS(FLOAT(NOX)-XNI(ioxygen)).GT.1.0d-8)) THEN
                   CALL OXYBCUINT(KI,KJ,XNI(icarb)+XNI(ihyd)-1.0d0,
     &                               XNI(ioxygen),NTCH-1,NOX,EXNIJ
     &                              ,DEXNI(icarb),DEXNI(ioxygen) )
                        DEXNI(ihyd) = DEXNI(icarb)
c$$$                      IF ( EXNIJ.LT.pij_min) THEN
c$$$                         EXNIJ = pij_min
c$$$                      ENDIF
                  ELSE
                      EXNIJ  = XHO(KI,KJ,NTCH-1,NOX)
                      DEXNI(icarb) = XHO1(KI,KJ,NTCH-1,NOX)
                      DEXNI(ioxygen) = XHO2(KI,KJ,NTCH-1,NOX)
                      DEXNI(ihyd) = DEXNI(icarb)
                  ENDIF
                ELSE
                  WRITE(*,*) 'Spline call error'
                  STOP
                ENDIF
C
C  Pji(I,NC,NH,NO)
C  -travisk
                EXNJI=0.0d0
                DEXNJ(:)=0.0d0
                EXNJ2I = 0.0d0
                DEXNJ2(:) = 0.0d0
C Carbon - Carbon spline PCC(NCO,NH)
C add oxygen as another carbon??
C same as AIREBO_CHO for now
C  -travis


                NH=INT(XNJ(ihyd)+1.0D-12)
                NC=INT(XNJ(icarb)+1.0D-12)
                NF=INT(XNJ(iflor)+1.0D-12)
                IF (NF.NE.1) THEN
                  WRITE(*,*) 'NF not equal 1, something wrong'
                  CALL write_error
                ENDIF
                NOX=INT(XNJ(ioxygen)+1.0D-12)
                NTCH =INT(XNJ(ihyd)+XNJ(icarb) +1.0D-12)
                IF(KJ.NE.ioxygen.and.KI.NE.ioxygen) THEN !1
                  IF(XNJ(ioxygen).LE.pno_min ) THEN
                      IF( (ABS(FLOAT(NH)-XNJ(ihyd)).GT.1.0d-8).OR.
     &                    (ABS( FLOAT(NC+NOX)-XNJ(icarb)-XNJ(ioxygen))
     &                    .GT.1.0d-8)) THEN 
                           CALL TCUINT(I,JN,KJ,KI,XNJ(ihyd),XNJ(icarb)
     &                 +XNJ(ioxygen)-1.0d0,XNJ(iflor) ,NH,NC+NOX - 1
     &                ,NF,EXNJI,DEXNJ(ihyd),DEXNJ(icarb), DEXNJ(iflor))
                      ELSE
                           EXNJI=XH(KI,NH,NC,NF)
                           DEXNJ(ihyd)=XH1(KI,NH,NC,NF)
                           DEXNJ(icarb)=XH2(KI,NH,NC,NF)
                           DEXNJ(iflor)=XH3(KI,NH,NC,NF)
                      ENDIF
                          DEXNJ(ioxygen)=DEXNJ(icarb)
                  ELSEIF ( XNJ(ioxygen).GE.pno_max) THEN 
                    IF((ABS(FLOAT(NTCH)-XNJ(icarb)-XNJ(ihyd)).GT.1.0d-8)
     &             .or.(ABS(FLOAT(NOX) - XNJ(ioxygen)).GT.1.0d-8)) THEN
                       CALL OXYBCUINT (KJ,KI,XNJ(icarb)+XNJ(ihyd)-1.0d0,
     &                               XNJ(ioxygen),NTCH-1,NOX,EXNJI
     &                              ,DEXNJ(icarb),DEXNJ(ioxygen) )
                        DEXNJ(ihyd) = DEXNJ(icarb)
C                      IF ( EXNJI.LT.pij_min) THEN
C                         EXNJI = pij_min
C                      ENDIF
                    ELSE
                      EXNJI  = XHO(KJ,KI,NTCH-1,NOX)
                      DEXNJ(icarb) = XHO1(KJ,KI,NTCH-1,NOX)
                      DEXNJ(ioxygen) = XHO2(KJ,KI,NTCH-1,NOX)
                      DEXNJ(ihyd) = DEXNJ(icarb)
                    ENDIF
                  ELSE !switching region 
                      IF( (ABS(FLOAT(NH)-XNJ(2)).GT.1.0d-8).OR.
     &                    (ABS( FLOAT(NC+NOX)-XNJ(icarb)-XNJ(ioxygen))
     &                    .GT.1.0d-8)) THEN 
                           CALL TCUINT(I,JN,KJ,KI,XNJ(ihyd),XNJ(icarb)
     &                 +XNJ(ioxygen)-1.0d0,XNJ(iflor) ,NH,NC+NOX - 1
     &                ,NF,EXNJI,DEXNJ(ihyd),DEXNJ(icarb), DEXNJ(iflor))
                      ELSE
                           EXNJI=XH(KI,NH,NC,NF)
                           DEXNJ(ihyd)=XH1(KI,NH,NC,NF)
                           DEXNJ(icarb)=XH2(KI,NH,NC,NF)
                           DEXNJ(iflor)=XH3(KI,NH,NC,NF)
                      ENDIF
                          DEXNJ(ioxygen)=DEXNJ(icarb)
                    IF((ABS(FLOAT(NTCH)-XNJ(icarb)-XNJ(ihyd)).GT.1.0d-8)
     &             .or.(ABS(FLOAT(NOX) - XNJ(ioxygen)).GT.1.0d-8)) THEN
                       CALL OXYBCUINT (KJ,KI,XNJ(icarb)+XNJ(ihyd)-1.0d0,
     &                               XNJ(ioxygen),NTCH-1,NOX,EXNJ2I
     &                              ,DEXNJ2(icarb),DEXNJ2(ioxygen) )
                        DEXNJ2(ihyd) = DEXNJ2(icarb)
C                      IF ( EXNJI.LT.pij_min) THEN
C                         EXNJ2I = pij_min
C                      ENDIF
                    ELSE
                      EXNJ2I  = XHO(KJ,KI,NTCH-1,NOX)
                      DEXNJ2(icarb) = XHO1(KJ,KI,NTCH-1,NOX)
                      DEXNJ2(ioxygen) = XHO2(KJ,KI,NTCH-1,NOX)
                      DEXNJ2(ihyd) = DEXNJ2(icarb)
                    ENDIF
C combine 2 spines
          fijmid = 0.5d0*(1.d0+cos(pie0*(XNJ(ioxygen)-pno_min)
     &              /delt_po))
          dfijmid = - pie0 * 0.5d0 *sin(pie0*(XNJ(ioxygen)-pno_min)
     &              /delt_po)/delt_po
          pijtmp = EXNJI*fijmid + EXNJ2I*(1-fijmid)
          DEXNJ(icarb)=DEXNJ(icarb)*fijmid+DEXNJ2(icarb)*(1-fijmid)
          DEXNJ(ihyd)=DEXNJ(ihyd)*fijmid+DEXNJ2(ihyd)*(1-fijmid)
          DEXNJ(ioxygen) = DEXNJ(ioxygen)*fijmid + 
     &                      DEXNJ2(ioxygen)*(1-fijmid) +
     &                      (EXNJI-EXNJ2I)*dfijmid
          EXNJI = pijtmp                                        
                  ENDIF
C Add sulfur bicubic spline
C  -travisk 
                ELSEIF( KJ.EQ.ioxygen.OR.KI.EQ.ioxygen) THEN
                  IF( (ABS(FLOAT(NTCH)-XNJ(icarb)-XNJ(ihyd)).GT.1.0d-8)
     &              .or.(ABS(FLOAT(NOX)-XNJ(ioxygen)).GT.1.0d-8)) THEN
                       CALL OXYBCUINT (KJ,KI,XNJ(icarb)+XNJ(ihyd)-1.0d0,
     &                               XNJ(ioxygen),NTCH-1,NOX,EXNJI
     &                              ,DEXNJ(icarb),DEXNJ(ioxygen) )
                        DEXNJ(ihyd) = DEXNJ(icarb)
c$$$                      IF ( EXNJI.LT.pij_min) THEN
c$$$                         EXNJI = pij_min
c$$$                      ENDIF
                    ELSE
                      EXNJI  = XHO(KJ,KI,NTCH-1,NOX)
                      DEXNJ(icarb) = XHO1(KJ,KI,NTCH-1,NOX)
                      DEXNJ(ioxygen) = XHO2(KJ,KI,NTCH-1,NOX)
                      DEXNJ(ihyd) = DEXNJ(icarb)
                    ENDIF
                ELSE
                  WRITE(*,*) 'Spline call error'
                  STOP
                ENDIF
C
C

                DIJ=(1.0d0+EXNIJ+SSUMK)
                BIJ=DIJ**(-0.50D0)
                DJI=(1.0d0+EXNJI+SSUML)
                BJI=DJI**(-0.50D0)
                DBDZI=-0.50D0*BIJ/DIJ
                DBDZJ=-0.50D0*BJI/DJI
                VATT=EXX1(J)

C Conjugate term
C Equ.14 of Brenner(2002)
C pi^{RC} = Fij(Nti,Ntj,Nconj)
C Nconj = 1 + (Ncarbon on I) + (Ncarbon on J)
                DRADI=0.0d0
                DRADJ=0.0d0
                DRDC=0.0d0
                CONJUG = 1.0D0 + (CONK**2) + (CONL**2)

C Loop over all neighbor types
C Nti = NN on I, Ntj = NN on J
C -travisk
                XNT1 = 1.0d0 - 1.0d0*FLOAT(RTYPES)
                XNT2 = 1.0d0 - 1.0d0*FLOAT(RTYPES)
                DO NN = 1,RTYPES
                  XNT1 = XNT1 + XNI(NN)
                  XNT2 = XNT2 + XNJ(NN)
                ENDDO

C Account for posible round off errors
                IF ( XNT1.GT.4.999d0) THEN
                  XNT1 = 4.999d0
                ELSE IF ( XNT1.LT.1.0d0) THEN
                  XNT1 = 1.0d0
                ENDIF
                IF ( XNT2.GT.4.999d0) THEN
                  XNT2 = 4.999d0
                ELSE IF ( XNT2.LT.1.0d0) THEN
                  XNT2 = 1.0d0
                ENDIF

            IF (PRNPR.and.NA(I).EQ.PR1.AND.NA(JN).EQ.PR2) THEN
               WRITE(23,*)  'RAD',XNT1,XNT2,CONJUG,RAD
            ENDIF


                IF (CONJUG.GT.9.999d0) THEN
                   CONJUG = 4.999d0
                ELSE IF (CONJUG.LT.1.0d0) THEN
                  CONJUG = 1.0d0
                ENDIF

                IF ( KI.EQ.ioxygen.OR.KJ.EQ.ioxygen) THEN
                  RAD = 0.0d0
                  DRADI = 0.0d0
                  DRADJ = 0.0d0
                  DRDC  = 0.0d0
                ELSE
C Switching function needs to be added
C -travisk
                 IF (XNI(ioxygen).GE.pno_max.OR.XNJ(ioxygen).GE.pno_max)
     &           THEN
                   RAD = 0.0d0
                   DRADI = 0.0d0
                   DRADJ = 0.0d0
                   DRDC  = 0.0d0
                 ELSEIF (XNI(ioxygen).LE.pno_min.AND.
     &                   XNJ(ioxygen).LE.pno_min) THEN  
                   CALL RADIC(KI,KJ,XNT1,XNT2,CONJUG,
     &                        RAD,DRADI,DRADJ,DRDC)
                 ELSE
                   CALL RADIC(KI,KJ,XNT1,XNT2,CONJUG,
     &                        RAD,DRADI,DRADJ,DRDC)
          fijmid = 0.5d0*(1.d0+cos(pie0*(XNI(ioxygen)-pno_min)
     +              /delt_po))
          dfijmid = - pie0 * 0.5d0 *sin(pie0*(XNI(ioxygen)-pno_min)
     +              /delt_po)/delt_po

          fjimid = 0.5d0*(1.d0+cos(pie0*(XNJ(ioxygen)-pno_min)
     +              /delt_po))
          dfjimid = - pie0 * 0.5d0 *sin(pie0*(XNJ(ioxygen)-pno_min)
     +              /delt_po)/delt_po

          rad = rad * fijmid * fjimid
          dradi = dradi * fijmid * fjimid
          dradj = dradj * fijmid * fjimid
          drdc = drdc * fijmid * fjimid
!dong, we use dpijdn for a temperory solution
          DEXNI(ioxygen) = DEXNI(ioxygen) + rad*dfijmid*fjimid
          DEXNJ(ioxygen) = DEXNJ(ioxygen) + rad*dfjimid*fijmid
                  ENDIF
                ENDIF



                BTOT=(BJI+BIJ+RAD)
C           IF (mynode.EQ.0) THEN
C           WRITE(295,192) mynode,NA(i),NA(jn),EXNIJ,SSUMK,EXNJI,SSUML
C           ENDIF
C           IF (mynode.EQ.1) THEN
C           WRITE(296,192) mynode,NA(i),NA(jn),EXNIJ,SSUMK,EXNJI,SSUML
C           ENDIF
C192        FORMAT (I2,2(1x,I4),4(F15.8))
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
ccccc                     IF(KTYPE_node(KN).EQ.2) THEN
                          IF(KTYPE_node(KN).EQ.2) THEN
                               dmin=1.3d0
                               dmax=1.6d0
c$$$                              else
c$$$                               dmin=1.7d0
c$$$                               dmax=2.0d0
c$$$                              endif
                            FCK=1.0D0
                            DFCK=0.0D0
                            IF(RCK.GE.dmax) GO TO 220
                            IF(RCK.GE.dmin) THEN
                               DTEMP=PIDT*(RCK-dmin)
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
ccccccccccccc                       IF(KTYPE_node(LN).EQ.2) THEN

C                                    IF(KTYPE_node(LN).GE.2) THEN
                                    IF (KTYPE_node(LN).EQ.2) THEN
                                      FCL=1.0D0
                                      DFCL=0.0D0
                                      dmin=1.3D0
                                      dmax=1.6d0
C                                     else
C                                      dmin=1.7D0
C                                      dmax=2.0d0
C                                     endif
                                      IF(RCL.GE.dmax) GO TO 210
                                      IF(RCL.GE.dmin) THEN
                                         DTEMP=PIDT*(RCL-dmin)
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
C           IF (mynode.EQ.0) THEN
C             WRITE(288,*) NA(i),NA(jn),NA(kn),NA(ln)
C             WRITE(288,199) RP1,RP2,RP3,RP4,RP5
C           ENDIF
C           IF (mynode.EQ.1) THEN
C             WRITE(290,*) NA(i),NA(jn),NA(kn),NA(ln)
C             WRITE(290,199) RP1,RP2,RP3,RP4,RP5
C           ENDIF
C199        FORMAT (5(F16.12))
                                DO 406 MM=1,3
                                  REP=RP1*CJ(MM) + AA*DT2DIJ(MM)
                                  RNP_node(MM,I)=RNP_node(MM,I) + REP
                                  RNP_node(MM,JN)=RNP_node(MM,JN) - REP
C
                                  REP=RP2*CK(MM) + AA*DT2DIK(MM)
                                  RNP_node(MM,I)=RNP_node(MM,I) + REP
                                  RNP_node(MM,KN)=RNP_node(MM,KN) - REP
C
                                  REP=RP3*CL(MM) + AA*DT2DJL(MM)
                                  RNP_node(MM,JN)=RNP_node(MM,JN) + REP
                                  RNP_node(MM,LN)=RNP_node(MM,LN) - REP
C
                                  REP=RP4*XK(NK,MM)
                                  RNP_node(MM,JN)=RNP_node(MM,JN) + REP
                                  RNP_node(MM,KN)=RNP_node(MM,KN) - REP
C
                                  REP=RP5*XL(NL,MM)
                                  RNP_node(MM,I)=RNP_node(MM,I) + REP
                                  RNP_node(MM,LN)=RNP_node(MM,LN) - REP
406                             CONTINUE
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

C                eatom(i) = eatom(i) - btot*vatt/2.0d0
C                eatom(jn)= eatom(jn) - btot*vatt/2.0d0
C

C
C Print pair interaction
C -travisk
C          IF ( mynode.EQ.0) THEN
            IF (PRNPR.and.NA(I).EQ.PR1.AND.NA(JN).EQ.PR2) THEN
              WRITE(23,*) 'I, J, BIJ, BJI,Gij,Gji, BTOT'
              WRITE(23,2302) NA(I),NA(JN),BIJ,BJI,SSUMK,SSUML,BTOT
              WRITE(23,*) GANGLE,EXNIJ,EXNJI
C             WRITE(23,*) XNI(ioxygen)
              PREN =  PREN - BTOT*VATT
            ENDIF
!
                VDBDI=VATT*DBDZI
                VDBDJ=VATT*DBDZJ
                VDRDC=VATT*DRDC
                VDRDI=VATT*DRADI
                VDRDJ=VATT*DRADJ
C
                RP= VDBDI*XSIJ + VDBDJ*XSJI + BTOT*DEXX1(J)
C           IF (mynode.EQ.0) THEN
C           WRITE(277,133) mynode,I,JN,NA(i),NA(jn),RP,BTOT,DEXX1(J)
C           ENDIF
C           IF (mynode.EQ.1) THEN
C           WRITE(278,133) mynode,i,JN,NA(i),NA(jn),RP,BTOT,DEXX1(J)
C           ENDIF
C133        FORMAT (I2,4(I4,2x),3(F15.8))
                DO 407 MM=1,3
                     REP=RP*CJ(MM)
                     RNP_node(MM,I)=RNP_node(MM,I)+REP
                     RNP_node(MM,JN)=RNP_node(MM,JN)-REP
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
                     KK=KTYPE_node(KN)
                     DWR=DWW(K)/RCOR(K)
                     NK=NK+1
C
C First Neighbors
C
                     RP1=VDBDI*(XSIK(NK)+DWR*DEXNI(KK))
     &                  +DWR*(VDRDI+VDRDC*CFUNI(NK))
     &                  +VDBDI*DWR*SDALIK
                     RP2=VDBDI*XSJK(NK)
C           IF (mynode.EQ.0) THEN
C             WRITE(279,121) NA(i),NA(jn),NA(kn),RP1,RP2
C           ENDIF
C           IF (mynode.EQ.1) THEN
C             WRITE(280,121) NA(i),NA(jn),NA(kn),RP1,RP2
C           ENDIF
C121        FORMAT (3(2x,I4),2(F15.8))
                     DO 408 MM=1,3
                          REP=RP1*COR(K,MM)
                          RNP_node(MM,I)=RNP_node(MM,I)+REP
                          RNP_node(MM,KN)=RNP_node(MM,KN)-REP
C
C Angular Forces
C
                          REP=RP2*XK(NK,MM)
                          RNP_node(MM,JN)=RNP_node(MM,JN)+REP
                          RNP_node(MM,KN)=RNP_node(MM,KN)-REP
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
C           IF (mynode.EQ.0) THEN
C             WRITE(281,122) NA(i),NA(jn),NA(kn),NA(mn),RP
C           ENDIF
C           IF (mynode.EQ.1) THEN
C             WRITE(282,122) NA(i),NA(jn),NA(kn),NA(mn),RP
C           ENDIF
C122        FORMAT (4(2x,I4),E23.15)
                          DO 409 MM=1,3
                               REP=RP*COR(M,MM)
                               RNP_node(MM,KN)=RNP_node(MM,KN)+REP
                               RNP_node(MM,MN)=RNP_node(MM,MN)-REP
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
                     KL=KTYPE_node(LN)
                     DWR=DWW(L)/RCOR(L)
                     NL=NL+1
C
C First Neighbors
C
                     RP1=VDBDJ*(XSJL(NL)+DWR*DEXNJ(KL))
     &                  +DWR*(VDRDJ+VDRDC*CFUNJ(NL))
     &                  +VDBDJ*DWR*SDALJL
                     RP2=VDBDJ*XSIL(NL)
C          IF (mynode.EQ.0) THEN
C             WRITE(283,103) I,JN,LN,NA(i),NA(jn),NA(ln),RP1,RP2
C           ENDIF
C           IF (mynode.EQ.1) THEN
C             WRITE(284,103) i,JN,LN,NA(i),NA(jn),NA(ln),RP1,RP2
C           ENDIF
C103        FORMAT (6(2x,I4),2(F15.8))
                     DO 410 MM=1,3
                          REP=RP1*COR(L,MM)
                          RNP_node(MM,JN)=RNP_node(MM,JN)+REP
                          RNP_node(MM,LN)=RNP_node(MM,LN)-REP
C
C Angular Forces
C
                          REP=RP2*XL(NL,MM)
                          RNP_node(MM,I)=RNP_node(MM,I)+REP
                          RNP_node(MM,LN)=RNP_node(MM,LN)-REP
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
C           IF (mynode.EQ.0) THEN
C             WRITE(285,106) NA(i),NA(jn),NA(ln),NA(nn),RP
C           ENDIF
C           IF (mynode.EQ.1) THEN
C             WRITE(286,106) NA(i),NA(jn),NA(ln),NA(nn),RP
C           ENDIF
C106        FORMAT (4(2x,I4),E25.15)
                          DO 411 MM=1,3
                          REP=RP*COR(N,MM)
                          RNP_node(MM,LN)=RNP_node(MM,LN)+REP
                          RNP_node(MM,NN)=RNP_node(MM,NN)-REP
411                  CONTINUE
18                   CONTINUE
12              CONTINUE
30         CONTINUE
40    CONTINUE
C      IF (mynode.EQ.0 ) WRITE(*,*) 'Pibond finished',RNP(2,2),tote

      RETURN
 2302 FORMAT(2I8,5F14.6)
      END

