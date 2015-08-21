      SUBROUTINE sili_germ
C
c     silicon: e = 4.6297255 eV/atom
c              rnn = 2.3521697 A
C
      USE dyn_array
      IMPLICIT REAL*8(A-H,O-Z)
c
      include 'common_files.inc'

      dimension xt(250,3) ,xslj(250),xsij(250)
      DO 30 L=1,NP
           JBEGIN=NABORS(L)
           JEND=NABORS(L+1)-1
           IF(JBEGIN.GT.JEND) GO TO 30
           kl = ktype(l)
           DO 40 I=JBEGIN,JEND
                IF(LCHECK(I).NE.2) GO TO 40
                NJ=0
                ki = ktype(i)
                XSLI=0.0
                SSUM=0.0
C  THIS IS THE IL TERM
                IN=LIST(I)
C  RSQ1 IS THE LI TERM
                S1=RCOR(I)
                RSQ1=S1*S1
                IF(JBEGIN.EQ.JEND) GO TO 51
C
                      DO 50 J=JBEGIN,JEND
                           IF(J.EQ.I) GO TO 50
                           IF(LCHECK(J).NE.2) GO TO 50
                           JN=LIST(J)
                           kj = ktype(jn)
                           NJ=NJ+1
C  RSQ3 IS THE LJ TERM
                           S3=RCOR(J)
                           RSQ3 = S3*S3
C   RSQ2 IS THE IJ TERM
                           RSQ2 = 0.0d0
                           do m=1,3
                                XT(NJ,m)=COR(J,M)-COR(I,M)
                                RSQ2 = RSQ2 + XT(NJ,m)*XT(NJ,m)
                           enddo
  640                      CONTINUE
                           SS=2.0*S1*S3
                           RR=RSQ1-RSQ3
                           COSTH=(RSQ1+RSQ3-RSQ2)/SS
                           AARG = (HDB(kl)+COSTH)
                           AARG2 = DDB2(kl) + AARG*AARG
                           GANGLE = 1.0d0 + CDB2(kl) * (1.0d0/DDB2(kl)
     $                             -1.0d0/AARG2)
C  DG / DCOSTHETA
                           DGDTHET=2.*CDB2(kl)*AARG/AARG2**2
C  DCOSTHETA / DRIJ  * (1./RIJ)
                           DCTDIJ = -2.0d0/SS
                           DCTDLI = (RR+RSQ2)/(SS*RSQ1)
                           DCTDLJ = (-RR+RSQ2)/(SS*RSQ3)
                           EXX = ADB(kl)* WW(J)
                           SSUM = SSUM +GANGLE*EXX
C  THE XSIJ  TERMS ARE GOING TO BE  DZ/DRIJ
                           XTEMP = EXX*DGDTHET
                           XSLI=XSLI + XTEMP*DCTDLI
C
                           XSLJ(NJ)= EXX * GANGLE * DWW(J) /S3
     &                              +XTEMP*DCTDLJ
C
                           XSIJ(NJ)= XTEMP*DCTDIJ
   50 CONTINUE
C
   51 CONTINUE
                      DLI=(1.0d0+(SSUM**XTN2(kl)))
                      BLI=DLI**(-XTN1(kl))
C
                      DBDZ=0.0
                      IF(SSUM.EQ.0.0) GO TO 33
                      DBDZ = -XTN1(kl)*BLI/DLI*XTN2(kl)
     &                         *(SSUM**(XTN2(kl)-1.0))
   33 CONTINUE
                      VATT =  EXX1(I)
                      VV=-VATT*BLI
                      TOTE=TOTE+VV
C                      eatom(l) = eatom(l) + vv/2.0d0
C                      eatom(in)= eatom(in) + vv/2.0d0

  650 CONTINUE
C  LI:
                      RP= VATT* DBDZ*XSLI + BLI * DEXX1(I)
                      do m=1,3
                           REP=RP*COR(I,m)
                           RNP(L,M)=RNP(L,M)+REP
                           RNP(IN,M)=RNP(IN,M)-REP
                      enddo
                      IF(JBEGIN.EQ.JEND) GO TO 40
C
                      NJ=0
C
                      DO 55 J=JBEGIN,JEND
                           IF(I.EQ.J) GO TO 55
                           IF(LCHECK(J).ne.2) GO TO 55
                           JN=LIST(J)
                           NJ=NJ+1
C  LJ:
                           RP= VATT*DBDZ*XSLJ(NJ)
                           do m=1,3
                                REP=RP*COR(J,M)
                                RNP(L,M)=RNP(L,M)+REP
                                RNP(JN,M)=RNP(JN,M)-REP
                           enddo
C  IJ:
                           RP= VATT*DBDZ*XSIJ(NJ)
                           do m=1,3
                                REP=RP*XT(NJ,M)
                                RNP(IN,M)=RNP(IN,M)+REP
                                RNP(JN,M)=RNP(JN,M)-REP
                           enddo
55                   CONTINUE
40              CONTINUE
30     CONTINUE
C
      RETURN
      END

