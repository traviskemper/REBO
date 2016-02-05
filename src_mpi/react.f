!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Programer:
!     Travis Kemper
!     Department of Materials Science and Engineering
!     University of FLorida
!     traviskemper@ufl.edu
!
!     Version 1.0 1/213/11 T. W. Kemper
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Dedect and print reactions
!
      SUBROUTINE react
!
      USE MDSTP
      USE ANALYSIS
      USE STRUCTR
      USE POTS
      USE MPIvars
      USE BEAM 
      USE SPECIF
!
      IMPLICIT none
!
      INCLUDE 'mpif.h'
!
      INTEGER :: MXPRN
      PARAMETER ( MXPRN = 3) !Maximum # of molecules to print 
      INTEGER :: I,J,NBJ,JBEGIN,JEND,Ir,Ii,Oi,E,Ni,No,DELN
     & ,NBJi,JBEGINi,JENDi,NBJo,JBEGINo,JENDo,Io,Jo,Oo,Ji
     & ,MRCNT,MPCNT,ARCNT,APCNT,NNo,NNi
     & ,NMo,ECNTo,IMo
     & ,NMi,ECNTi,IMi
     & ,ELCNT(RTYPES)
     & ,RMOL(MMX),PMOL(MMX),NNLIST(NMX),NNCNT
     & ,KBEGINo,KENDo,NBKo,Ko
     & ,KBEGINi,KENDi,NBKi,Ki
     & ,KTI,KTJ,SPCNT,SP2CNT,SP3CNT
     & ,ID,D,M
      REAL*8 :: COM(3),MOLM,DP,DPTH,DPTHDIV,VSQR,MVEL(3)
      CHARACTER(30), DIMENSION(MMX):: RECTNT,PRODCT 
      CHARACTER(4), DIMENSION(2) :: PREF 
      CHARACTER(30*MXPRN) :: PRNP,PRNR
      LOGICAL :: FINDR,PBX,MMATCH
      CHARACTER(30) :: PREL
      CHARACTER(12*SDIV) :: BRDS,FRDS,EXDS
!
!
! begin debug
!      WRITE(*,*) 'Starting reaction check, previous recs: ',REC
! end debug       
      IF ( mynode.EQ.0 )THEN
!
      BRB(:) = 0
      BRBDP(:,:) = 0
      BRBGS(:) =0
      FRB(:) = 0
      FRBDP(:,:) = 0
      FRBGS(:) =0
      EXB(:) = 0
      EXBDP(:,:) = 0
      EXBGS(:) =0
!      
      SPCNT = 0
      SP2CNT = 0
      SP3CNT = 0
      ELCOMP(:,:) = 0 
      GASKE(:) = 0.0d0  
      DO 1001 I=1,NP
        JBEGINi=NABORSi(I)
        JENDi=NABORSi(I+1)-1
        NNi =JENDi - JBEGINi   + 1
        KI = KTYPEi(I)
        IF( KTYPEi(I).EQ.icarb ) THEN
          IF( NNi.EQ.4 ) THEN
            SP3CNT = SP3CNT + 1
          ELSEIF( NNi.EQ.3) THEN
            SP2CNT = SP2CNT + 1
          ELSEIF(NNi.EQ.2 ) THEN
            SPCNT= SPCNT + 1
          ENDIF 
        ENDIF 
!       Record atoms in gas and sub and gas KE
        IF ( R0i(DEPD,I).GT.SUBMN.AND.R0i(DEPD,I).LT.SUBCT ) THEN
          ELCOMP(1,KI) = ELCOMP(1,KI) +1
        ELSE
          ELCOMP(2,KI) = ELCOMP(2,KI) +1
          VSQR = 0.0d0
          DO M = 1,3
            VSQR = VSQR + R1i(M,I)*R1i(M,I)/DELTA/DELTA
          ENDDO 
          GASKE(KI) = GASKE(KI) + XMASS(KI)*VSQR*ECON
        ENDIF
        IF( RECTi(I).EQ.0 ) THEN 
          JBEGINo=NABORSo(I)
          JENDo=NABORSo(I+1)-1
          NNo = JENDo - JBEGINo  + 1
          DELN = NNo - NNi
          IF(KTYPEi(I).EQ.ihyd) THEN 
            IF ( NNo.GT.1 ) GOTO 1001
            IF ( NNi.GT.1 ) GOTO 1001
          ENDIF
          NNCNT = 0 
          IF ( NNo.GT.NNi ) THEN !broken bond          
!           Check for second nearest neighbors
            DO 2001  NBJ=JBEGINo,JENDo
                Jo = LISTo(NBJ)
                DO NBJi=JBEGINi,JENDi
                  Ji = LISTi(NBJi)
                  IF( Jo .EQ. Ji ) GOTO 2001
                ENDDO
! begin debug
c$$$           WRITE(27,*)' Step',LSTEP
c$$$           WRITE(27,*)' Initial atom',I,KTYPEi(I)
c$$$           WRITE(27,*)' N old/N new',NNo,NNi
c$$$           WRITE(27,*)' Neibor',Jo,'not found in new Nlist'
! end debug 
                NNCNT = NNCNT + 1
                NNLIST(NNCNT) = Jo            
                DO  NBJi=JBEGINi,JENDi
                  Ji = LISTi(NBJi)
                  KBEGINi = NABORSi(Ji)
                  KENDi = NABORSi(Ji+1)-1
                  DO NBKi = KBEGINi,KENDi
                    Ki = LISTi(NBKi)
                    IF( Jo.EQ.Ki ) GOTO 1001
                  ENDDO
                ENDDO 
 2001       CONTINUE
            RECTi(I) = 1
            REC = REC  + 1
            KTI = KTYPEi(I)
            KTJ = KTYPEi(Jo)
            Ir = I 
            MRCNT = 0
            MPCNT = 0
            ARCNT = 0
            APCNT = 0
            FINDR = .TRUE.
            DO WHILE ( FINDR) 
              RCLIST(Ir) = REC
              IF( .NOT.ASUBo(Ir) ) THEN
! begin debug
!             WRITE(*,*) 'Atom',Ir,'not part of sub'
! end debug
                NMo = IMOLo(Ir)
                MMATCH =.TRUE.
                DO No=1,MRCNT
                  IF(NMo.EQ.RMOL(No) ) MMATCH = .FALSE.
                ENDDO
                IF( MMATCH ) THEN
                  Io=MPNTo(NMo)
                  Jo=MPNTo(NMo+1)-1
                  ECNTo  = 0
                  DO 2002 Oo=Io,Jo
                    IMo = MOLSTo(Oo)
                    DO No=1,ARCNT
                     IF( IMo.EQ.ARLIST(No)) GOTO 2002
                    ENDDO
                    ECNTo = ECNTo + 1
                    ELISTi(ECNTo) = IMo
                    RECTi(IMo) = RECTi(I)
                    RCLIST(IMo) = REC
                    ARCNT = ARCNT + 1
                    ARLIST(ARCNT) = IMo
 2002             CONTINUE
                  IF( ECNTo.NE.0 ) THEN
                    MRCNT = MRCNT + 1
                    CALL PRNEL(ECNTo,PREL)
!                    CALL COMAS(ECNTo,COM,MOLM,MVEL)
                    RECTNT(MRCNT) = PREL
                    RMOL(MRCNT) = NMo
                  ENDIF
! begin debug
!              WRITE(*,*) 'Reactant Molecule ',NMo
!              WRITE(*,*) 'of atoms'(1:ECNTo)
! end debug                 
                ENDIF
              ELSE
! begin debug
!                WRITE(*,*) 'Atom',Ir,'part of sub'
! end debug
                RECTi(Ir) =RECTi(I)
                RCLIST(Ir) = REC
                ARCNT = ARCNT + 1
                ARLIST(ARCNT) = Ir
c$$$              DO 1202 NBJo=NABORSo(Ir),NABORSo(Ir+1)-1
c$$$                 Jo = LISTo(NBJo)
c$$$                 DO No=1,ARCNT
c$$$                   IF(Jo.EQ.ARLIST(No)) GOTO 1202
c$$$                 ENDDO
c$$$                 RECTi(Jo) =RECTi(I)
c$$$                 RCLIST(Jo) = REC
c$$$                 ARCNT = ARCNT + 1
c$$$                 ARLIST(ARCNT) = Jo
c$$$ 1202         CONTINUE
              ENDIF
              IF(.NOT.ASUBi(Ir) ) THEN
                NMi = IMOLi(Ir)
                MMATCH = .TRUE.
                DO Ni=1,MPCNT
                  IF(NMi.EQ.PMOL(Ni) ) MMATCH = .FALSE.
                ENDDO
                IF( MMATCH ) THEN
                  Ii=MPNTi(NMi)
                  Ji=MPNTi(NMi+1)-1
                  ECNTi  = 0
                  DO 2003 Oi=Ii,Ji
                    IMi = MOLSTi(Oi)
                    DO Ni=1,APCNT
                     IF( IMi.EQ.APLIST(Ni)) GOTO 2003
                    ENDDO
                    ECNTi = ECNTi + 1
                    ELISTi(ECNTi) = IMi
                    RECTi(IMi) =RECTi(I)
                    RCLIST(IMi) = REC
                    APCNT = APCNT + 1
                    APLIST(APCNT) = IMi
 2003             CONTINUE
                  IF ( ECNTi.GT.0 ) THEN 
                    MPCNT = MPCNT + 1
                    CALL PRNEL(ECNTi,PREL)
!                    CALL COMAS(ECNTi,COM,MOLM,MVEL)
                    PRODCT(MPCNT) = PREL
                    PMOL(MPCNT) = NMi
                  ENDIF
! begin debug
!                  WRITE(*,*) 'Product Molecule ',NMi
!                WRITE(*,*) 'of atoms',ELISTi(1:ECNTi)
c$$$
! end debug                 
                ENDIF
              ELSE
                RECTi(Ir) =RECTi(I)
                RCLIST(Ir) = REC
                APCNT = APCNT + 1
                APLIST(APCNT) = Ir
c$$$                DO 1201 NBJi=NABORSi(Ir),NABORSi(Ir+1)-1
c$$$                   Ji = LISTi(NBJo)
c$$$                   DO Ni=1,APCNT
c$$$                     IF(Ji.EQ.APLIST(Ni)) GOTO 1201
c$$$                   ENDDO
c$$$                   RECTi(Ji) =RECTi(I)
c$$$                   RCLIST(Ji) = REC
c$$$                   APCNT = APCNT + 1
c$$$                   APLIST(APCNT) = Ji
c$$$ 1201           CONTINUE 
              ENDIF
!             Check for unacounted atoms
              IF ( NNCNT.NE.0 ) THEN
                Ir = NNLIST(NNCNT)
! begin debug
!            WRITE(27,*)'   Check new n',NNCNT,Ir
! end debug                   
                NNCNT = NNCNT - 1
              ELSE
                FINDR = .FALSE.
                DO 2004 Ni = 1,APCNT
                  Ji = APLIST(Ni)
                  DO No = 1,ARCNT
                    Jo = ARLIST(No)
                    IF( Ji.EQ.Jo) GOTO 2004
                  ENDDO            
                  FINDR = .TRUE.
                  Ir = Ji
 2004           CONTINUE 
                DO 2005 No = 1,ARCNT
                  Jo = ARLIST(No)
                  DO Ni = 1,APCNT
                    Ji = APLIST(Ni)
                    IF( Ji.EQ.Jo) GOTO 2005
                  ENDDO            
                  FINDR  = .TRUE.
                  Ir = Jo
 2005           CONTINUE 
              ENDIF
            ENDDO            
!           All reactants and products found
    !       CALL PRSTR(ECNT,ELIST,REC,)
            IF ( MRCNT.LE.MXPRN )THEN
              WRITE( PRNR,*) RECTNT(1:MRCNT)
            ELSE
              WRITE( PRNR,*) '# of mol in reactents >MXPRN'
            ENDIF 
            IF ( MPCNT.LE.MXPRN )THEN
              WRITE( PRNP,*) PRODCT(1:MPCNT)
            ELSE
              WRITE( PRNP,*) '# of mol in products >MXPRN'
            ENDIF 
            WRITE(27,271)  TTIME,LSTEP,REC
     &       ,ATYPE(KTI),ATYPE(KTJ),I,Jo,MRCNT,PRNR,MPCNT,PRNP
!           Store bond depth profile information
            ID = BID(KTI,KTJ)
            DP = R0i(DEPD,I)
            IF ( DP.GT.SUBMN.AND. DP.LT.SUBCT ) THEN
              DPTH =  SUBCT - DP
              DPTHDIV = DPTH/DPDIV + 1.0d0
              D = INT(DPTHDIV) 
              BRB(ID) = BRB(ID) + 1
              BRBDP(D,ID) = BRBDP(D,ID) + 1
            ELSE
              BRBGS(ID) = BRBGS(ID) + 1
            ENDIF
          ELSEIF ( NNi.GT.NNo) THEN !new bond
!           Check for second nearest neighbors
            DO 3001 NBJi=JBEGINi,JENDi
              Ji = LISTi(NBJi)
              DO  NBJo=JBEGINo,JENDo
                Jo = LISTo(NBJo)
                IF ( Ji.EQ.Jo ) GOTO 3001
              ENDDO
! begin debug
c$$$           WRITE(27,*)' Step',LSTEP
c$$$           WRITE(27,*)' Initial atom',I,KTYPEi(I)
c$$$           WRITE(27,*)' N old/N new',NNo,NNi
c$$$           WRITE(27,*)' Neibor',Ji,'not found in old  Nlist'
! end debug                   
              NNCNT = NNCNT + 1
              NNLIST(NNCNT) = Ji       
              DO  NBJ=JBEGINo,JENDo
                Jo = LISTo(NBJ)
                KBEGINo = NABORSo(Jo)
                KENDo = NABORSo(Jo+1)-1
                DO NBKo = KBEGINo,KENDo
                  Ko = LISTo(NBKo)
                  IF( Ji.EQ.Ko ) GOTO 1001
                ENDDO
              ENDDO
 3001       CONTINUE
            RECTi(I) = 2 
            REC = REC  + 1
            KTI = KTYPEi(I)
            KTJ = KTYPEi(Ji)  
            Ir = I 
            MRCNT = 0
            MPCNT = 0
            ARCNT = 0
            APCNT = 0
            FINDR = .TRUE.
            DO WHILE ( FINDR) 
              RCLIST(Ir) = REC
              IF( .NOT.ASUBo(Ir) ) THEN
! begin debug
!          WRITE(*,*) 'Atom',Ir,'not part of sub'
! end debug
                NMo = IMOLo(Ir)
                MMATCH =.TRUE.
                DO No=1,MRCNT
                  IF(NMo.EQ.RMOL(No) ) MMATCH = .FALSE.
                ENDDO
                IF( MMATCH ) THEN
                  Io=MPNTo(NMo)
                  Jo=MPNTo(NMo+1)-1
                  ECNTo  = 0
                  DO 3002 Oo=Io,Jo
                    IMo = MOLSTo(Oo)
                    DO No=1,ARCNT
                     IF( IMo.EQ.ARLIST(No)) GOTO 3002
                    ENDDO
                    ECNTo = ECNTo + 1
                    ELISTi(ECNTo) = IMo
                    RECTi(IMo) = RECTi(I)
                    RCLIST(IMo) = REC
                    ARCNT = ARCNT + 1
                    ARLIST(ARCNT) = IMo
 3002             CONTINUE
                  IF( ECNTo.NE.0 ) THEN
                    MRCNT = MRCNT + 1
                    CALL PRNEL(ECNTo,PREL)
!                    CALL COMAS(ECNTo,COM,MOLM,MVEL)
                    RECTNT(MRCNT) = PREL
                    RMOL(MRCNT) = NMo
                  ENDIF
! begin debug
!              WRITE(*,*) 'Reactant Molecule ',NMo
!              WRITE(*,*) 'of atoms',ELISTi(1:ECNTo)
! end debug                 
                ENDIF
              ELSE
! begin debug
!                WRITE(*,*) 'Atom',Ir,'part of sub'
! end debug
                RECTi(Ir) =RECTi(I)
                RCLIST(Ir) = REC
                ARCNT = ARCNT + 1
                ARLIST(ARCNT) = Ir
c$$$              DO 1202 NBJo=NABORSo(Ir),NABORSo(Ir+1)-1
c$$$                 Jo = LISTo(NBJo)
c$$$                 DO No=1,ARCNT
c$$$                   IF(Jo.EQ.ARLIST(No)) GOTO 1202
c$$$                 ENDDO
c$$$                 RECTi(Jo) =RECTi(I)
c$$$                 RCLIST(Jo) = REC
c$$$                 ARCNT = ARCNT + 1
c$$$                 ARLIST(ARCNT) = Jo
c$$$ 1202         CONTINUE
              ENDIF
              IF(.NOT.ASUBi(Ir) ) THEN
                NMi = IMOLi(Ir)
                MMATCH = .TRUE.
                DO Ni=1,MPCNT
                  IF(NMi.EQ.PMOL(Ni) ) MMATCH = .FALSE.
                ENDDO
                IF( MMATCH ) THEN
                  Ii=MPNTi(NMi)
                  Ji=MPNTi(NMi+1)-1
                  ECNTi  = 0
                  DO 3003 Oi=Ii,Ji
                    IMi = MOLSTi(Oi)
                    DO Ni=1,APCNT
                     IF( IMi.EQ.APLIST(Ni)) GOTO 3003
                    ENDDO
                    ECNTi = ECNTi + 1
                    ELISTi(ECNTi) = IMi
                    RECTi(IMi) =RECTi(I)
                    RCLIST(IMi) = REC
                    APCNT = APCNT + 1
                    APLIST(APCNT) = IMi
 3003             CONTINUE
                  IF ( ECNTi.GT.0 ) THEN 
                    MPCNT = MPCNT + 1
                    CALL PRNEL(ECNTi,PREL)
!                    CALL COMAS(ECNTi,COM,MOLM,MVEL)
                    PRODCT(MPCNT) = PREL
                    PMOL(MPCNT) = NMi
                  ENDIF
! begin debug
!                  WRITE(*,*) 'Product Molecule ',NMi
!                WRITE(*,*) 'of atoms',ELISTi(1:ECNTi)
c$$$
! end debug                 
                ENDIF
              ELSE
                RECTi(Ir) =RECTi(I)
                RCLIST(Ir) = REC
                APCNT = APCNT + 1
                APLIST(APCNT) = Ir
c$$$                DO 1201 NBJi=NABORSi(Ir),NABORSi(Ir+1)-1
c$$$                   Ji = LISTi(NBJo)
c$$$                   DO Ni=1,APCNT
c$$$                     IF(Ji.EQ.APLIST(Ni)) GOTO 1201
c$$$                   ENDDO
c$$$                   RECTi(Ji) =RECTi(I)
c$$$                   RCLIST(Ji) = REC
c$$$                   APCNT = APCNT + 1
c$$$                   APLIST(APCNT) = Ji
c$$$ 1201           CONTINUE 
              ENDIF
!             Check for unacounted atoms
              IF ( NNCNT.NE.0 ) THEN
                Ir = NNLIST(NNCNT)
! begin debug
!            WRITE(27,*)'   Check new n',NNCNT,Ir
! end debug                   
                NNCNT = NNCNT - 1
              ELSE
                FINDR = .FALSE.
                DO 3004 Ni = 1,APCNT
                  Ji = APLIST(Ni)
                  DO No = 1,ARCNT
                    Jo = ARLIST(No)
                    IF( Ji.EQ.Jo) GOTO 3004
                  ENDDO            
                  FINDR = .TRUE.
                  Ir = Ji
 3004           CONTINUE 
                DO 3005 No = 1,ARCNT
                  Jo = ARLIST(No)
                  DO Ni = 1,APCNT
                    Ji = APLIST(Ni)
                    IF( Ji.EQ.Jo) GOTO 3005
                  ENDDO            
                  FINDR  = .TRUE.
                  Ir = Jo
 3005           CONTINUE 
              ENDIF
            ENDDO            
!           All reactants and products found
!           CALL PRSTR(ECNT,ELIST,REC,)
            IF ( MRCNT.LE.MXPRN )THEN
              WRITE( PRNR,*) RECTNT(1:MRCNT)
            ELSE
              WRITE( PRNR,*) '# of mol in reactents >MXPRN'
            ENDIF 
            IF ( MPCNT.LE.MXPRN )THEN
              WRITE( PRNP,*) PRODCT(1:MPCNT)
            ELSE
              WRITE( PRNP,*) '# of mol in products >MXPRN'
            ENDIF 
            WRITE(27,272)  TTIME,LSTEP,REC
     &        ,ATYPE(KTI),ATYPE(KTJ),I,Ji,MRCNT,PRNR,MPCNT,PRNP
!           Store bond depth profile information
            ID = BID(KTI,KTJ)
            DP = R0i(DEPD,I)
            IF ( DP.GT.SUBMN.AND. DP.LT.SUBCT ) THEN
              DPTH =  SUBCT - DP
              DPTHDIV = DPTH/DPDIV + 1.0d0
              D = INT(DPTHDIV) 
              FRB(ID) = FRB(ID) + 1
              FRBDP(D,ID) = FRBDP(D,ID) + 1
            ELSE
              FRBGS(ID) = FRBGS(ID) + 1
            ENDIF
!         Check for exchange
          ELSE
            DO 4001  NBJo=JBEGINo,JENDo
              Jo = LISTo(NBJo)
              DO NBJi=JBEGINi,JENDi
                Ji = LISTi(NBJi)
                IF( Jo .EQ. Ji ) GOTO 4001
              ENDDO
! begin debug
c$$$          WRITE(27,*)' Exchange at Step',LSTEP
c$$$          WRITE(27,*)' Initial atom',I,KTYPEi(I) 
c$$$          WRITE(27,*)' Neibor',Jo,'not found in new Nlist'
! end debug                   
              NNCNT = NNCNT + 1
              NNLIST(NNCNT) = Jo
              DO  NBJi=JBEGINi,JENDi
                  Ji = LISTi(NBJi)
                  KBEGINi = NABORSi(Ji)
                  KENDi = NABORSi(Ji+1)-1
                  DO NBKi = KBEGINi,KENDi
                    Ki = LISTi(NBKi)
                    IF( Jo.EQ.Ki ) GOTO 1001
                  ENDDO
              ENDDO 
 4001       CONTINUE
!           Check for second nearest neighbors
            DO 4011 NBJi=JBEGINi,JENDi
              Ji = LISTi(NBJi)
              DO  NBJo=JBEGINo,JENDo
                Jo = LISTo(NBJo)
                IF ( Ji.EQ.Jo ) GOTO 4011
              ENDDO
! begin debug
c$$$           WRITE(27,*)' Step',LSTEP
c$$$           WRITE(27,*)' Initial atom',I,KTYPEi(I)
c$$$           WRITE(27,*)' N old/N new',NNo,NNi
c$$$           WRITE(27,*)' Neibor',Ji,'not found in old  Nlist'
! end debug                   
              NNCNT = NNCNT + 1
              NNLIST(NNCNT) = Ji       
              DO  NBJ=JBEGINo,JENDo
                Jo = LISTo(NBJ)
                KBEGINo = NABORSo(Jo)
                KENDo = NABORSo(Jo+1)-1
                DO NBKo = KBEGINo,KENDo
                  Ko = LISTo(NBKo)
                  IF( Ji.EQ.Ko ) GOTO 1001
                ENDDO
              ENDDO
 4011       CONTINUE
            IF( NNCNT.EQ.0 ) GOTO 1001
            REC = REC  + 1
            KTI = KTYPEi(I)
!            WRITE(27,*) 'A ',ATYPE(KTI)
!     &     ,' atom has exhanged some nieghbors at step',LSTEP
!            WRITE(27,*)' N old/N new',NNo,NNi
            Ir = I 
            MRCNT = 0
            MPCNT = 0
            ARCNT = 0
            APCNT = 0
            FINDR = .TRUE.
            DO WHILE ( FINDR) 
              RCLIST(Ir) = REC
              IF( .NOT.ASUBo(Ir) ) THEN
! begin debug
!          WRITE(*,*) 'Atom',Ir,'not part of sub'
! end debug
                NMo = IMOLo(Ir)
                MMATCH =.TRUE.
                DO No=1,MRCNT
                  IF(NMo.EQ.RMOL(No) ) MMATCH = .FALSE.
                ENDDO
                IF( MMATCH ) THEN
                  Io=MPNTo(NMo)
                  Jo=MPNTo(NMo+1)-1
                  ECNTo  = 0
                  DO 4002 Oo=Io,Jo
                    IMo = MOLSTo(Oo)
                    DO No=1,ARCNT
                     IF( IMo.EQ.ARLIST(No)) GOTO 4002
                    ENDDO
                    ECNTo = ECNTo + 1
                    ELISTi(ECNTo) = IMo
                    RECTi(IMo) = RECTi(I)
                    RCLIST(IMo) = REC
                    ARCNT = ARCNT + 1
                    ARLIST(ARCNT) = IMo
 4002             CONTINUE
                  IF( ECNTo.NE.0 ) THEN
                    MRCNT = MRCNT + 1
                    CALL PRNEL(ECNTo,PREL)
!                    CALL COMAS(ECNTo,COM,MOLM,MVEL)
                    RECTNT(MRCNT) = PREL
                    RMOL(MRCNT) = NMo
                  ENDIF
! begin debug
!              WRITE(*,*) 'Reactant Molecule ',NMo
!              WRITE(*,*) 'of atoms',ELISTi(1:ECNTo)
! end debug                 
                ENDIF
              ELSE
! begin debug
!                WRITE(*,*) 'Atom',Ir,'part of sub'
! end debug
                RECTi(Ir) =RECTi(I)
                RCLIST(Ir) = REC
                ARCNT = ARCNT + 1
                ARLIST(ARCNT) = Ir
c$$$              DO 1202 NBJo=NABORSo(Ir),NABORSo(Ir+1)-1
c$$$                 Jo = LISTo(NBJo)
c$$$                 DO No=1,ARCNT
c$$$                   IF(Jo.EQ.ARLIST(No)) GOTO 1202
c$$$                 ENDDO
c$$$                 RECTi(Jo) =RECTi(I)
c$$$                 RCLIST(Jo) = REC
c$$$                 ARCNT = ARCNT + 1
c$$$                 ARLIST(ARCNT) = Jo
c$$$ 1202         CONTINUE
              ENDIF
              IF(.NOT.ASUBi(Ir) ) THEN
                NMi = IMOLi(Ir)
                MMATCH = .TRUE.
                DO Ni=1,MPCNT
                  IF(NMi.EQ.PMOL(Ni) ) MMATCH = .FALSE.
                ENDDO
                IF( MMATCH ) THEN
                  Ii=MPNTi(NMi)
                  Ji=MPNTi(NMi+1)-1
                  ECNTi  = 0
                  DO 4003 Oi=Ii,Ji
                    IMi = MOLSTi(Oi)
                    DO Ni=1,APCNT
                     IF( IMi.EQ.APLIST(Ni)) GOTO 4003
                    ENDDO
                    ECNTi = ECNTi + 1
                    ELISTi(ECNTi) = IMi
                    RECTi(IMi) =RECTi(I)
                    RCLIST(IMi) = REC
                    APCNT = APCNT + 1
                    APLIST(APCNT) = IMi
 4003             CONTINUE
                  IF ( ECNTi.GT.0 ) THEN 
                    MPCNT = MPCNT + 1
                    CALL PRNEL(ECNTi,PREL)
!                    CALL COMAS(ECNTi,COM,MOLM,MVEL)
                    PRODCT(MPCNT) = PREL
                    PMOL(MPCNT) = NMi
                  ENDIF
! begin debug
!                  WRITE(*,*) 'Product Molecule ',NMi
!                WRITE(*,*) 'of atoms',ELISTi(1:ECNTi)
c$$$
! end debug                 
                ENDIF
              ELSE
                RECTi(Ir) =RECTi(I)
                RCLIST(Ir) = REC
                APCNT = APCNT + 1
                APLIST(APCNT) = Ir
c$$$                DO 1201 NBJi=NABORSi(Ir),NABORSi(Ir+1)-1
c$$$                   Ji = LISTi(NBJo)
c$$$                   DO Ni=1,APCNT
c$$$                     IF(Ji.EQ.APLIST(Ni)) GOTO 1201
c$$$                   ENDDO
c$$$                   RECTi(Ji) =RECTi(I)
c$$$                   RCLIST(Ji) = REC
c$$$                   APCNT = APCNT + 1
c$$$                   APLIST(APCNT) = Ji
c$$$ 1201           CONTINUE 
              ENDIF
!             Check for unacounted atoms
              IF ( NNCNT.NE.0 ) THEN
                Ir = NNLIST(NNCNT)
! begin debug
!            WRITE(27,*)'   Check new n',NNCNT,Ir
! end debug                   
                NNCNT = NNCNT - 1
              ELSE
                FINDR = .FALSE.
                DO 4004 Ni = 1,APCNT
                  Ji = APLIST(Ni)
                  DO No = 1,ARCNT
                    Jo = ARLIST(No)
                    IF( Ji.EQ.Jo) GOTO 4004
                  ENDDO            
                  FINDR = .TRUE.
                  Ir = Ji
 4004           CONTINUE 
                DO 4005 No = 1,ARCNT
                  Jo = ARLIST(No)
                  DO Ni = 1,APCNT
                    Ji = APLIST(Ni)
                    IF( Ji.EQ.Jo) GOTO 4005
                  ENDDO            
                  FINDR  = .TRUE.
                  Ir = Jo
 4005           CONTINUE 
              ENDIF
            ENDDO            
!           All reactants and products found
    !       CALL PRSTR(ECNT,ELIST,REC,)
            IF ( MRCNT.LE.MXPRN )THEN
              WRITE( PRNR,*) RECTNT(1:MRCNT)
            ELSE
              WRITE( PRNR,*) '# of mol in reactents >MXPRN'
            ENDIF 
            IF ( MPCNT.LE.MXPRN )THEN
              WRITE( PRNP,*) PRODCT(1:MPCNT)
            ELSE
              WRITE( PRNP,*) '# of mol in products >MXPRN'
            ENDIF 
            WRITE(27,273) TTIME,LSTEP,REC,MRCNT,PRNR,MPCNT,PRNP
!           Store bond depth profile information
            ID = BID(KTI,KTJ)
            DP = R0i(DEPD,I)
            IF ( DP.GT.SUBMN.AND. DP.LT.SUBCT ) THEN
              DPTH =  SUBCT - DP
              DPTHDIV = DPTH/DPDIV + 1.0d0
              D = INT(DPTHDIV) 
              EXB(ID) = EXB(ID) + 1
              EXBDP(D,ID) = EXBDP(D,ID) + 1
            ELSE
              EXBGS(ID) = EXBGS(ID) + 1
            ENDIF
          ENDIF
!
c$$$          IF( .NOT.ASUBi(I).AND.(.NOT.ASUBo(I))) THEN
c$$$            NMi = IMOLi(Ir)
c$$$            NMo = IMOLo(Ir)
c$$$            Io=MPNTo(NMo)
c$$$            Jo=MPNTo(NMo+1)-1
c$$$            Ii=MPNTi(NMi)
c$$$            Ji=MPNTi(NMi+1)-1
c$$$            MMATCH = .FALSE.
c$$$            DO 1003 Oo=Io,Jo
c$$$                IMo = MOLSTo(Oo)
c$$$                DO Oi=Ii,Ji
c$$$                  IMi = MOLSTi(Oi)
c$$$                  IF( IMo.EQ.IMi) GOTO 1003
c$$$                ENDDO
c$$$                MMATCH = .TRUE.
c$$$! begin debug
c$$$          WRITE(27,*)' Molecule R and P the same',LSTEP
c$$$          WRITE(27,*)' Initial atom',I,KTYPEi(I)
c$$$          WRITE(27,*)' Neibors o',LISTo(NABORSo(I):NABORSo(I+1)-1)
c$$$          WRITE(27,*)' Neibors i',LISTi(NABORSi(I):NABORSi(I+1)-1)
c$$$          WRITE(27,*) 'Reaction type',RECTi(I)
c$$$          WRITE(27,*) KTYPEi(I),DELN,NNo
c$$$          WRITE(27,*) 
c$$$!          STOP
c$$$! end debug                   
c$$$ 1003       CONTINUE
c$$$             
c$$$            IF(MMATCH) GOTO 1001
c$$$          ENDIF
! begin debug
!          WRITE(27,*)'   Check new n',NNCNT
! end debug                   
          IF ( NNCNT.GT.NMX) THEN
            WRITE(*,*) 'Atom ',I,'has',NNCNT,' change in nieghbors'
     &  ,'which is more than the ',NMX,' alowed'
     &  ,'change NMX in module analysis'
            STOP
          ENDIF

! begin debug
c$$$          WRITE(27,*)' Initial atom',I,KTYPEi(I)
c$$$          WRITE(27,*)' Asub', ASUBi(I),ASUBo(I)
c$$$          WRITE(27,*)' Neibors o',LISTo(NABORSo(I):NABORSo(I+1)-1)
c$$$          WRITE(27,*)' Neibors i',LISTi(NABORSi(I):NABORSi(I+1)-1)
c$$$          WRITE(27,*)' R Molec cnt',MRCNT
c$$$          WRITE(27,*)' P Molec cnt',MPCNT
c$$$          WRITE(27,*)' R atomic cnt',ARCNT
c$$$          WRITE(27,*)  ARLIST(1:ARCNT)
c$$$          WRITE(27,*)' P atomic cnt',APCNT
c$$$          WRITE(27,*)  APLIST(1:APCNT)
c$$$          WRITE(27,*) 'Reaction type',RECTi(I)
c$$$          WRITE(27,*) KTYPEi(I),DELN,NNo
c$$$          WRITE(27,*) 
!          STOP
! end debug         
        ENDIF
 1001 CONTINUE
!     Print hybridization and other element properties
      WRITE(31,311) TTIME,LSTEP,SPCNT,SP2CNT,SP3CNT
     & ,ELCOMP(1,icarb),ELCOMP(2,icarb),GASKE(icarb)
     & ,ELCOMP(1,ihyd),ELCOMP(2,ihyd),GASKE(ihyd)
     & ,ELCOMP(1,isulfur),ELCOMP(2,isulfur),GASKE(isulfur)
     & ,ELCOMP(1,ioxygen),ELCOMP(2,ioxygen),GASKE(ioxygen)
     & ,ELCOMP(1,iflor),ELCOMP(2,iflor),GASKE(iflor)
!     Print bond depth profile 
        DO ID = 1,BTYPES
          WRITE(BRDS,*) BRBDP(:,ID)
          WRITE(FRDS,*) FRBDP(:,ID)
          WRITE(EXDS,*) EXBDP(:,ID)
          WRITE(33,331) ' Broken ',TTIME,LSTEP,ID,BRBGS(ID)
     &                  ,BRB(ID),BRDS
          WRITE(33,331) ' Formed ',TTIME,LSTEP,ID,FRBGS(ID)
     &                  ,FRB(ID),FRDS
          WRITE(33,331) ' Exchange ',TTIME,LSTEP,ID,EXBGS(ID)
     &                  ,EXB(ID),EXDS
        ENDDO
!
!     Save previous list as o
      NABORSo(:) = NABORSi(:)
      LISTo(:)   = LISTi(:)
!
! begin debug
!      WRITE(*,*) ' Reaction check finished current recs: ',REC
! end debug  
      ENDIF
      RETURN
 271  FORMAT( F16.1,2I10,' Broken ',2A4,2I8
     &       ,I8,' Reactants: ',A,I8,' Products: ',A )
 272  FORMAT( F16.1,2I10,' Formed ',2A4,2I8
     &         ,I8,' Reactants: ',A,I8,' Products: ',A)
 273  FORMAT( F16.1,2I10,' Exchange '
     &         ,I8,' Reactants: ',A,I8,' Products: ',A)
 311  FORMAT(F16.1,I12,3I10,5(2I9,F16.6))
 331  FORMAT(A,F16.2,I10,I4,2I8,A) 
!
      END SUBROUTINE react 
