! Create cfg files for atomeye
!
      SUBROUTINE setcfg
!
      USE STRUCTR
      USE SPECIF
!
      IMPLICIT none
!
      INTEGER :: N  
      REAL*8  :: DM
!
      N = KVC/NSTR
!
      IF ( N.GT.999999) THEN
        WRITE(6,*) 'Warning more than 999999 cfg files are needed'
        WRITE(6,*) ' change filename format in cfg.f'
        STOP
      ENDIF
!     Set box size 
      DO N = 1,3 
        CBOX(N) = CUBE(N)
        DM = RMN(2,N) - RMN(1,N) 
        IF(CBOX(N).GT.DM+50.d0 ) THEN
          cbox(N) = DM*2.0d0
        ENDIF
        IF ( CBOX(N).LT.10.d0 ) THEN
          CBOX(N) = 10.d0
        ENDIF
      ENDDO
!     Create Cfgs folder
      CALL SYSTEM('sh -c "if [ ! -d Cfgs ]; then mkdir Cfgs ; fi"')
!      
      RETURN
      END SUBROUTINE setcfg 
! 
! Create cfg files for atomeye 
! 
      SUBROUTINE cfg
!
      USE STRUCTR 
      USE POTS 
      USE MDSTP
      USE SPECIF
      USE PARAMS 
!
      IMPLICIT none 
!
      INTEGER :: N,I,K,LP,D,J
      REAL*8, DIMENSION(3) :: R0F,R1V,R2A
      REAL*8 :: RNPSQ,PRNRNP,VSQ,ASQ,MFRAC,MBF,XX,TEMPKi
      CHARACTER *79 FILENAME
!
      N = LSTEP / NSTR + 1
      LP = 19
!     Write output
      WRITE (FILENAME, '("Cfgs/", I6.6, ".cfg")'),N
      OPEN(UNIT=LP,FILE=FILENAME,STATUS='UNKNOWN')
      write(LP,10)'Number of particles = ',NP
      write(LP,11)'A = 1.0'
      write(LP,12)'H0(1,1) =', cbox(1)
      write(LP,12)'H0(1,2) =', 0.
      write(LP,12)'H0(1,3) =', 0.
      write(LP,12)'H0(2,1) =', 0.
      write(LP,12)'H0(2,2) =', cbox(2)
      write(LP,12)'H0(2,3) =', 0.
      write(LP,12)'H0(3,1) =', 0.
      write(LP,12)'H0(3,2) =', 0.
      write(LP,12)'H0(3,3) =', cbox(3)
      write(LP,11)'.NO_VELOCITY.'
      write(LP,11)'entry_count = 13'
      write(LP,11)'auxiliary[0]= Atom #'
      write(LP,11)'auxiliary[1]= Thermostat'
      write(LP,11)'auxiliary[2]= Temperature [K]'
      write(LP,11)'auxiliary[3]= Force [A/fs]'
      write(LP,11)'auxiliary[4]= Ave_velocity [A/fs]'
      write(LP,11)'auxiliary[5]= Ave_acceleration [A/fs]'
!     Loop over all atom types 
      DO K = 1,NTYPES
        IF( NOA(K).GT.0 ) THEN
          WRITE(LP,*) XMASS(K)
          WRITE(LP,*) ATYPE(K)
          DO I=1,NP
            IF ( KTYPE(I).EQ.K) THEN
!             Calculate fractional coordinates  
              VSQ = 0.0d0
              ASQ = 0.0d0         
              DO D =1,3
                R0F(D) = R0(I,D)/CBOX(D) + .50d0
                R1V(D) = R1(I,D)/DELTA
                VSQ = VSQ + R1V(D)**2
                R2A(D) = R2(I,D)/DELTSQ
                ASQ = ASQ + R2A(D)**2
              ENDDO 
!             Caclulate force 
              RNPSQ = 0.0d0
              DO  J=1,3
                RNPSQ = RNPSQ + RNP(I,J)**2
              ENDDO
              PRNRNP = SQRT(RNPSQ)
              XX = VSQ*XMASS(KTYPE(I))
              TEMPKi = ENPR*XX/6.0d0/DELTSQ*ECONV
              WRITE(LP,13)R0F(1:3),I,ITR(I),TEMPKi
     &          ,PRNRNP,SQRT(VSQ),SQRT(ASQ)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      CLOSE(LP)
!
10    format(a,1x,i6)
11    format(a)
12    format(a,1x,f26.6)
13    format(3(E12.4,1x),I8,I5,E12.4,1x,3(E12.4,1x))
14    format(3E20.11)
      RETURN
      END SUBROUTINE cfg


