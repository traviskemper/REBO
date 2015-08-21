C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C Dynamic Simulation                                                          C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C

      SUBROUTINE MODEL

      USE SPECIF
      USE MDSTP
      USE STRUCTR
      USE PARAMS
      USE POTS

      IMPLICIT none 

      INTEGER :: I,J,N,M,L,K
      REAL*8 :: RNPSQ,PRNRNP

C Main MD loop
      DO LSTEP=1,KVC
        IF ( INTG.EQ.1 ) THEN
          call cpred 
        ELSEIF( INTG.EQ.2) THEN
          CALL VVa
        ELSE
          WRITE(*,*) 'Integrator not specified'
          STOP
        ENDIF
C
c calculate energy and forces 
c set forces to zero
C
        DO I=1,NP
           eatom(i) = 0.0d0
           DO J=1,3
                rnp(i,j) = 0.0d0
           ENDDO
        ENDDO
C
C set potential energy to zero
        tote = 0.0d0
C
C       LCHK =1
        CALL nborlist
C
        CALL CAGUTS

        IF ( SYTP.EQ.1) THEN
            CALL PIBONDCH
        ELSEIF ( SYTP.EQ.2) THEN
            CALL PIBONDCHO
        ELSEIF ( SYTP.EQ.3) THEN
            CALL PIBONDCHS
        ELSEIF ( SYTP.EQ.4) THEN
            CALL PIBONDCHF
        ELSE
          WRITE(*,*) 'Potential type SYTP=',SYTP,' is undefined'
        ENDIF  
C
        if(ILJ) call ljguts
C
C
        call thermos 
C
        IF ( INTG.EQ.1 ) THEN
          call ccorr
        ELSEIF( INTG.EQ.2) THEN
          CALL VVb
        ELSE
          WRITE(*,*) 'Integrator not specified'
          STOP
        ENDIF
C
        IF ( mod(LSTEP,NSTR).eq.0 ) THEN
            IF(PXMOL) call xmol
            IF(PCFG) call CFG
        ENDIF
        IF(KFLAG.EQ.5) CALL BERE
        if(mod(LSTEP,maxkb).eq.0) then
           call write_data
           CALL WRITE_COORD 
         endif 
      ENDDO 
C
 2304 FORMAT(2F14.6)
 2501 FORMAT(2I6,3F14.6)
 2601 FORMAT(3F14.6,5I6,3F14.6,3I6,5F14.6)
C
      return
      end
C
