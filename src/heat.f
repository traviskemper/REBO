      SUBROUTINE HEAT

      USE SPECIF
      USE STRUCTR
      USE PARAMS
      USE POTS
      USE MDSTP

      IMPLICIT none 

      INTEGER :: I

      IF ( TEM.LT.TMAX) THEN
        TEM = TEM + DELTMP
      ENDIF

      CALL setgle ! Recalc Thermostate param
      
      RETURN
      END SUBROUTINE HEAT


