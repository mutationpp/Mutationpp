C-----------------------------------------------------------------------
      SUBROUTINE READINPUT (RULE, PMIN, PMAX, DELTAP, TMIN, TMAX, 
     &                      DELTAT, TOL, PATH, MIXTURE, U, NMIXMAX, 
     &                      NMIX, MIX, NSPCMAX, NSPC, SPC, IHSONINE, 
     &                      IESONINE, IMOD)
C-----------------------------------------------------------------------
C     This subroutine reads the parameters for mutation stand-alone.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER RULE, NMIXMAX, NMIX, MIX(1:NMIXMAX), NSPCMAX, NSPC, 
     &        SPC(1:NSPCMAX), IHSONINE, IESONINE, IMOD
      DOUBLE PRECISION PMIN, PMAX, DELTAP, TMIN, TMAX, DELTAT, TOL, U
      CHARACTER(10) MIXTURE
      CHARACTER(100) PATH 
C-----------------------------------------------------------------------
      INTEGER I, IN, LPATH, LCHAR
      PARAMETER (IN = 30)
      CHARACTER(4) COM1
      CHARACTER(80) FULLCOM1
C-----------------------------------------------------------------------
      IMOD = 0
      LPATH = LCHAR(PATH)
      U = 0.D0
      OPEN(UNIT=IN,FILE=PATH(1:LPATH)//'/input/parameters',STATUS='OLD')
      COM1 = '   '
      DO WHILE (COM1(1:4)/= 'STOP')
        READ(IN,*) FULLCOM1
        COM1 = FULLCOM1(1:4)
        IF (COM1(1:4) == 'Name') THEN
          READ(IN,*) MIXTURE
        ELSEIF (COM1(1:4) == 'Pres') THEN
          READ(IN,*) PMIN 
          READ(IN,*) PMAX
          READ(IN,*) DELTAP 
        ELSEIF (COM1(1:4) == 'Temp') THEN
          READ(IN,*) TMIN 
          READ(IN,*) TMAX
          READ(IN,*) DELTAT 
        ELSEIF (COM1(1:4) == 'Tole') THEN
          READ(IN,*) TOL
        ELSEIF (COM1(1:4) == 'Meth') THEN
          READ(IN,*) RULE 
        ELSEIF (COM1(1:4) == 'Velo') THEN
          READ(IN,*) U
        ELSEIF (COM1(1:4) == 'Heav') THEN
          READ(IN,*) IHSONINE
        ELSEIF (COM1(1:4) == 'Elec') THEN
          READ(IN,*) IESONINE
        ELSEIF (COM1(1:4) == 'Mixt') THEN
          READ(IN,*) NMIX
          READ(IN,*)  
          DO I = 1, NMIX
            READ(IN,*) MIX(I)
            IF (MIX(I) >= 100)  THEN
              IMOD = 1
            ENDIF
          ENDDO
        ELSEIF (COM1(1:4) == 'Spec') THEN
          READ(IN,*) NSPC
          READ(IN,*)  
          DO I = 1, NSPC
            READ(IN,*) SPC(I)
            IF (SPC(I) >= 100)  THEN
              IMOD = 1
            ENDIF
          ENDDO
        ENDIF 
      ENDDO
      CLOSE(IN)

      END SUBROUTINE READINPUT
C-----------------------------------------------------------------------

