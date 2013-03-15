C-----------------------------------------------------------------------
      SUBROUTINE EGSCG1 ( NG, G, DMI, AN, ZN, RN, TEMP, ITERMX )
C-----------------------------------------------------------------------
C     Conjugate gradient with preconditioning (EGSLIB)
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION G(*), DMI(NG),  AN(NG), ZN(NG),  RN(NG), TEMP(NG)
C-----------------------------------------------------------------------
      NITER = 0
      DO I = 1, NG
         III = NG*(I-1) - (I*(I-1))/2 + I
         AN(I) = 0.0D0
         ZN(I) = 0.0D0
         DMI(I) = 1.0D0 / G(III)
      ENDDO
      BETAN = 0.0D0
      AAA = 0.0D0
      DO I = 1, NG
         AAA = AAA + DMI(I) * RN(I)*RN(I)
      ENDDO

 100  CONTINUE
      NITER = NITER + 1
      DO I = 1, NG
         ZN(I) = DMI(I)*RN(I) + BETAN*ZN(I)
      ENDDO
      CALL EGSAXS(NG, G, ZN, TEMP)
      BBB = DDDOT (NG, ZN, 1, TEMP, 1)
      DO I = 1, NG
         AN(I) = AN(I) + AAA/BBB*ZN(I)
         RN(I) = RN(I) - AAA/BBB*TEMP(I)
      ENDDO
      CCC = 0.0D0
      DO I = 1, NG
         CCC = CCC + DMI(I) * RN(I)*RN(I)
      ENDDO
      BETAN = CCC/AAA
      AAA   = CCC
      IF ( NITER .LT. ITERMX ) GO TO 100

      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE CG ( NG, G, D, AN, RN, TOLRES, ITERMAX, NITER)
C-----------------------------------------------------------------------
C     Conjugate gradient with preconditioning.
C     Iterative process stopped when residual/solution tolerance or 
C     maximum number of iterations is reached.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION G(*), D(NG), DMI(NG), AN(NG), ZN(NG), RN(NG), TEMP(NG),
     &          ANOLD(NG)
      LOGICAL LOOPIN
C-----------------------------------------------------------------------
      TOL = TOLRES *TOLRES
      NITER = 0 
      DO I = 1, NG
         AN(I)  = 0.0D0
         ZN(I)  = 0.0D0
         DMI(I) = 1.D0 /D(I)
      ENDDO
      BETAN = 0.0D0
      AAA = 0.0D0
      DO I = 1, NG
         AAA = AAA + DMI(I) * RN(I)*RN(I)
      ENDDO
      RESINI = AAA
      IF (RESINI < 1.D-32) THEN
        WRITE(*,*) 'RESINI < 1.D-32 in CG...'
        LOOPIN = .FALSE.
      ELSE
        LOOPIN = .TRUE.
      ENDIF
      RES    = 1.D0
      ERROR  = 1.D0

      DO WHILE ( (NITER < ITERMAX) .AND. (LOOPIN) .AND.
C     &           ( (RES > TOL) .OR. (ERROR > TOL) ) )
     &           (RES > TOL) )
        NITER = NITER + 1
        DO I = 1, NG
          ANOLD(I) = AN(I)
        ENDDO
        DO I = 1, NG
           ZN(I) = DMI(I)*RN(I) + BETAN*ZN(I)
        ENDDO
        CALL EGSAXS(NG, G, ZN, TEMP)
        BBB = DDDOT (NG, ZN, 1, TEMP, 1)
        IF (BBB <= 0.D0 ) THEN
          IF (ABS(BBB) > 1.D-32 ) THEN
            WRITE(*,*) BBB,' <= 0.D0 in CG...'
          ENDIF
          NITER = NITER -1
          LOOPIN = .FALSE.
        ELSE
          DO I = 1, NG
            AN(I) = AN(I) + AAA/BBB*ZN(I)
            RN(I) = RN(I) - AAA/BBB*TEMP(I)
          ENDDO
        ENDIF
        CCC = 0.0D0
        DO I = 1, NG
           CCC = CCC + DMI(I) * RN(I)*RN(I)
        ENDDO
        RES = CCC /RESINI
        IF (RES > TOL) THEN
          BETAN = CCC/AAA
        ENDIF
        AAA   = CCC
        ERROR = 0.D0
        DO I = 1, NG
          DIF = AN(I)-ANOLD(I)
          ERROR = ERROR + DIF *DIF
        ENDDO
      ENDDO

      END SUBROUTINE CG
C-----------------------------------------------------------------------
