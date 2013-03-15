C-----------------------------------------------------------------------
      SUBROUTINE MINRES (NG, G, D, X, VHAT, TOLRES, ITERMAX, NITER)
C-----------------------------------------------------------------------
C     MINRES with preconditioner P = L^{-1} such that L L^T \approx A.
C     Beware, D = diag(A) is given as input. The iterative process 
C     stopped when the residual tolerance or maximum number of 
C     iterations is reached. The initial guess is zero.
C
C     Taken from "Polynomial based iteration methods for symmetric 
C     linear systems", Bernd Fischer, Wiley & Teubner, New York, 
C     Leipzig, 1996.
C-----------------------------------------------------------------------
      IMPLICIT NONE 
      INTEGER NG, ITERMAX, NITER
      DOUBLE PRECISION G(1:(NG*NG+1)/2), D(NG), X(NG), VHAT(NG), 
     &                 TOLRES
C-----------------------------------------------------------------------
      INTEGER I
      DOUBLE PRECISION V(NG), VOLD(NG), W(NG), WOLD(NG), WOOLD(NG), 
     &                 TEMP(NG), TOL, BETA, BETAOLD, ETA, C, COLD, 
     &                 COOLD, S, SOLD, SOOLD, RES, ALPHA, DELTA, RHO1, 
     &                 RHO2, RHO3, DIF, ERROR, P(NG)
C-----------------------------------------------------------------------
      TOL = TOLRES *TOLRES

C     Initialization
      NITER = 0 
      BETA = 0.D0
      DO I = 1, NG
        P(I)    = 1.D0 /DSQRT(D(I))
        V(I)    = 0.D0
        VHAT(I) = P(I) *VHAT(I)
        BETA    = BETA +VHAT(I) *VHAT(I)
        W(I)    = 0.D0
        WOLD(I) = 0.D0
        X(I)    = 0.0D0
      ENDDO
      BETA = DSQRT(BETA) ; ETA = BETA
      C = 1.D0; COLD = 1.D0; S = 0.D0; SOLD = 0.D0

      RES = 1.D0; ERROR = 1.D0
      DO WHILE ( (NITER < ITERMAX) .AND. 
     &           ( (RES > TOL) .OR. (ERROR > TOL) ) )
        NITER = NITER + 1

C     The Lanczos recurrence
        DO I = 1, NG
          VOLD(I) = V(I) 
          V(I)    = VHAT(I) /BETA
          TEMP(I) = V(I) *P(I)
        ENDDO
        CALL SYMMATVEC(NG, G, TEMP, VHAT)
        ALPHA = 0.D0  
        DO I = 1, NG
          VHAT(I) = VHAT(I) *P(I)
          ALPHA   = ALPHA +VHAT(I) *V(I)
        ENDDO
        DO I = 1, NG
          VHAT(I) = VHAT(I) -ALPHA *V(I) - BETA *VOLD(I)
        ENDDO

        BETAOLD = BETA; BETA = 0.D0
        DO I    = 1, NG
          BETA  = BETA +VHAT(I) *VHAT(I)
        ENDDO
        BETA    = DSQRT(BETA) 

C     QR factorization:
C     A. Old given rot's on new colum of T
        COOLD = COLD; COLD = C; SOOLD = SOLD; SOLD = S

        DELTA = COLD *ALPHA -COOLD *SOLD *BETAOLD
        RHO1  = DSQRT(DELTA *DELTA +BETA *BETA)
        RHO2  = SOLD *ALPHA +COOLD *COLD *BETAOLD
        RHO3  = SOOLD *BETAOLD

C     B. New given rotation for subdiag elt
        C = DELTA /RHO1
        S = BETA  /RHO1

C     Update of solution (with W = V R^{-1})
        ERROR = 0.D0
        DO I = 1, NG
          WOOLD(I) = WOLD(I); WOLD(I) = W(I)
          W(I)     = (V(I) -RHO3 *WOOLD(I) -RHO2 *WOLD(I)) /RHO1
          DIF      = C *ETA *W(I) *P(I)
          ERROR    = ERROR +DIF *DIF
          X(I)     = X(I)  +DIF
        ENDDO

        ETA = -S *ETA
        RES = RES *S *S

      ENDDO

      END SUBROUTINE MINRES 
C-----------------------------------------------------------------------
