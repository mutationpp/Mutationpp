C--------------------------- -------------------------------------------
      SUBROUTINE COMPOTOL (X, TOL, XTOL)
C-----------------------------------------------------------------------
C     This subroutine adds a small number to the composition respecting 
C     the mass constraint.
C     BEWARE! Don't used the modified mole fractions as initial guess to
C     compute the new composition, otherwise the Newton convergence will
C     be killed.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(1:NS), XTOL(1:NS), TOL
C-----------------------------------------------------------------------
      INTEGER I
      DOUBLE PRECISION SUM 
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      SUM = 0.D0
      DO I = 1, NS
        XTOL(I) = X(I) +TOL
        SUM = SUM +XTOL(I)
      ENDDO
      DO I = 1, NS
        XTOL(I) = XTOL(I) /SUM 
      ENDDO

C-----------------------------------------------------------------------
      END SUBROUTINE COMPOTOL 
C--------------------------- -------------------------------------------
C--------------------------- -------------------------------------------
      SUBROUTINE COMPOSITION (WR1, LWR1, WI, LWI, TT, P, XN, XINI, X)
C-----------------------------------------------------------------------
C     This subroutine computes the composition (molar fractions) of a 
C     mixture in thermo-chemical equilibrium for prescribed pressure
C     and temperature.
C     Due to the ill-conditioning of the system at low temperatures, 
C     the minimal temperature is set to TMIN. TMIN = 250 K is a sound 
C     value. 
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), P, TT, XN(1:NC), XINI(1:NS), 
     &                 X(1:NS)
C-----------------------------------------------------------------------
      INTEGER ITER, ITERMAX, INDX(1:NC+1)
      DOUBLE PRECISION LNK(1:NR), V(1:NS+1), RES2(1:NR), VAR, TV(1:NV),
     &                 A(1:NC+1,1:NC+1), B(1:NC+1), MU(1:NS), D
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
C     Minimal temperature
      IF (TT < WR1(ITMIN)) THEN
        RESMIN = 1.D-6
        T = WR1(ITMIN)
      ELSE
        RESMIN = 1.D-14
        T = TT
      ENDIF
      TV(1) = T
      DO I = 1, NV
        TV(I) = T
      ENDDO

C     Equilibrium constants
      CALL GIBBS (WR1, LWR1, WI, LWI, T, T, T, TV, 1.D0, MU)
      CALL EQUILIBRIUM (NS, NR, NC, P, T, MU, WI(INUIK), WR1(IUR), LNK)

C     Change of variables for initial composition
      DO I = 1, NS
        IF (XINI(I) <= 0.D0) THEN
          XINI(I) = 1.D0
        ENDIF
      ENDDO
      DO I = 1, NC
        V(I) = DLOG(XINI(I))
      ENDDO
      V(NC+1) = 1.D0
      DO I = 1, NR
        V(NC+1+I) = DLOG(XINI(NC+I))
      ENDDO

C     Newton iterative procedure
C     A. Initialization
      ITER = 0 
      CALL SYST (NS, NC, NR, WI(INUIK), LNK, XN, V, RES2, A, B)
      RES = 0.D0
      DO I = 1, NC+1
        RES = RES + B(I)*B(I)
      ENDDO

C     B. Test (residual and maximum number of iterations)
C      RESMIN = RESMIN *RESMIN; ITERMAX = 100
      RESMIN = RESMIN *RESMIN; ITERMAX = 40
      
      DO WHILE ( (ITER < ITERMAX) .AND.  (RES > RESMIN) )
       ITER =ITER +1

C     C. Schur complement variables
        CALL LUDCMP(A, NC+1, NC+1, INDX, D)
        CALL LUBKSB(A, NC+1, NC+1, INDX, B)
        DO I = 1, NC+1
          V(I) = V(I) +B(I)
        ENDDO       

C     D. Other variables
        DO I = 1, NR
          V(NC+1+I) = V(NC+1+I) -RES2(I) /WI(INUIK-1+(NC+I-1)*NS+NC+I)
          DO J = 1, NC+1
            V(NC+1+I) = V(NC+1+I) -B(J) *WR1(IJ22M1J21-1+(I-1)*(NC+1)+J)
          ENDDO          
        ENDDO

C     E. Residual
        CALL SYST (NS, NC, NR, WI(INUIK), LNK, XN, V, RES2, A, B)

        RES = 0.D0
        DO I = 1, NC+1
          RES = RES + B(I)*B(I)
        ENDDO
      ENDDO     

      IF (ITER == ITERMAX) THEN
        WRITE(*,*) 'Maximum number of iterations exceeded...'
        WRITE(*,*) '   RESMIN', RESMIN
        WRITE(*,*) '   RES', RES
      ENDIF

C     Back to molar fractions.
      DO I = 1, NC
        X(I) = DEXP(V(I))
      ENDDO
      DO I = 1, NR
        X(NC+I) = DEXP(V(NC+1+I))
      ENDDO

      END SUBROUTINE COMPOSITION 
C----------------------------------------------------------------------- 
C-----------------------------------------------------------------------
      SUBROUTINE SYST (NS, NC, NR, NU, LNK, XN, V, RES2, A, B) 
C-----------------------------------------------------------------------
C     This subroutine computes the system obtained after linearization 
C     of the nuclear conservation, chemical reaction equilibrium and 
C     mass constraint equations in thermal equilibrium (molar 
C     fractions).
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NS, NC, NR, NU(NS*NS)
      DOUBLE PRECISION LNK(1:NR), V(1:NS+1), XN(1:NC), A(1:NC+1,1:NC+1),
     &                 B(1:NC+1), RES2(1:NR)
C-----------------------------------------------------------------------
      INTEGER I, J, K
      DOUBLE PRECISION X(1:NS), J11(1:(NC+1),1:NC+1), TOL,
     &                 J12INVJ22(1:NC+1,NR), RES1(1:NC+1)
C-----------------------------------------------------------------------
C     Molar fractions
      TOL = DLOG(1.D-99)
      DO I = 1, NC
        IF (V(I) > 0.D0) THEN
          V(I) = 0.D0
        ENDIF
        IF (V(I) < TOL) THEN
          V(I) = TOL
        ENDIF
        X(I) = DEXP(V(I))
      ENDDO
      DO I = 1, NR 
        IF (V(NC+1+I) > 0.D0) THEN
          V(NC+1+I) = 0.D0
        ENDIF
        IF (V(NC+1+I) < -228.D0) THEN
          V(NC+1+I) = -228.D0
        ENDIF
        X(NC+I) = DEXP(V(NC+1+I))
      ENDDO
      
C     Residuals
      DO I = 1, NC
        RES1(I) = NU((I-1)*NS+I) *(X(I) -XN(I) *V(NC+1))     
        DO J = 1 ,NR
          RES1(I) = RES1(I) + NU((I-1)*NS+NC+J) *X(NC+J)
        ENDDO
      ENDDO
      RES1(NC+1) = -1.D0
      DO I = 1, NS
        RES1(NC+1) = RES1(NC+1) +X(I)
      ENDDO

      DO I = 1, NR
        RES2(I) = NU((NC+I-1)*NS+NC+I) *V(NC+I+1) -LNK(I)
        DO J = 1, NC
          RES2(I) =  RES2(I) +NU((NC+I-1)*NS+J) *V(J)
        ENDDO
      ENDDO

C     Jacobian matrix
      DO I = 1, NC+1
        DO J = 1, NC +1
          J11(I,J) = 0.D0
        ENDDO
      ENDDO
      DO I = 1, NC
        J11(I,I) = NU((I-1)*NS+I)*X(I)
        J = NC +1
          J11(I,J) = -XN(I)*NU((I-1)*NS+I)
      ENDDO
      I = NC +1
        DO J = 1, NC
          J11(I,J) = X(J)
        ENDDO

      DO I = 1, NC
        DO K = 1, NR
          J12INVJ22(I,K)= NU((I-1)*NS+NC+K) *X(NC+K)
     &                    /NU((NC+K-1)*NS+NC+K)
        ENDDO
      ENDDO
      I = NC+1
        DO K = 1, NR
          J12INVJ22(I,K)= X(NC+K) /NU((NC+K-1)*NS+NC+K)
        ENDDO

C     Schur complements
      DO I = 1, NC+1
        DO J = 1, NC
          A(I,J) = J11(I,J) 
          DO K = 1, NR
            A(I,J) = A(I,J) -J12INVJ22(I,K)*NU((NC+K-1)*NS+J)
          ENDDO
        ENDDO
        J = NC +1
          A(I,J) = J11(I,J)
      ENDDO

      DO I = 1, NC+1
        B(I) = -RES1(I) 
        DO K = 1, NR
          B(I) = B(I) + J12INVJ22(I,K) *RES2(K)
        ENDDO 
      ENDDO

      END SUBROUTINE SYST
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE EQUILIBRIUM (NS, NR, NC, P, T, MU, NU, UR, LNK)
C-----------------------------------------------------------------------
C     This subroutine computes the logarithm of the equilibrium constant
C     in thermal equilibrium (molar fractions).
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NS, NR, NC, NU(1:NS*NS)
      DOUBLE PRECISION P, T, UR, LNK(1:NR), MU(1:NS)
C-----------------------------------------------------------------------
      DOUBLE PRECISION MUSUM, NUSUM, LNP
      INTEGER I, J
C-----------------------------------------------------------------------
      LNP = DLOG(P)
      DO I = 1, NR
        J = I
        MUSUM = NU((NC+I-1)*NS+NC+J) *MU(NC+J)
        NUSUM = NU((NC+I-1)*NS+NC+J)
        DO J = 1, NC
          MUSUM = MUSUM + NU((NC+I-1)*NS+J) *MU(J)
          NUSUM = NUSUM + NU((NC+I-1)*NS+J)
        ENDDO
        LNK(I) = -MUSUM /(UR *T) -LNP *NUSUM
      ENDDO

      END SUBROUTINE EQUILIBRIUM
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE NUCLEAR (WR1, LWR1, XN) 
C-----------------------------------------------------------------------
C     This subroutine provides the default nuclear fraction read in the 
C     chemistry file, if this quantity is not provided by the user. 
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      INTEGER LWR1
      DOUBLE PRECISION  WR1(1:LWR1), XN(1:NC) 
C-----------------------------------------------------------------------
      INTEGER I 
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      DO I = 1, NC
        XN(I) = WR1(IXN+I-1) 
      ENDDO

      END SUBROUTINE NUCLEAR
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE MASSCOMPOSITION (WR1, LWR1, WI, LWI, TT, RHO, XN, 
     &                            XINI, X)
C-----------------------------------------------------------------------
C     This subroutine computes the composition (mass fractions) of a 
C     mixture in thermo-chemical equilibrium for prescribed density
C     and temperature.
C     Due to the ill-conditioning of the system at low temperatures, 
C     the minimal temperature is set to TMIN. TMIN = 250 K is a sound 
C     value. 
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), RHO, TT, XN(1:NC), XINI(1:NS), 
     &                 X(1:NS)
C-----------------------------------------------------------------------
      INTEGER ITER, ITERMAX, INDX(1:NC+1)
      DOUBLE PRECISION LNK(1:NR), V(1:NS+1), RES2(1:NR), VAR, TV(1:NV),
     &                 A(1:NC+1,1:NC+1), B(1:NC+1), MU(1:NS), D
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
C     Minimal temperature
      IF (TT < WR1(ITMIN)) THEN
        T = WR1(ITMIN)
        RESMIN = 1.D-4
      ELSE
        T = TT
        RESMIN = 1.D-10
      ENDIF
      TV(1) = T
      DO I = 1, NV
        TV(I) = T
      ENDDO

C     Equilibrium constants
      CALL GIBBS (WR1, LWR1, WI, LWI, T, T, T, TV, 1.D0, MU)
      CALL MASSEQUILIBRIUM (NS, NR, NC, RHO, T, MU, WI(INUIK), WR1(IUR),
     &                      WR1(IMI), LNK)

C     Change of variables for initial composition
      DO I = 1, NS
        IF (XINI(I) <= 0.D0) THEN
          XINI(I) = 1.D0
        ENDIF
      ENDDO
      DO I = 1, NC
        V(I) = DLOG(XINI(I))
      ENDDO
      V(NC+1) = 1.D0
      DO I = 1, NR
        V(NC+1+I) = DLOG(XINI(NC+I))
      ENDDO

C     Newton iterative procedure
C     A. Initialization
      ITER = 0 
      CALL MASSSYST (NS, NC, NR, WI(INUIK), LNK, WR1(IMI), XN, V, RES2, 
     &               A, B)
      RES = 0.D0
      DO I = 1, NC+1
        RES = RES + B(I)*B(I)
      ENDDO

C     B. Test (residual and maximum number of iterations)
      RESMIN = RESMIN *RESMIN; ITERMAX = 40 
      DO WHILE ( (ITER < ITERMAX) .AND.  (RES > RESMIN) )
       ITER =ITER +1

C     C. Schur complement variables
        CALL LUDCMP(A, NC+1, NC+1, INDX, D)
        CALL LUBKSB(A, NC+1, NC+1, INDX, B)
        DO I = 1, NC+1
          V(I) = V(I) +B(I)
        ENDDO       

C     D. Other variables
        DO I = 1, NR
          V(NC+1+I) = V(NC+1+I) -RES2(I) /WI(INUIK-1+(NC+I-1)*NS+NC+I)
          DO J = 1, NC+1
            V(NC+1+I) = V(NC+1+I) -B(J) *WR1(IJ22M1J21-1+(I-1)*(NC+1)+J)
          ENDDO          
        ENDDO

C     E. Residual
        CALL MASSSYST (NS, NC, NR, WI(INUIK), LNK, WR1(IMI), XN, V, 
     &                 RES2, A, B)
        RES = 0.D0
        DO I = 1, NC+1
          RES = RES + B(I)*B(I)
        ENDDO

      ENDDO     
 
C     Back to mass fractions. 
      DO I = 1, NC
        X(I) = DEXP(V(I)) 
      ENDDO
      DO I = 1, NR
        X(NC+I) = DEXP(V(NC+1+I))
      ENDDO
     
      END SUBROUTINE MASSCOMPOSITION 
C----------------------------------------------------------------------- 
C-----------------------------------------------------------------------
      SUBROUTINE MASSSYST (NS, NC, NR, NU, LNK, MI, XN, V, RES2, A, B) 
C-----------------------------------------------------------------------
C     This subroutine computes the system obtained after linearization 
C     of the nuclear conservation, chemical reaction equilibrium and 
C     mass constraint equations in thermal equilibrium (mass fractions).
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NS, NC, NR, NU(NS*NS)
      DOUBLE PRECISION LNK(1:NR), V(1:NS+1), MI(1:NS), XN(1:NC), 
     &                 A(1:NC+1,1:NC+1), B(1:NC+1), RES2(1:NR)
C-----------------------------------------------------------------------
      INTEGER I, J, K
      DOUBLE PRECISION X(1:NS), J11(1:(NC+1),1:NC+1), TOL,
     &                 J12INVJ22(1:NC+1,NR), RES1(1:NC+1)
C-----------------------------------------------------------------------
C     Mass fractions
      TOL = DLOG(1.D-99)
      DO I = 1, NC
        IF (V(I) > 0.D0) THEN
          V(I) = 0.D0
        ENDIF
        IF (V(I) < TOL) THEN
          V(I) = TOL
        ENDIF
        X(I) = DEXP(V(I))
      ENDDO
      DO I = 1, NR 
        IF (V(NC+1+I) > 0.D0) THEN
          V(NC+1+I) = 0.D0
        ENDIF
        IF (V(NC+1+I) < -228.D0) THEN
          V(NC+1+I) = -228.D0
        ENDIF
        X(NC+I) = DEXP(V(NC+1+I))
      ENDDO
      
C     Residuals
      DO I = 1, NC
        RES1(I) = NU((I-1)*NS+I) *(X(I) /MI(I) -XN(I) *V(NC+1))     
        DO J = 1 ,NR
          RES1(I) = RES1(I) + NU((I-1)*NS+NC+J) *X(NC+J) /MI(NC+J)
        ENDDO
      ENDDO
      RES1(NC+1) = -1.D0
      DO I = 1, NS
        RES1(NC+1) = RES1(NC+1) +X(I)
      ENDDO

      DO I = 1, NR
        RES2(I) = NU((NC+I-1)*NS+NC+I) *V(NC+I+1) -LNK(I)
        DO J = 1, NC
          RES2(I) =  RES2(I) +NU((NC+I-1)*NS+J) *V(J)
        ENDDO
      ENDDO

C     Jacobian matrix
      DO I = 1, NC+1
        DO J = 1, NC +1
          J11(I,J) = 0.D0
        ENDDO
      ENDDO
      DO I = 1, NC
        J11(I,I) = NU((I-1)*NS+I) *X(I) /MI(I)
        J = NC +1
          J11(I,J) = -XN(I) *NU((I-1)*NS+I)
      ENDDO
      I = NC +1
        DO J = 1, NC
          J11(I,J) = X(J)
        ENDDO

      DO I = 1, NC
        DO K = 1, NR
          J12INVJ22(I,K)= NU((I-1)*NS+NC+K) *X(NC+K)
     &                    /(NU((NC+K-1)*NS+NC+K) *MI(NC+K))
        ENDDO
      ENDDO
      I = NC+1
        DO K = 1, NR
          J12INVJ22(I,K)= X(NC+K) /NU((NC+K-1)*NS+NC+K)
        ENDDO

C     Schur complements
      DO I = 1, NC+1
        DO J = 1, NC
          A(I,J) = J11(I,J) 
          DO K = 1, NR
            A(I,J) = A(I,J) -J12INVJ22(I,K)*NU((NC+K-1)*NS+J)
          ENDDO
        ENDDO
        J = NC +1
          A(I,J) = J11(I,J)
      ENDDO

      DO I = 1, NC+1
        B(I) = -RES1(I) 
        DO K = 1, NR
          B(I) = B(I) + J12INVJ22(I,K) *RES2(K)
        ENDDO 
      ENDDO

      END SUBROUTINE MASSSYST
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE MASSEQUILIBRIUM (NS, NR, NC, RHO, T, MU, NU, UR, MI,
     &                            LNK)
C-----------------------------------------------------------------------
C     This subroutine computes the logarithm of the equilibrium constant
C     in thermal equilibrium (mass fractions).
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NS, NR, NC, NU(1:NS*NS)
      DOUBLE PRECISION RHO, T, UR, LNK(1:NR), MU(1:NS), MI(1:NS)
C-----------------------------------------------------------------------
      DOUBLE PRECISION MUSUM, NUSUM, MASSSUM, LNRHORT, LNMI(1:NS)
      INTEGER I, J
C-----------------------------------------------------------------------
      LNRHORT = DLOG(RHO *UR *T)
      DO I = 1, NS
        LNMI(I) = DLOG(MI(I))
      ENDDO
      DO I = 1, NR
        J = I
        MUSUM   = NU((NC+I-1)*NS+NC+J) *MU(NC+J)
        NUSUM   = NU((NC+I-1)*NS+NC+J)
        MASSSUM = NU((NC+I-1)*NS+NC+J) *LNMI(NC+J)
        DO J = 1, NC
          MUSUM   = MUSUM   + NU((NC+I-1)*NS+J) *MU(J)
          NUSUM   = NUSUM   + NU((NC+I-1)*NS+J)
          MASSSUM = MASSSUM +NU((NC+I-1)*NS+J) *LNMI(J)
        ENDDO
        LNK(I) = -MUSUM /(UR *T) -LNRHORT *NUSUM
     &           +MASSSUM
      ENDDO

      END SUBROUTINE MASSEQUILIBRIUM
C-----------------------------------------------------------------------
