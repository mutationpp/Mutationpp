C-----------------------------------------------------------------------
      SUBROUTINE BERNOULLI (P, T, RHO, U, PTB)
C-----------------------------------------------------------------------
C     This subroutine computes the total pressure of a mixture using
C     Bernoulli's formula
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION P, T, RHO, U, PTB
C-----------------------------------------------------------------------
      PTB = P +0.5D0 *RHO *U *U

      END SUBROUTINE BERNOULLI
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE TOTAL (WR1, LWR1, WI, LWI, T, P, U, XN, X, EPS,  
     &                  TT, PT)
C-----------------------------------------------------------------------
C     This subroutine computes the total temperature and pressure of a
C     mixture in thermo-chemical equilibrium from the conservation of
C     entropy and total enthalpy. 
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      INTEGER LWR1, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), T, P, U, XN(1:NC), X(1:NS), 
     &                 EPS, TT, PT
C-----------------------------------------------------------------------
      DOUBLE PRECISION MM, R0(1:2), C(1:2), ND, XINI(1:NS), DELTA(1:2),
     &                 R(1:2), J(1:2,1:2), EPSP1, TV(1:NV), TTV(1:NV), GAMMAI, 
     &                 SOUNDE, GAMMAE, ALPHA
      EXTERNAL ENTHALPY      
C-----------------------------------------------------------------------
      EPSP1 = EPS +1.D0
      TV(1) = 1.D0  !? why do you need this?
      DO I = 1, NV
        TV(I)  = T
      ENDDO

C     Mixture molar mass
      MM = 0.D0
      DO I = 1, NS
        MM  = MM  +X(I)  *WR1(IMI+I-1)
        XINI(I) = X(I) 
      ENDDO

C     Mixture enthalpy at P and T (per unit mass)
      CALL MIXPROPERTY (WR1, LWR1, WI, LWI, T, T, T, TV, DUM, X,
     &                  H1, H2, H3, H4, H5, H6, ENTHALPY)
      H1 = H1 /MM             
      WRITE(*,*)' CP1 ', H1/T
      WRITE(*,*)'H:', H1, H2, H3, H4, H5, H6
C     Mixture entropy at P and T (per unit mass)
      CALL MIXENTROPY (WR1, LWR1, WI, LWI, T, T, T, TV, P, X,
     &                  S1, S2, S3, S4, S5, S6)
      S1 = S1 /MM

C     Newton iterative procedure
C     A. Initialization 
      C(1) = H1 +0.5D0 *U *U;  C(2) = S1
      WRITE(*,*) "From static"
      WRITE(*,*) C
      CALL NUMBERD (WR1, LWR1, P, T, T, X, ND)
      CALL DENSITY (WR1, LWR1, X, ND, RHO)

C     Bernoulli total pressure
C      PT = P +0.5D0 *RHO *U *U; TT = T

C     Saint Venant total pressure
      CALL SAINTVENANTPT (WR1, LWR1, WI, LWI, T, P, U, XN, X, EPS,
     &                    TT, PT, GAMMAE, GAMMAI, SOUNDE)
     
      DO I = 1,NV
        TTV(I) = TT
      END DO
     
      CALL MIXPROPERTY (WR1, LWR1, WI, LWI, TT, TT, TT, TV, DUM, X,
     &                  H1, H2, H3, H4, H5, H6, ENTHALPY)
      H1 = H1 /MM               

C     Mixture entropy at P and T (per unit mass)
      CALL MIXENTROPY (WR1, LWR1, WI, LWI, TT, TT, TT, TV, PT, X,
     &                  S1, S2, S3, S4, S5, S6)
      S1 = S1 /MM
      
      WRITE(*,*)' CP2 ', H1/TT
      
      WRITE(*,*)'H:', H1, H2, H3, H4, H5, H6
      C(1) = H1
      C(2) = S1
      WRITE(*,*) "From total"
      WRITE(*,*) C 

      CALL RESIDUAL (WR1, LWR1, WI, LWI, TT, PT, XN, XINI, C, R0)
      RES = R0(1) *R0(1) +R0(2) *R0(2)
      RESINI = RES
      RES = 1.D0
      
      IF (RESINI < RESMIN) RETURN

C     B. Test (residual and maximum number of iterations)
      RESMIN = 1.D-16; ITERMAX = 100; ITER = 0
      RESMIN = RESMIN*RESMIN
      ALPHA = 1.D0
      DO WHILE ( (ITER < ITERMAX) .AND.  (RES > RESMIN) )
        ITER =ITER +1
        WRITE(*,*) RES

C     C. Jacobian
        CALL RESIDUAL (WR1, LWR1, WI, LWI, TT, PT*EPSP1, XN, XINI,
     &                 C, R)
        J(1,1) = (R(1) -R0(1)) /(PT *EPS)
        J(2,1) = (R(2) -R0(2)) /(PT *EPS)
        CALL RESIDUAL (WR1, LWR1, WI, LWI, TT*EPSP1, PT, XN, XINI,
     &                 C, R)
        J(1,2) = (R(1) -R0(1)) /(TT *EPS)
        J(2,2) = (R(2) -R0(2)) /(TT *EPS)

C     D. Linear solution 
        DET = J(1,1) *J(2,2) -J(1,2) *J(2,1)
        IF (DET < 1.D-32) THEN
C          WRITE(*,*) 'Singular system'
          DELTA(1) = 0.D0
          DELTA(2) = 0.D0
        ELSE
          DELTA(1) = (-J(2,2) *R0(1) +J(1,2) *R0(2)) /DET
          DELTA(2) = (-J(1,1) *R0(2) +J(2,1) *R0(1)) /DET
        ENDIF

        PT = PT + ALPHA *DELTA(1)
        TT = TT + ALPHA *DELTA(2)

C     E. Residual
        CALL RESIDUAL (WR1, LWR1, WI, LWI, TT, PT, XN, XINI, C, R0)
        RES = R0(1) *R0(1) +R0(2) *R0(2)
        RES = RES /RESINI
      ENDDO
      WRITE(*,*) RES, ITER

      END SUBROUTINE TOTAL 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE RESIDUAL (WR1, LWR1, WI, LWI, TT, PT, XN, XINI, 
     &                     C, R)
C-----------------------------------------------------------------------
C     This subroutine computes the residual of the system for total 
C     temperature and pressure.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      INTEGER LWR1, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), TT, PT, XN(1:NC), XINI(1:NS), 
     &                 C(1:2), R(1:2)
C-----------------------------------------------------------------------
      DOUBLE PRECISION MM, X(1:NS), TTV(1:NV)
      EXTERNAL ENTHALPY 
C-----------------------------------------------------------------------
      TTV(1) = 1.D0
      DO J = 1, NV
        TTV(J)  = TT
      ENDDO

C     Composition at PT and TT
      CALL COMPOSITION (WR1, LWR1, WI, LWI, TT, PT, XN, XINI, X)
      MM = 0.D0
      DO I = 1, NS
        MM      = MM  +X(I)  *WR1(IMI+I-1)
        XINI(I) = X(I)
      ENDDO

C     Mixture enthalpy at PT and TT (per unit mass)
      CALL MIXPROPERTY (WR1, LWR1, WI, LWI, TT, TT, TT, TTV, DUM, X,
     &                  H1, H2, H3, H4, H5, H6, ENTHALPY)
      H1 = H1 /MM 

C     Mixture entropy at PT and TT (per unit mass)
      CALL MIXENTROPY (WR1, LWR1, WI, LWI, TT, TT, TT, TTV, PT, X,
     &                  S1, S2, S3, S4, S5, S6)
      S1 = S1 /MM

C     Residual
      R(1) = H1 -C(1)
      R(2) = S1 -C(2)

      END SUBROUTINE RESIDUAL
C---------------------------------------------------------------------- 
C-----------------------------------------------------------------------
      SUBROUTINE SAINTVENANT (P, MACH, GAMMA, PT)
C-----------------------------------------------------------------------
C     This subroutine computes the total pressure of a mixture using
C     Saint-Venant Wantzel's formula.
C     This formula is not exact for reacting flows.
C     Choice of gamma: the isentropic exponent give more accurate 
C     results than the ratio of specific heats.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION P, MACH, GAMMA, PT
C-----------------------------------------------------------------------
      PT = P *(1.0D0 +0.5D0*(GAMMA -1.D0) *MACH *MACH)
     &     **(GAMMA /(GAMMA -1.D0))

      END SUBROUTINE SAINTVENANT
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SAINTVENANTPT (WR1, LWR1, WI, LWI, T, P, U, XN, X, EPS,
     &                          TT, PT, GAMMAE, GAMMAI, SOUNDE)
C-----------------------------------------------------------------------
C     This subroutine computes the total pressure of a mixture using
C     Saint-Venant Wantzel's formula based on the isentropic exponent.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      INTEGER LWR1, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), T, P, U, XN(1:NC), X(1:NS), 
     &                 EPS, TT, PT, GAMMAE, GAMMAI, SOUNDE
C-----------------------------------------------------------------------
      DOUBLE PRECISION MM, EPSP1, TV(1:NV), DRHODP, MACH, ND, RHO 
      EXTERNAL ENTHALPY 
C-----------------------------------------------------------------------
      EPSP1 = EPS +1.D0
      TV(1) = 1.D0
      DO I = 1, NV
        TV(I)  = T
      ENDDO

C     Mixture molar mass
      MM = 0.D0
      DO I = 1, NS
        MM  = MM  +X(I)  *WR1(IMI+I-1)
      ENDDO

C     Mixture enthalpy at P and T (per unit mass)
      CALL MIXPROPERTY (WR1, LWR1, WI, LWI, T, T, T, TV, DUM, X,
     &                  H1, H2, H3, H4, H5, H6, ENTHALPY)
      H1 = H1 /MM               

      CALL NUMBERD (WR1, LWR1, P, T, T, X, ND)
      CALL DENSITY (WR1, LWR1, X, ND, RHO)


      CALL EQUIGAMMA (WR1, LWR1, WI, LWI, T, P, RHO, XN, X, EPS, 
     &                GAMMAE, DRHODP)
      SOUNDE = DSQRT(GAMMAE /DRHODP)
      GAMMAI = RHO /P *GAMMAE /DRHODP
      MACH   = U /SOUNDE
      
      WRITE(*,*) "Gammas: ", GAMMAE, GAMMAI
      WRITE(*,*) "sounde: ", SOUNDE
      WRITE(*,*) "MACH:   ", MACH
      
      CALL SAINTVENANT (P, MACH, GAMMAI, PT)
      CALL SAINTVENANTTEMP (T, MACH, GAMMAI, TT)

      END SUBROUTINE SAINTVENANTPT
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SAINTVENANTTEMP (T, MACH, GAMMA, TT)
C-----------------------------------------------------------------------
C     This subroutine computes the total temperature of a mixture using
C     Saint-Venant Wantzel's formula.
C     This formula is not exact for reacting flows.
C     Choice of gamma: the isentropic exponent give more accurate 
C     results than the ratio of specific heats.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION T, MACH, GAMMA, TT
C-----------------------------------------------------------------------
      TT = T *(1.0D0 +0.5D0*(GAMMA -1.D0) *MACH *MACH)

      END SUBROUTINE SAINTVENANTTEMP
C-----------------------------------------------------------------------
