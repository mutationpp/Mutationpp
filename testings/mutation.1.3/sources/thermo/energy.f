C-----------------------------------------------------------------------
      SUBROUTINE ENERGY (WR1, LWR1, WI, LWI, TH, TE, TR, TV, P, E1, 
     &                     E2, E3, E4, E5, E6)
C-----------------------------------------------------------------------
C     This subroutine computes the species energy per unit
C     mole in thermal non-equilibrium.
C     -E1(I): total species energy (including the formation enthalpy)
C     -E2(I): translational species energy
C     -E3(I): electronic species energy
C     -E4(I): rotational species energy
C     -E5(I): vibrational species energy
C     -E6(I): formation species energy
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      INTEGER LWR1, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), TH, TE, TR, TV(1:NV), P, E1(1:NS), 
     &                 E2(1:NS), E3(1:NS), E4(1:NS), E5(1:NS), E6(1:NS)
C-----------------------------------------------------------------------
      IC1 = 0; IC2 = 0
      DO I = 1, NS
        IATOMICITY = WI(IATOMI+I-1)
        SELECT CASE(IATOMICITY)
C         Electron
          CASE(0)
            E2(I) = 1.5D0 *WR1(IUR) *TE
            E3(I) = 0.D0
            E4(I) = 0.D0
            E5(I) = 0.D0
C         Atom
          CASE(1)
            E2(I) = 1.5D0 *WR1(IUR) *TH
            SUM1 = 0.D0; SUM2 = 0.D0; NELE = WI(IELEI+I-1)
            IF (NELE/=0) THEN
              DO J = 1, NELE
                EXPO = DEXP(-WR1(IEEIK+IC1+J-1) /TE)
                SUM1 = SUM1 +WI(IEGIK+IC1+J-1) *WR1(IEEIK+IC1+J-1) *EXPO
                SUM2 = SUM2 +WI(IEGIK+IC1+J-1) *EXPO
              ENDDO
            ELSE
              SUM2 = 1.D0
            ENDIF
            IC1 = IC1 +NELE
            E3(I) = WR1(IUR) *SUM1 /SUM2
            E4(I) = 0.D0
            E5(I) = 0.D0

C         Polyatomic molecule
          CASE DEFAULT
            E2(I) = 1.5D0 *WR1(IUR) *TH
            SUM1 = 0.D0; SUM2 =0.D0; NELE = WI(IELEI+I-1)
            IF (NELE/=0) THEN
              DO J = 1, NELE
                EXPO = DEXP(-WR1(IEEIK+IC1+J-1) /TE)
                SUM1 = SUM1 +WI(IEGIK+IC1+J-1) *WR1(IEEIK+IC1+J-1) *EXPO
                SUM2 = SUM2 +WI(IEGIK+IC1+J-1) *EXPO
              ENDDO
            ELSE
              SUM2 = 1.D0
            ENDIF
            IC1 = IC1 +NELE
            E3(I) = WR1(IUR) *SUM1 /SUM2
            IF (WI(ILINI+I-1)==1) THEN
              E4(I) = WR1(IUR) *TR
            ELSE
              E4(I) = 1.5D0 *WR1(IUR) *TR
            ENDIF
            SUM1 = 0.D0; NVIBMODE = WI(IVIBI+I-1)
            DO J = 1, NVIBMODE
              SUM1 = SUM1 + WR1(ITVIK+IC2+J-1) 
     &               /(DEXP(WR1(ITVIK+IC2+J-1) /TV(IC2+J)) -1.D0)
            ENDDO          
            IC2 = IC2 + NVIBMODE
            E5(I) = SUM1 *WR1(IUR)
        END SELECT
        E6(I) = WR1(IHFORI+I-1)
        E1(I) = E2(I) +E3(I) +E4(I) +E5(I) +E6(I)
      ENDDO
 
      END SUBROUTINE ENERGY 
C----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE ENERGYMASS (WR1, LWR1, WI, LWI, TH, TE, TR, TV, P, E1, 
     &                       E2, E3, E4, E5, E6)
C-----------------------------------------------------------------------
C     This subroutine computes the species energy per unit
C     mole in thermal non-equilibrium.
C     -E1(I): total species energy (including the formation enthalpy)
C     -E2(I): translational species energy
C     -E3(I): electronic species energy
C     -E4(I): rotational species energy
C     -E5(I): vibrational species energy
C     -E6(I): formation species energy
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      INTEGER LWR1, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), TH, TE, TR, TV(1:NV), P, E1(1:NS), 
     &                 E2(1:NS), E3(1:NS), E4(1:NS), E5(1:NS), E6(1:NS)
      DOUBLE PRECISION MASS 
C-----------------------------------------------------------------------
      IC1 = 0; IC2 = 0
      DO I = 1, NS
        MASS = WR1(IMI+I-1)
        IATOMICITY = WI(IATOMI+I-1)
        SELECT CASE(IATOMICITY)
C         Electron
          CASE(0)
            E2(I) = 1.5D0 *WR1(IUR) *TE /MASS
            E3(I) = 0.D0
            E4(I) = 0.D0
            E5(I) = 0.D0
C         Atom
          CASE(1)
            E2(I) = 1.5D0 *WR1(IUR) *TH /MASS
            SUM1 = 0.D0; SUM2 = 0.D0; NELE = WI(IELEI+I-1)
            IF (NELE/=0) THEN
              DO J = 1, NELE
                EXPO = DEXP(-WR1(IEEIK+IC1+J-1) /TE)
                SUM1 = SUM1 +WI(IEGIK+IC1+J-1) *WR1(IEEIK+IC1+J-1) *EXPO
                SUM2 = SUM2 +WI(IEGIK+IC1+J-1) *EXPO
              ENDDO
            ELSE
              SUM2 = 1.D0
            ENDIF
            IC1 = IC1 +NELE
            E3(I) = WR1(IUR) /MASS *SUM1 /SUM2
            E4(I) = 0.D0
            E5(I) = 0.D0

C         Polyatomic molecule
          CASE DEFAULT
            E2(I) = 1.5D0 *WR1(IUR) *TH /MASS
            SUM1 = 0.D0; SUM2 =0.D0; NELE = WI(IELEI+I-1)
            IF (NELE/=0) THEN
              DO J = 1, NELE
                EXPO = DEXP(-WR1(IEEIK+IC1+J-1) /TE)
                SUM1 = SUM1 +WI(IEGIK+IC1+J-1) *WR1(IEEIK+IC1+J-1) *EXPO
                SUM2 = SUM2 +WI(IEGIK+IC1+J-1) *EXPO
              ENDDO
            ELSE
              SUM2 = 1.D0
            ENDIF
            IC1 = IC1 +NELE
            E3(I) = WR1(IUR) /MASS *SUM1 /SUM2
            IF (WI(ILINI+I-1)==1) THEN
              E4(I) = WR1(IUR) *TR /MASS
            ELSE
              E4(I) = 1.5D0 *WR1(IUR) *TR /MASS
            ENDIF
            SUM1 = 0.D0; NVIBMODE = WI(IVIBI+I-1)
            DO J = 1, NVIBMODE
              SUM1 = SUM1 + WR1(ITVIK+IC2+J-1) 
     &               /(DEXP(WR1(ITVIK+IC2+J-1) /TV(IC2+J)) -1.D0)
            ENDDO          
            IC2 = IC2 + NVIBMODE
            E5(I) = SUM1 *WR1(IUR) /MASS
        END SELECT
        E6(I) = WR1(IHFORI+I-1) /MASS
        E1(I) = E2(I) +E3(I) +E4(I) +E5(I) +E6(I)
      ENDDO
 
      END SUBROUTINE ENERGYMASS 
C----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE TGUESS (WR1, LWR1, WI, LWI, RHOE, RHOI, TINI) 
C-----------------------------------------------------------------------
C     This subroutine computes an inital temperature guess for the 
C     Newton-Raphson iterative procedure used to compute T as a function 
C     of RHOE and RHOI. 
C-----------------------------------------------------------------------
      IMPLICIT NONE 
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      INTEGER LWR1, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), RHOE, RHOI(1:NS), TINI
      INTEGER I, IATOMICITY 
      DOUBLE PRECISION MASS, CVTR, SUM1, SUM2 
C-----------------------------------------------------------------------
      SUM1 = 0.D0; SUM2 = 0.D0
      DO I = 1, NS
        MASS = WR1(IMI+I-1)
        IATOMICITY = WI(IATOMI+I-1)
        SELECT CASE(IATOMICITY)
C         Electron
          CASE(0)
            CVTR = 1.5D0 *WR1(IUR) /MASS
C         Atom
          CASE(1)
            CVTR = 1.5D0 *WR1(IUR) /MASS
C         Polyatomic molecule
          CASE DEFAULT
            CVTR = 2.5D0 *WR1(IUR) /MASS
        END SELECT
        SUM1 = SUM1 + RHOI(I) * WR1(IHFORI+I-1) /MASS
        SUM2 = SUM2 + RHOI(I) * CVTR
      ENDDO
      TINI = (RHOE -SUM1) /SUM2
 
      END SUBROUTINE TGUESS 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE TCNEQNEWTON (WR1, LWR1, WI, LWI, RHOE, RHOI, TINI, T) 
C-----------------------------------------------------------------------
C     This subroutine computes the temperature as a function of RHOE and 
C     and RHOI based on a Newton-Raphson iterative procedure with an 
C     initial temperature TINI.
C-----------------------------------------------------------------------
      IMPLICIT NONE 
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      INTEGER LWR1, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), RHOE, RHOI(1:NS), TINI, T 
      INTEGER N, NMAX, IS, IV
      DOUBLE PRECISION INVINC, FTN, FTNP1, TV(1:NV), TP, TVP(1:NV)
      DOUBLE PRECISION EPS, EPSP1, TOL, SUM1, P, EM1(1:NS), EM2(1:NS),
     &                 EM3(1:NS), EM4(1:NS), EM5(1:NS), EM6(1:NS)
C-----------------------------------------------------------------------
C     EPS for finite difference computation
      EPS = 1.D-2; EPSP1 = EPS +1.D0
C     TOL for Newton method
      TOL = 1.D-8

C     Maximum number of iterations
      NMAX = 10
C     Dummy variable!
      P = 0.D0
      
      T = TINI
      N = 1
      DO WHILE (N < NMAX)
        DO IV = 1, NV
          TV(IV) = T
        ENDDO
        CALL ENERGYMASS (WR1, LWR1, WI, LWI, T, T, T, TV, P, EM1, EM2,
     &                   EM3, EM4, EM5, EM6)
        SUM1 = 0.D0
        DO IS = 1, NS
          SUM1 = SUM1 + RHOI(IS) *EM1(IS)
        ENDDO
        FTN = RHOE - SUM1

        TP = T *EPSP1
        DO IV = 1, NV
          TVP(IV) = TP
        ENDDO
        CALL ENERGYMASS (WR1, LWR1, WI, LWI, TP, TP, TP, TVP, P,
     &                   EM1, EM2, EM3, EM4, EM5, EM6)
        SUM1 = 0.D0
        DO IS = 1, NS
          SUM1 = SUM1 + RHOI(IS) *EM1(IS)
        ENDDO
        FTNP1 = RHOE - SUM1

        IF (ABS(FTN) < TOL*RHOE) THEN
          NMAX = 1
        ELSE
          INVINC = (FTN - FTNP1) / (EPS *FTN) 
          T = T *(1.D0 +1.D0 /INVINC)
        ENDIF
        N = N +1
      ENDDO

      END SUBROUTINE TCNEQNEWTON
C----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE TCEQNEWTON (WR1, LWR1, WI, LWI, RHOE, RHO, XN, YINI,
     &                       TINI, T, Y)
C-----------------------------------------------------------------------
C     This subroutine computes the temperature as a function of RHOE and 
C     and RHO in chemical equilibrium based on a Newton-Raphson 
C     iterative procedure with an initial temperature TINI.
C-----------------------------------------------------------------------
      IMPLICIT NONE 
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      INTEGER LWR1, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), RHOE, RHO, XN(1:NC), YINI(1:NS),  
     &                 TINI, T 
      INTEGER N, NMAX, IS, IV
      DOUBLE PRECISION INVINC, FTN, FTNP1, TV(1:NV), TP, TVP(1:NV)
      DOUBLE PRECISION EPS, EPSP1, TOL, SUM1, P, EM1(1:NS), EM2(1:NS),
     &                 EM3(1:NS), EM4(1:NS), EM5(1:NS), EM6(1:NS), 
     &                 Y(1:NS), YP(1:NS)
C-----------------------------------------------------------------------
C     EPS for finite difference computation
      EPS = 1.D-2; EPSP1 = EPS +1.D0
C     TOL for Newton method
      TOL = 1.D-8
C     Maximum number of iterations
      NMAX = 10
C     Dummy variable!
      P = 0.D0
      
      T = TINI
      DO IS = 1, NS
        Y(IS) = YINI(IS) 
      ENDDO
      N = 1
      DO WHILE (N < NMAX)
        DO IV = 1, NV
          TV(IV) = T
        ENDDO
        CALL ENERGYMASS (WR1, LWR1, WI, LWI, T, T, T, TV, P, EM1, EM2,
     &                   EM3, EM4, EM5, EM6)
        DO IS = 1, NS
          YINI(IS) = Y(IS) 
        ENDDO
        CALL MASSCOMPOSITION (WR1, LWR1, WI, LWI, T, RHO, XN, YINI, Y)
        SUM1 = 0.D0
        DO IS = 1, NS
          SUM1 = SUM1 + Y(IS) *EM1(IS)
        ENDDO
        FTN = RHOE - RHO *SUM1

        TP = T *EPSP1
        DO IV = 1, NV
          TVP(IV) = TP
        ENDDO
        CALL ENERGYMASS (WR1, LWR1, WI, LWI, TP, TP, TP, TVP, P,
     &                   EM1, EM2, EM3, EM4, EM5, EM6)
        DO IS = 1, NS
          YINI(IS) = Y(IS) 
        ENDDO
        CALL MASSCOMPOSITION (WR1, LWR1, WI, LWI, TP, RHO, XN, YINI, YP)
        SUM1 = 0.D0
        DO IS = 1, NS
          SUM1 = SUM1 + YP(IS) *EM1(IS)
        ENDDO
        FTNP1 = RHOE - RHO *SUM1

        IF (ABS(FTN) < TOL*RHOE) THEN
          NMAX = 1
        ELSE
          INVINC = (FTN - FTNP1) / (EPS *FTN) 
          T = T *(1.D0 +1.D0 /INVINC)
        ENDIF
        N = N +1
      ENDDO

      END SUBROUTINE TCEQNEWTON
C----------------------------------------------------------------------
