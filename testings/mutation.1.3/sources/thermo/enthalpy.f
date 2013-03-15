C-----------------------------------------------------------------------
      SUBROUTINE ENTHALPY (WR1, LWR1, WI, LWI, TH, TE, TR, TV, P, H1, 
     &                     H2, H3, H4, H5, H6)
C-----------------------------------------------------------------------
C     This subroutine computes the species enthalpy per unit
C     mole in thermal non-equilibrium.
C     -H1(I): total species enthalpy (including the formation enthalpy)
C     -H2(I): translational species enthalpy
C     -H3(I): electronic species enthalpy
C     -H4(I): rotational species enthalpy
C     -H5(I): vibrational species enthalpy
C     -H6(I): formation species enthalpy
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      INTEGER LWR1, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), TH, TE, TR, TV(1:NV), P, H1(1:NS), 
     &                 H2(1:NS), H3(1:NS), H4(1:NS), H5(1:NS), H6(1:NS)
C-----------------------------------------------------------------------
      IC1 = 0; IC2 = 0
      DO I = 1, NS
        IATOMICITY = WI(IATOMI+I-1)
        SELECT CASE(IATOMICITY)
C         Electron
          CASE(0)
            H2(I) = 2.5D0 *WR1(IUR) *TE
            H3(I) = 0.D0
            H4(I) = 0.D0
            H5(I) = 0.D0
C         Atom
          CASE(1)
            H2(I) = 2.5D0 *WR1(IUR) *TH
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
            H3(I) = WR1(IUR) *SUM1 /SUM2
            H4(I) = 0.D0
            H5(I) = 0.D0

C         Polyatomic molecule
          CASE DEFAULT
            H2(I) = 2.5D0 *WR1(IUR) *TH
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
            H3(I) = WR1(IUR) *SUM1 /SUM2
            IF (WI(ILINI+I-1)==1) THEN
              H4(I) = WR1(IUR) *TR
            ELSE
              H4(I) = 1.5D0 *WR1(IUR) *TR
            ENDIF
            SUM1 = 0.D0; NVIBMODE = WI(IVIBI+I-1)
            DO J = 1, NVIBMODE
              SUM1 = SUM1 + WR1(ITVIK+IC2+J-1) 
     &               /(DEXP(WR1(ITVIK+IC2+J-1) /TV(IC2+J)) -1.D0)
            ENDDO          
            IC2 = IC2 + NVIBMODE
            H5(I) = SUM1 *WR1(IUR)
        END SELECT
        H6(I) = WR1(IHFORI+I-1)
        H1(I) = H2(I) +H3(I) +H4(I) +H5(I) +H6(I)
      ENDDO
 
      END SUBROUTINE ENTHALPY 
C---------------------------------------------------------------------- 
C-----------------------------------------------------------------------
      SUBROUTINE ENTHALPYMASS (WR1, LWR1, WI, LWI, TH, TE, TR, TV, P, 
     &                         H1, H2, H3, H4, H5, H6)
C-----------------------------------------------------------------------
C     This subroutine computes the species enthalpy per unit
C     mass in thermal non-equilibrium.
C     -H1(I): total species enthalpy (including the formation enthalpy)
C     -H2(I): translational species enthalpy
C     -H3(I): electronic species enthalpy
C     -H4(I): rotational species enthalpy
C     -H5(I): vibrational species enthalpy
C     -H6(I): formation species enthalpy
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      INTEGER LWR1, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), TH, TE, TR, TV(1:NV), P, H1(1:NS), 
     &                 H2(1:NS), H3(1:NS), H4(1:NS), H5(1:NS), H6(1:NS)
C-----------------------------------------------------------------------
      DOUBLE PRECISION MASS 
C-----------------------------------------------------------------------
      IC1 = 0; IC2 = 0
      DO I = 1, NS
        MASS = WR1(IMI+I-1)
        IATOMICITY = WI(IATOMI+I-1)
        SELECT CASE(IATOMICITY)
C         Electron
          CASE(0)
            H2(I) = 2.5D0 *WR1(IUR) *TE /MASS
            H3(I) = 0.D0
            H4(I) = 0.D0
            H5(I) = 0.D0
C         Atom
          CASE(1)
            H2(I) = 2.5D0 *WR1(IUR) *TH /MASS
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
            H3(I) = WR1(IUR) /MASS *SUM1 /SUM2 
            H4(I) = 0.D0
            H5(I) = 0.D0

C         Polyatomic molecule
          CASE DEFAULT
            H2(I) = 2.5D0 *WR1(IUR) *TH /MASS
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
            H3(I) = WR1(IUR) /MASS *SUM1 /SUM2
            IF (WI(ILINI+I-1)==1) THEN
              H4(I) = WR1(IUR) *TR /MASS
            ELSE
              H4(I) = 1.5D0 *WR1(IUR) *TR /MASS
            ENDIF
            SUM1 = 0.D0; NVIBMODE = WI(IVIBI+I-1)
            DO J = 1, NVIBMODE
              SUM1 = SUM1 + WR1(ITVIK+IC2+J-1) 
     &               /(DEXP(WR1(ITVIK+IC2+J-1) /TV(IC2+J)) -1.D0)
            ENDDO          
            IC2 = IC2 + NVIBMODE
            H5(I) = SUM1 *WR1(IUR) /MASS
        END SELECT
        H6(I) = WR1(IHFORI+I-1) /MASS
        H1(I) = H2(I) +H3(I) +H4(I) +H5(I) +H6(I)
      ENDDO
 
      END SUBROUTINE ENTHALPYMASS
C---------------------------------------------------------------------- 
