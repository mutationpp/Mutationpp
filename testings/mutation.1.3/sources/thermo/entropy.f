C-----------------------------------------------------------------------
      SUBROUTINE ENTROPY (WR1, LWR1, WI, LWI, TH, TE, TR, TV, P, S1, S2,
     &                    S3, S4, S5, S6)
C-----------------------------------------------------------------------
C     This subroutine computes the species entropy per unit
C     mole in thermal non-equilibrium.
C     -S1(I): total species entropy 
C     -S2(I): translational species entropy
C     -S3(I): electronic species entropy
C     -S4(I): rotational species entropy
C     -S5(I): vibrational species entropy
C     -S6(I): THE FORMATION ENTHALPY is null at 0K... 
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      INTEGER LWR1, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), TH, TE, TR, TV(1:NV), P, S1(1:NS),
     &                 S2(1:NS), S3(1:NS), S4(1:NS), S5(1:NS), S6(1:NS)
C-----------------------------------------------------------------------
      DOUBLE PRECISION LNTE, LNTH 
C-----------------------------------------------------------------------
      LNTH = DLOG(TH); LNTE = DLOG(TE)
      IC1 = 0; IC2 = 0; IC3 = 0
      DO I = 1, NS
        IATOMICITY = WI(IATOMI+I-1)
        SELECT CASE(IATOMICITY)
C         Electron Beware! Don't forget the spin of the electron...
          CASE(0)
            S2(I) = WR1(IUR) *(2.5D0 *(1.D0 +LNTE) -DLOG(P) 
     &              +WR1(IUQT) +1.5D0*DLOG(WR1(IMI+I-1))) 
            S3(I) = WR1(IUR) *DLOG(2.D0)
            S4(I) = 0.D0
            S5(I) = 0.D0
C         Atom
          CASE(1)
            S2(I) = WR1(IUR) *(2.5D0 *(1.D0 +LNTH) -DLOG(P) 
     &              +WR1(IUQT) +1.5D0*DLOG(WR1(IMI+I-1))) 
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
            S3(I) = WR1(IUR) *(DLOG(SUM2) +SUM1 /(SUM2 *TE))
            S4(I) = 0.D0
            S5(I) = 0.D0

C         Polyatomic molecule
          CASE DEFAULT
            S2(I) = WR1(IUR) *(2.5D0 *(1.D0 +LNTH) -DLOG(P) 
     &              +WR1(IUQT) +1.5D0*DLOG(WR1(IMI+I-1))) 
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
            S3(I) = WR1(IUR) *(DLOG(SUM2) +SUM1 /(SUM2 *TE))
            IF (WI(ILINI+I-1)==1) THEN
              S4(I) = WR1(IUR) *(1.D0 +DLOG(TR /WR1(ITRIK+IC3))
     &                -DLOG(1.D0 *WI(ISYMI+I-1)))
            ELSE
              S4(I) = WR1(IUR) *(1.5D0 +1.5 *DLOG(TR /WR1(ITRIK+IC3))
     &                -DLOG(1.D0 *WI(ISYMI+I-1)))
            ENDIF
            IC3 =IC3 +1
            SUM1 = 0.D0; NVIBMODE = WI(IVIBI+I-1)
            DO J = 1, NVIBMODE
              SUM1 = SUM1 + WR1(ITVIK+IC2+J-1) /TV(IC2+J)  
     &               /(DEXP(WR1(ITVIK+IC2+J-1) /TV(IC2+J)) -1.D0)
     &               -DLOG(1.D0 -DEXP(-WR1(ITVIK+IC2+J-1) /TV(IC2+J)))
            ENDDO          
            IC2 = IC2 + NVIBMODE
            S5(I) = SUM1 *WR1(IUR)
        END SELECT
        S6(I) = 0.D0 
        S1(I) = S2(I) +S3(I) +S4(I) +S5(I)
      ENDDO
 
      END SUBROUTINE ENTROPY 
C---------------------------------------------------------------------- 
