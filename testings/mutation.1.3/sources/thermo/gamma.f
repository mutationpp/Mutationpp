C-----------------------------------------------------------------------
      SUBROUTINE FROZENGAMMA(WR1, LWR1, WI, LWI, T, P, X, EPS, GAMMA)
C-----------------------------------------------------------------------
C     This subroutine computes the ratio of the mixture frozen specific 
C     heat in thermal equilibrium. 
C     Two remarks valid for the frozen case:
C     -This ratio is identical computed per unit mole or unit mass.
C     -The specific heat at constant volume can be evaluated at constant 
C     pressure since the species energy only depends on temperature.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE '../general/memory.cmn'      
C-----------------------------------------------------------------------
      INTEGER LWR1, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), T, P, X(1:NS), EPS, 
     &                 GAMMA
C-----------------------------------------------------------------------
      INTEGER J
      DOUBLE PRECISION TV(1:NV), TP, TVP(1:NV) 
C-----------------------------------------------------------------------
      EXTERNAL ENERGY, ENTHALPY
C-----------------------------------------------------------------------
      DUM = 0.D0
      TP = T *(EPS +1.D0)
      TV(1) = 1.D0; TVP(1) = 1.D0
      DO J = 1, NV
        TV(J)  = T
        TVP(J) = TP
      ENDDO

C     Specific heat at constant pressure (per unit mole)
      CALL MIXPROPERTY (WR1, LWR1, WI, LWI, T, T, T, TV, DUM, X,
     &                  H1, H2, H3, H4, H5, H6, ENTHALPY)               

      CALL MIXPROPERTY (WR1, LWR1, WI, LWI, TP, TP, TP, TVP, DUM, X,
     &                  H1P, H2P, H3P, H4P, H5P, H6P, ENTHALPY)

C     Specific heat at constant volume (per unit mole)
      CALL MIXPROPERTY (WR1, LWR1, WI, LWI, T, T, T, TV, DUM, X,
     &                  E1, E2, E3, E4, E5, E6, ENERGY)               

      CALL MIXPROPERTY (WR1, LWR1, WI, LWI, TP, TP, TP, TVP, DUM, X,
     &                  E1P, E2P, E3P, E4P, E5P, E6P, ENERGY)

C     Ration of specific heats
      GAMMA = (H1P -H1) /(E1P -E1)

      END SUBROUTINE FROZENGAMMA 
C---------------------------------------------------------------------- 
C-----------------------------------------------------------------------
      SUBROUTINE EQUIGAMMA(WR1, LWR1, WI, LWI, T, P, RHO, XN, X, EPS, 
     &                     GAMMA, DRHODP)
C-----------------------------------------------------------------------
C     This subroutine computes the ratio of the mixture equilibrium 
C     specific heat in thermal equilibrium. 
C     Remark : This ratio is computed per unit mass 
C     (different per unit mole).
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE '../general/memory.cmn'      
C-----------------------------------------------------------------------
      INTEGER LWR1, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), T, P, RHO, X(1:NS), XN(1:NC), EPS, 
     &                 GAMMA, DRHODP
C-----------------------------------------------------------------------
      INTEGER J
      DOUBLE PRECISION TV(1:NV), TP, TVP(1:NV), XP(1:NS), MM, MMP, MMYP,
     &                 Y(1:NS), YP(1:NS), XOFYP(1:NS), XPP(1:NS), ND, PP
      EXTERNAL ENERGY, ENTHALPY
C-----------------------------------------------------------------------
      DUM = 0.D0
      TP = T *(EPS +1.D0)
      TV(1) = 1.D0; TVP(1) = 1.D0
      DO J = 1, NV
        TV(J)  = T
        TVP(J) = TP
      ENDDO
      PP = P *(EPS +1.D0)

C     Perturbed composition at constant pressure (molar fractions)
      CALL COMPOSITION (WR1, LWR1, WI, LWI, TP, P, XN, X, XP)
      MM = 0.D0; MMP = 0.D0 
      DO I = 1, NS
        MM  = MM  +X(I)  *WR1(IMI+I-1)
        MMP = MMP +XP(I) *WR1(IMI+I-1)
      ENDDO

C     Perturbed composition at constant temperature (molar fractions)
      CALL COMPOSITION (WR1, LWR1, WI, LWI, T, PP, XN, X, XPP)
      CALL NUMBERD (WR1, LWR1, PP, T, T, XPP, ND)
      CALL DENSITY (WR1, LWR1, XPP, ND, RHOP)

C     Perturbed composition at constant pressure (mass fractions)
      DO I = 1, NS
        Y(I) = X(I) *WR1(IMI+I-1) /MM 
      ENDDO
      CALL MASSCOMPOSITION (WR1, LWR1, WI, LWI, TP, RHO, XN, Y, YP)
      SUM = 0.D0
      DO I = 1, NS
        SUM = SUM +YP(I) /WR1(IMI+I-1)
      ENDDO
      MMYP = 0.D0 
      DO I = 1, NS
        XOFYP(I) = YP(I) /(SUM *WR1(IMI+I-1))
        MMYP = MMYP +XOFYP(I) *WR1(IMI+I-1)
      ENDDO

C     Specific heat:
C     -at constant pressure (per unit mass)
      CALL MIXPROPERTY (WR1, LWR1, WI, LWI, T, T, T, TV, DUM, X,
     &                  H1, H2, H3, H4, H5, H6, ENTHALPY)
      H1 = H1 /MM               
      CALL MIXPROPERTY (WR1, LWR1, WI, LWI, TP, TP, TP, TVP, DUM, XP,
     &                  H1P, H2P, H3P, H4P, H5P, H6P, ENTHALPY)
      H1P = H1P /MMP               

C     -at constant volume (per unit mass)
      CALL MIXPROPERTY (WR1, LWR1, WI, LWI, T, T, T, TV, DUM, X,
     &                  E1, E2, E3, E4, E5, E6, ENERGY)               
      E1 = E1 /MM
      CALL MIXPROPERTY (WR1, LWR1, WI, LWI, TP, TP, TP, TVP, DUM, XOFYP,
     &                  E1P, E2P, E3P, E4P, E5P, E6P, ENERGY)
      E1P = E1P /MMYP

C     Ratio of specific heats
      GAMMA = (H1P -H1) /(E1P -E1)

C     Partial derivative of density with respect to pressure at 
C     constant temperature
      DRHODP = (RHOP -RHO) /(P *EPS)

      END SUBROUTINE EQUIGAMMA 
C---------------------------------------------------------------------- 
