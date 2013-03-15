C-----------------------------------------------------------------------
      SUBROUTINE DENSITY (WR1, LWR1, X, ND, RHO)
C-----------------------------------------------------------------------
C     This subroutine computes the density of a mixture possibly
C     in thermal non-equilibrium for given molar fraction and number 
C     density. 
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1
      DOUBLE PRECISION WR1(1:LWR1), X(1:NS), ND, RHO
C-----------------------------------------------------------------------
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      RHO = 0.D0
      DO I = 1, NS
        RHO = RHO + X(I) *WR1(IMI+I-1) /WR1(IUNA)
      ENDDO
      
      RHO = RHO *ND

      END SUBROUTINE DENSITY 
C---------------------------------------------------------------------- 
C-----------------------------------------------------------------------
      SUBROUTINE NUMBERD (WR1, LWR1, P, TH, TE, X, ND)
C-----------------------------------------------------------------------
C     This subroutine computes the number density of a mixture possibly
C     in thermal non-equilibrium for given translational temperatures
C     and pressure.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1
      DOUBLE PRECISION WR1(1:LWR1), P, TH, TE, X(1:NS), ND
C-----------------------------------------------------------------------
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      ND = P /(WR1(IUKB) *TH *(1.D0 +X(1) *(TE /TH -1.D0)))     
 
      END SUBROUTINE NUMBERD 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE MOLARMASS (WR1, LWR1, RHO, ND, MM)
C-----------------------------------------------------------------------
C     This subroutine computes the molar mass of a mixture possibly
C     in thermal non-equilibrium for given mass and number density. 
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1
      DOUBLE PRECISION WR1(1:LWR1), RHO, ND, MM
C-----------------------------------------------------------------------
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      MM = RHO *WR1(IUNA) /ND

      END SUBROUTINE MOLARMASS 
C---------------------------------------------------------------------- 
C-----------------------------------------------------------------------
      SUBROUTINE MOLE2MASS (WR1, LWR1, ND, X, RHOI)
C-----------------------------------------------------------------------
C     This subroutine computes the particle mass density of the species
C     of a mixture possibly in thermal non-equilibrium for given molar 
C     fractions and mixture number number. 
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1
      DOUBLE PRECISION WR1(1:LWR1), ND, X(1:NS), RHOI(1:NS)
C-----------------------------------------------------------------------
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      DO IS = 1, NS
        RHOI(IS) = X(IS) *ND *WR1(IMI+IS-1)/WR1(IUNA)
      ENDDO

      END SUBROUTINE MOLE2MASS 
C----------------------------------------------------------------------- 
C-----------------------------------------------------------------------
      SUBROUTINE MASS2MOLE (WR1, LWR1, Y, X)
C-----------------------------------------------------------------------
C     This subroutine computes the mole fractions based on the mass
C     fractions.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1
      DOUBLE PRECISION WR1(1:LWR1), X(1:NS), Y(1:NS)
      DOUBLE PRECISION SUM1 
C-----------------------------------------------------------------------
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      SUM1 = 0.D0
      DO IS = 1, NS
        SUM1 = SUM1 +Y(IS) /WR1(IMI+IS-1)
      ENDDO
      DO IS = 1, NS
        X(IS) = Y(IS) /(WR1(IMI+IS-1) *SUM1)
      ENDDO
     

      END SUBROUTINE MASS2MOLE 
C----------------------------------------------------------------------- 
C-----------------------------------------------------------------------
      SUBROUTINE MASS2PRESSURE (WR1, LWR1, RHO, TH, TE, Y, P)
C-----------------------------------------------------------------------
C     This subroutine computes the pressure of a mixture possibly in 
C     in thermal non-equilibrium for given mass fractions. 
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1
      DOUBLE PRECISION WR1(1:LWR1), RHO, TH, TE, Y(1:NS), P
      DOUBLE PRECISION SUM1 
      INTEGER IS 
C-----------------------------------------------------------------------
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      SUM1 = 0.D0
      DO IS = 1, NS
        SUM1 = SUM1 + Y(IS) / WR1(IMI+IS-1)
      ENDDO
      P = (SUM1 +Y(1)/WR1(IMI) *(TE/TH -1.D0)) *RHO *TH *WR1(IUR)

      END SUBROUTINE MASS2PRESSURE
C---------------------------------------------------------------------- 
C-----------------------------------------------------------------------
      SUBROUTINE SETMASS (WR1, LWR1, MI)
C-----------------------------------------------------------------------
C     This subroutine retrieves the molar mass of the species
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1
      DOUBLE PRECISION WR1(1:LWR1), MI(1:NS)
      INTEGER IS 
C-----------------------------------------------------------------------
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      DO IS = 1, NS
        MI(IS) = WR1(IMI+IS-1)
      ENDDO

      END SUBROUTINE SETMASS
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE MOLE2ELEMENT (WI, LWI, X, XN)
C-----------------------------------------------------------------------
C     This subroutine computes the element mole fractions from the 
C     species mole fractions 
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWI, WI(1:LWI)
      DOUBLE PRECISION X(1:NS), XN(1:NC)
      DOUBLE PRECISION SUM1 
      INTEGER I, J
C-----------------------------------------------------------------------
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      SUM1 = 0.D0
      DO I = 1, NC
        XN(I) = 0.D0
        DO J = 1, NS
          XN(I) = XN(I) + WI(INUIK+(I-1)*NS+J-1) *X(J)
        ENDDO
        SUM1 = SUM1 +XN(I)
      ENDDO
      DO I = 1, NC
          XN(I) = XN(I) /SUM1
      ENDDO
      

      END SUBROUTINE MOLE2ELEMENT 
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
      SUBROUTINE MOLE2MASSFRAC (WR1, LWR1, X, Y)
C-----------------------------------------------------------------------
C     This subroutine computes the species mass fractions from the 
C     species mole fractions 
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1
      DOUBLE PRECISION WR1(1:LWR1), X(1:NS), Y(1:NS)
      DOUBLE PRECISION SUM1 
      INTEGER IS
C-----------------------------------------------------------------------
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      SUM1 = 0.D0
      DO IS = 1, NS
        SUM1 = SUM1 +X(IS) *WR1(IMI+IS-1)
      ENDDO
      DO IS = 1, NS
          Y(IS) = X(IS) *WR1(IMI+IS-1) /SUM1 
      ENDDO
      

      END SUBROUTINE MOLE2MASSFRAC 
C-----------------------------------------------------------------------


