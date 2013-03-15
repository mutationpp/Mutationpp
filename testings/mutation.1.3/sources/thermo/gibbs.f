C-----------------------------------------------------------------------
      SUBROUTINE GIBBS (WR1, LWR1, WI, LWI, TH, TE, TR, TV, P, G)
C-----------------------------------------------------------------------
C     This subroutine computes the species Gibbs free energy per unit
C     mole in thermal non-equilibrium. 
C     The chemical potential may be obtained by setting P = 1.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      INTEGER LWR1, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), TH, TE, TR, TV(1:NV), P, G(1:NS)
C-----------------------------------------------------------------------
      DOUBLE PRECISION GTRA, GROT, GVIB, GELE, LNTH, LNTE, LNP
C-----------------------------------------------------------------------
      LNTH = DLOG(TH)
      LNTE = DLOG(TE)
      LNP  = DLOG(P)
      IC1 = 0; IC2 = 0; IC3 = 0
      DO I = 1, NS
        IATOMICITY = WI(IATOMI+I-1)
        SELECT CASE(IATOMICITY)

C         Electron Beware! Don't forget the spin of the electron...
          CASE(0)
            GTRA = WR1(IUR)*TE*(-2.5D0*LNTE+LNP
     &           -1.5D0*DLOG(WR1(IMI+I-1))-WR1(IUQT))
            GROT = 0.D0; GVIB = 0.D0
            GELE = -DLOG(2.D0) *WR1(IUR) *TE
C         Atom
          CASE(1)
            GTRA = WR1(IUR)*TH*(-2.5D0*LNTH+LNP
     &             -1.5D0*DLOG(WR1(IMI+I-1))-WR1(IUQT))
            GELE = 0.D0; NELE = WI(IELEI+I-1)
            IF (NELE/=0) THEN
              DO J = 1, NELE
                GELE = GELE + WI(IEGIK+IC1+J-1)
     &                 *DEXP(-WR1(IEEIK+IC1+J-1) /TE)
              ENDDO
            ELSE
              GELE = 1.D0
            ENDIF
            GELE = -DLOG(GELE) *WR1(IUR) *TE
            IC1 = IC1 +NELE
            GROT = 0.D0; GVIB = 0.D0

C         Polyatomic molecule
          CASE DEFAULT
            GTRA = WR1(IUR)*TH*(-2.5D0*LNTH+LNP
     &           -1.5D0*DLOG(WR1(IMI+I-1))-WR1(IUQT))
            GELE = 0.D0; NELE = WI(IELEI+I-1)
            IF (NELE/=0) THEN
              DO J = 1, NELE
                GELE = GELE + WI(IEGIK+IC1+J-1)
     &                 *DEXP(-WR1(IEEIK+IC1+J-1) /TE)
              ENDDO
            ELSE
              GELE = 1.D0
            ENDIF
            GELE = -DLOG(GELE) *WR1(IUR) *TE
            IC1 = IC1 +NELE
            IF (WI(ILINI+I-1)==1) THEN
              GROT = -WR1(IUR) *TR *(DLOG(TR /WR1(ITRIK+IC3))
     &                -DLOG(1.D0 *WI(ISYMI+I-1)))
            ELSE
              GROT = -WR1(IUR) *TR *(1.5 *DLOG(TR /WR1(ITRIK+IC3)) 
     &                -DLOG(1.D0 *WI(ISYMI+I-1)))
            ENDIF
            IC3 =IC3 +1
            GVIB = 0.D0; NVIBMODE = WI(IVIBI+I-1)
            DO J = 1, NVIBMODE
              GVIB = GVIB +DLOG(1.D0 -DEXP(-WR1(ITVIK+IC2+J-1) 
     &               /TV(IC2+J))) *WR1(IUR) *TV(IC2+J) 
            ENDDO          
            IC2 = IC2 + NVIBMODE
        END SELECT
        G(I) = GTRA + GELE + GROT + GVIB +WR1(IHFORI+I-1)
      ENDDO

      END SUBROUTINE GIBBS
C---------------------------------------------------------------------- 
