C-----------------------------------------------------------------------
      SUBROUTINE PHIE (WR1, LWR1, WR2, LWR2, X, TE, FIJ, N)
C-----------------------------------------------------------------------
C     This subroutine computes the correction functions of the electron 
C     gas (N = Lagrange-Sonine polynomial order).
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2 
      DOUBLE PRECISION WR1(1:LWR1), WR2(1:LWR2), X(1:NS), TE,
     &                 FIJ(1:NS-NE)
C-----------------------------------------------------------------------
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      SELECT CASE(N)
        CASE(1)
          DO I = 1, NS -NE
            FIJ(I) = 1.D0
          ENDDO    

        CASE(2)
          DENO = WR2(IDQ11EE)
          DO I = 1, NS -NE
            FIJ(I) = 2.D0 /(3.D0 *X(I+1)) * WR2(IBINIJ+I)
     &               *DSQRT(WR1(IMI) /(2.D0 *WR1(IUPI) *WR1(IUR) *TE)) 
     &               *WR2(IDQ10EI+I-1) *WR2(IDQ10EE) /DENO
            FIJ(I) = 1.D0 /(1.D0 +FIJ(I))
          ENDDO    

        CASE(3)
          DENO = WR2(IDQ11EE) *WR2(IDQ22EE) -WR2(IDQ12EE)**2 
          DO I = 1, NS -NE
            BETA1 = (WR2(IDQ10EI+I-1) *WR2(IDQ22EE) -WR2(IDQ12EE) 
     &              *WR2(IDQ20EI+I-1)) /DENO
            BETA2 = (WR2(IDQ20EI+I-1) *WR2(IDQ11EE) -WR2(IDQ12EE) 
     &              *WR2(IDQ10EI+I-1)) /DENO
            FIJ(I) = 2.D0 /(3.D0 *X(I+1)) * WR2(IBINIJ+I)
     &               *DSQRT(WR1(IMI) /(2.D0 *WR1(IUPI) *WR1(IUR) *TE)) 
     &               *(WR2(IDQ10EE) *BETA1 +WR2(IDQ20EE) *BETA2)
            FIJ(I) = 1.D0 /(1.D0 +FIJ(I))
          ENDDO    

        CASE DEFAULT
            WRITE(*,*) 'Lagrange-Sonine order not implemented' 

      END SELECT  

      END SUBROUTINE PHIE
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE TDIFE (WR1, LWR1, WR2, LWR2, X, TE, CHI, N)
C-----------------------------------------------------------------------
C     This subroutine computes the thermal diffusion ratios of the 
C     electron gas (N = Lagrange-Sonine polynomial order).
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2 
      DOUBLE PRECISION WR1(1:LWR1), WR2(1:LWR2), X(1:NS), TE,
     &                 CHI(1:NS)
C-----------------------------------------------------------------------
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      CHI(1) = 0.D0
      SELECT CASE(N)
        CASE(1)
          DO I = 1, NS -NE
            CHI(I+1) = 0.D0 
          ENDDO    

        CASE(2)
          DENO = WR2(IDQ11EE)
          DO I = 1, NS -NE
            CHI(I+1) = 5.D0 /2.D0 *WR2(IDQ10EI+I-1)*X(1) /DENO
            CHI(1)   = CHI(1) -CHI(I+1)
          ENDDO    

        CASE(3)
          DENO = WR2(IDQ11EE) *WR2(IDQ22EE) -WR2(IDQ12EE)**2 
          DO I = 1, NS -NE
            BETA1 = (WR2(IDQ10EI+I-1) *WR2(IDQ22EE) -WR2(IDQ12EE) 
     &              *WR2(IDQ20EI+I-1)) /DENO
            BETA2 = (WR2(IDQ20EI+I-1) *WR2(IDQ11EE) -WR2(IDQ12EE) 
     &              *WR2(IDQ10EI+I-1)) /DENO
            CHI(I+1) = 5.D0 /2.D0 *X(1) *BETA1
            CHI(1)   = CHI(1) -CHI(I+1)
          ENDDO    

        CASE DEFAULT
            WRITE(*,*) 'Lagrange-Sonine order not implemented' 

      END SELECT  

      END SUBROUTINE TDIFE
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE LAMBDAE (WR1, LWR1, WR2, LWR2, X, TE, LAMBDA, N)
C-----------------------------------------------------------------------
C     This subroutine computes the translational thermal conductivity 
C     of the electron gas (N = Lagrange-Sonine polynomial order).
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2 
      DOUBLE PRECISION WR1(1:LWR1), WR2(1:LWR2), X(1:NS), TE, LAMBDA
C-----------------------------------------------------------------------
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      SELECT CASE(N)
        CASE(1)
            A = 0.D0 
        CASE(2)
            A = 1.D0 /WR2(IDQ11EE)
        CASE(3)
            A = 1.D0 /(WR2(IDQ11EE) - WR2(IDQ12EE)**2 /WR2(IDQ22EE))
        CASE DEFAULT
            WRITE(*,*) 'Lagrange-Sonine order not implemented' 
      END SELECT  

      LAMBDA = 75.D0 /8.D0 *WR1(IUKB) *X(1) *DSQRT(2.D0 *WR1(IUPI) 
     &         *WR1(IUR) *TE / WR1(IMI) ) *A

      END SUBROUTINE LAMBDAE
C-----------------------------------------------------------------------
