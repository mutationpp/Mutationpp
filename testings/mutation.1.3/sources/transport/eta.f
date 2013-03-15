C-----------------------------------------------------------------------
      SUBROUTINE ETAWILKE (WR1, LWR1, WR2, LWR2, X, ETA)
C-----------------------------------------------------------------------
C     This subroutine computes the mixture viscosity. The electron 
C     contribution is neglected. Wilke's mixture rule is used.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2 
      DOUBLE PRECISION  WR1(1:LWR1), WR2(1:LWR2), X(1:NS), ETA
C-----------------------------------------------------------------------
      DOUBLE PRECISION GI(1:NS-NE)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      N = NS-NE
      CALL SYSTWILKE (N, WR1(IMI+NE), X(NE+1), GI, WR2(IETAI))

      ETA = 0.D0
      DO I = 1, N
        ETA = ETA + X(NE+I) *WR2(IETAI+I-1) /GI(I)
      ENDDO

      END SUBROUTINE ETAWILKE
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SYSTWILKE (N, MI, X, GI, MUI)
C-----------------------------------------------------------------------
C     This subroutine computes the system for the mixture viscosity 
C     or heavy particle thermal conductivity using Wilke's rule. 
C-----------------------------------------------------------------------
      IMPLICIT NONE 
      INTEGER N
      DOUBLE PRECISION MI(1:N), X(1:N), GI(1:N), MUI(1:N)
C-----------------------------------------------------------------------
      INTEGER I, J, IJ
      DOUBLE PRECISION RATIO, A 
C-----------------------------------------------------------------------
C     3*N**2 Square roots to compute, very expensive...

      DO I = 1, N
        GI(I) = 0.D0
        DO J = 1, N
          RATIO = MI(I) /MI(J)
          A     = 1.D0 +DSQRT(MUI(I) /MUI(J) /DSQRT(RATIO))
          GI(I) = GI(I) +X(J) *A *A /DSQRT(8.D0 *(1.D0 +RATIO)) 
        ENDDO
      ENDDO

      END SUBROUTINE SYSTWILKE
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE ETAYOS (WR1, LWR1, WR2, LWR2, X, ETA)
C-----------------------------------------------------------------------
C     This subroutine computes the mixture viscosity. The electron 
C     contribution is neglected. Yos'mixture rule is used. 
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2 
      DOUBLE PRECISION  WR1(1:LWR1), WR2(1:LWR2), X(1:NS), ETA
C-----------------------------------------------------------------------
      INTEGER N
      DOUBLE PRECISION GIJ(1:(NS-NE+1)*(NS-NE)/2), GII(1:NS-NE)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      N = NS-NE
      CALL SYSTETAYOS (N, WR1(IMI+NE), WR2(IBINIJ+NE*NS), WR2(IAIJ),
     &                 WR1(IUNA), X(NE+1), GIJ, GII)


      SUM1 = 0.D0; SUM2 = 0.D0
      DO I = 1, N-1
        DO J = I+1, N
          IJ   = ((I-1)*(2*N-I)+2*J)/2
          FAC  = 2.D0 *X(NE+I) *X(NE+J) *(1.D0 /GII(I) -1.D0 /GII(J))**2
          SUM1 = SUM1 + FAC  
          SUM2 = SUM2 + FAC *GIJ(IJ)
        ENDDO
      ENDDO
      GAV = SUM2 / SUM1

      SUM = 0.D0
      DO I = 1, N
        SUM = SUM + X(NE+I) /(GII(I) +GAV)
      ENDDO

      ETA = SUM /(1.D0 -GAV *SUM)

      END SUBROUTINE ETAYOS
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SYSTETAYOS (N, MI, BINIJ, AIJ, NA, X, GIJ, GII)
C-----------------------------------------------------------------------
C     This subroutine computes the system for the mixture viscosity 
C     using Yos' rule. 
C-----------------------------------------------------------------------
      IMPLICIT NONE 
      INTEGER N
      DOUBLE PRECISION MI(1:N), BINIJ(1:(N+1)*N/2), AIJ(1:(N+1)*N/2), 
     &                 NA, X(1:N), GIJ(1:(N+1)*N/2), GII(1:N)
C-----------------------------------------------------------------------
      DOUBLE PRECISION B(1:(N+1)*N/2)
      INTEGER I, J, IJ
C-----------------------------------------------------------------------
      DO I = 1, N
        DO J = I, N
          IJ = ((I-1)*(2*N-I)+2*J)/2
          GIJ(IJ) = 2.D0 /(MI(I) +MI(J)) *NA /BINIJ(IJ) 
     &            *(1.D0 -0.6D0 *AIJ(IJ))
          B(IJ) = AIJ(IJ) /BINIJ(IJ)
        ENDDO
      ENDDO

      CALL SYMMATVEC (N, B, X, GII)
      DO I = 1, N
        GII(I) = GII(I) *1.2D0 *NA /MI(I)
      ENDDO

      END SUBROUTINE SYSTETAYOS
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE ETAD (WR1, LWR1, WR2, LWR2, X, ETA)
C-----------------------------------------------------------------------
C     This subroutine computes the mixture viscosity. The electron 
C     contribution is neglected. A direct method is used to solve the 
C     linear system.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2 
      DOUBLE PRECISION  WR1(1:LWR1), WR2(1:LWR2), X(1:NS), ETA
C-----------------------------------------------------------------------
      DOUBLE PRECISION G(1:(NS-NE+1)*(NS-NE)/2), BETA(1:NS-NE),
     &                 ALPHA(1:NS-NE), TEMP(1:NS-NE)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      CALL ETASYST (NS-NE, WR1(IMI+NE), WR2(IBINIJ+NE*NS), WR2(IAIJ), 
     &             WR2(IETAI), WR1(IUNA), X(NE+1), G, BETA)

      DO I = 1, NS-NE
        ALPHA(I) = BETA(I)
      ENDDO
      
      CALL EGSDEC (NS-NE, G, TEMP, IER)
      CALL EGSSOL (NS-NE, G, ALPHA)

      ETA = 0.D0
      DO I = 1, NS-NE
        ETA = ETA +BETA(I) *ALPHA(I)
      ENDDO

      END SUBROUTINE ETAD
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE ETACG (WR1, LWR1, WR2, LWR2, X, ETA)
C-----------------------------------------------------------------------
C     This subroutine computes the mixture viscosity. The electron 
C     contribution is neglected. The conjugate gradient with 
C     preconditioning is used to solve the linear system.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2 
      DOUBLE PRECISION  WR1(1:LWR1), WR2(1:LWR2), X(1:NS), ETA
C-----------------------------------------------------------------------
      DOUBLE PRECISION G(1:(NS-NE+1)*(NS-NE)/2), BETA(1:NS-NE),
     &                 ALPHA(1:NS-NE), RN(1:NS-NE), P(1:NS-NE)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      CALL ETASYST (NS-NE, WR1(IMI+NE), WR2(IBINIJ+NE*NS), WR2(IAIJ), 
     &             WR2(IETAI), WR1(IUNA), X(NE+1), G, BETA)

      DO I = 1, NS-NE
        RN(I) = BETA(I)
      ENDDO

      DO I = 1, NS-NE
        II = ((I-1)*(2*(NS-NE)-I)+2*I)/2
        P(I) = G(II)
      ENDDO
 
      ITERMAX = 2; TOLRES = 6.D-2
      CALL CG (NS-NE, G, P, ALPHA, RN, TOLRES, ITERMAX, ITER)

      ETA = 0.D0
      DO I = 1, NS-NE
        ETA = ETA +BETA(I) *ALPHA(I)
      ENDDO

      END SUBROUTINE ETACG
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE ETASYST (N, MI, BINIJ, AIJ, ETAI, NA, X, G, BETA)
C-----------------------------------------------------------------------
C     This subroutine computes the system for the mixture viscosity. 
C-----------------------------------------------------------------------
      IMPLICIT NONE 
      INTEGER N
      DOUBLE PRECISION MI(1:N), BINIJ(1:(N+1)*N/2), AIJ(1:(N+1)*N/2), 
     &                 ETAI(1:N), NA, X(1:N), G(1:(N+1)*N/2), BETA(1:N)
C-----------------------------------------------------------------------
      DOUBLE PRECISION FAC
      INTEGER I, J, II, IJ, JJ
C-----------------------------------------------------------------------
      DO I = 1, N
        BETA(I) =  X(I)
        II = ((I-1)*(2*N-I)+2*I)/2 
        G(II) = X(I) * X(I) /ETAI(I)
      ENDDO
      DO I = 1, N
        DO J = I+1, N
          IJ = ((I-1)*(2*N-I)+2*J)/2 
          II = ((I-1)*(2*N-I)+2*I)/2 
          JJ = ((J-1)*(2*N-J)+2*J)/2 
          FAC = 2.D0 *X(I) *X(J) /(BINIJ(IJ) *(MI(I) +MI(J)) /NA)
          G(IJ) = FAC *(-1.D0 +0.6D0 *AIJ(IJ))
          G(II) = G(II) + FAC* (1.D0 +0.6D0 * MI(J) /MI(I) *AIJ(IJ))
          G(JJ) = G(JJ) + FAC* (1.D0 +0.6D0 * MI(I) /MI(J) *AIJ(IJ))
        ENDDO
      ENDDO

      END SUBROUTINE ETASYST
C-----------------------------------------------------------------------

