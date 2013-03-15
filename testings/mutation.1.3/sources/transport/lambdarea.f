C-----------------------------------------------------------------------
      SUBROUTINE LAMBDAREAYOS (WI, LWI, WR1, LWR1, WR2, LWR2, H, T, X, 
     &                         LAMBDA)
C-----------------------------------------------------------------------
C     This subroutine computes the reactive thermal conductivity using 
C     Butler and Brokaw's formula in thermal equilibrium and Yos' rule.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2, LWI, WI(1:LWI) 
      DOUBLE PRECISION WR1(1:LWR1), WR2(1:LWR2), X(1:NS), LAMBDA,
     &                 H(1:NS), T
C-----------------------------------------------------------------------
      DOUBLE PRECISION  G(1:NR), BETA(1:NR)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      CALL LAMBDAREASYSTYOS (NS, NR, WI(INUIK), WR2(IBINIJ), WR1(IUNA),
     &                       X, H, G, BETA)

      LAMBDA = 0.D0
      DO I = 1, NR
        LAMBDA = LAMBDA +BETA(I) /G(I) *BETA(I) /(WR1(IUKB) *T *T)
      ENDDO

      END SUBROUTINE LAMBDAREAYOS
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE LAMBDAREASYSTYOS (NS, NR, NUIK, BINIJ, NA, X, H, G, 
     &                             BETA)
C-----------------------------------------------------------------------
C     This subroutine computes the system for the reactive thermal 
C     conductivity using Butler and Brokaw's formula in thermal 
C     equilibrium and Yos' rule.
C-----------------------------------------------------------------------
      IMPLICIT NONE 
      INTEGER NS, NR, NUIK(1:NS*NS)
      DOUBLE PRECISION BINIJ(1:(NS+1)*NS/2), NA, X(1:NS), H(1:NS),
     &                 G(1:NR), BETA(1:NR)
C-----------------------------------------------------------------------
      INTEGER NC, I, J, IJ, K, L, KL
C-----------------------------------------------------------------------
      NC = NS -NR
      DO I = 1, NR
        BETA(I) = 0.D0
        DO J = 1, NC
          BETA(I) = BETA(I) +NUIK((NC+I-1)*NS+J) *H(J) /NA
        ENDDO
        BETA(I)   = BETA(I) +NUIK((NC+I-1)*NS+NC+I) *H(NC+I) /NA
      ENDDO
      DO I = 1, NR
        G(I) = 0.D0
        DO K = 1, NS-1
          DO L = K+1, NS
            KL = ((K-1)*(2*NS-K)+2*L)/2 
            G(I) = G(I) +X(K) *X(L) /BINIJ(KL) 
     &             *( NUIK((NC+I-1)*NS+K) /X(K) 
     &               -NUIK((NC+I-1)*NS+L) /X(L))**2 
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE LAMBDAREASYSTYOS
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE LAMBDAREACG (WI, LWI, WR1, LWR1, WR2, LWR2, H, T, X, 
     &                       LAMBDA)
C-----------------------------------------------------------------------
C     This subroutine computes the reactive thermal conductivity 
C     using Butler and Brokaw's formula in thermal equilibrium.
C     The conjugate gradient with preconditioning is used to solve the 
C     linear system.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2, LWI, WI(1:LWI) 
      DOUBLE PRECISION  WR1(1:LWR1), WR2(1:LWR2), X(1:NS), LAMBDA,
     &                  H(1:NS), T
C-----------------------------------------------------------------------
      DOUBLE PRECISION  G(1:NR*(NR+1)/2), BETA(1:NR), ALPHA(1:NR),
     &                  TEMP(1:NR), DMI(1:NR), ZN(1:NR), RN(1:NR)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      CALL LAMBDAREASYST (NS, NR, WI(INUIK), WR2(IBINIJ), WR1(IUNA), X,
     &                    H, G, BETA)
      DO I = 1, NR
        RN(I) = BETA(I)
      ENDDO
      
      ITERMAX = 2 
      CALL EGSCG1 (NR, G, DMI, ALPHA, ZN, RN, TEMP, ITERMAX)    

      LAMBDA = 0.D0
      DO I = 1, NR
        LAMBDA = LAMBDA +BETA(I) *ALPHA(I) /(WR1(IUKB) *T *T)
      ENDDO

      END SUBROUTINE LAMBDAREACG
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE LAMBDAREAD (WI, LWI, WR1, LWR1, WR2, LWR2, H, T, X, 
     &                       LAMBDA)
C-----------------------------------------------------------------------
C     This subroutine computes the reactive thermal conductivity 
C     using Butler and Brokaw's formula in thermal equilibrium.
C     A direct method is used to solve the linear system.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2, LWI, WI(1:LWI) 
      DOUBLE PRECISION  WR1(1:LWR1), WR2(1:LWR2), X(1:NS), LAMBDA,
     &                  H(1:NS), T
C-----------------------------------------------------------------------
      DOUBLE PRECISION  G(1:NR*(NR+1)/2), BETA(1:NR), ALPHA(1:NR),
     &                  TEMP(1:NR)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      CALL LAMBDAREASYST (NS, NR, WI(INUIK), WR2(IBINIJ), WR1(IUNA), X,
     &                    H, G, BETA)
      DO I = 1, NR
        ALPHA(I) = BETA(I)
      ENDDO

      CALL EGSDEC (NR, G, TEMP, IER)
      CALL EGSSOL (NR, G, ALPHA)

      LAMBDA = 0.D0
      DO I = 1, NR
        LAMBDA = LAMBDA +BETA(I) *ALPHA(I) /(WR1(IUKB) *T *T)
      ENDDO

      END SUBROUTINE LAMBDAREAD
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE LAMBDAREASYST (NS, NR, NUIK, BINIJ, NA, X, H, G, BETA)
C-----------------------------------------------------------------------
C     This subroutine computes the system for the reactive thermal 
C     conductivity using Butler and Brokaw's formula in thermal 
C     equilibrium.
C-----------------------------------------------------------------------
      IMPLICIT NONE 
      INTEGER NS, NR, NUIK(1:NS*NS)
      DOUBLE PRECISION BINIJ(1:(NS+1)*NS/2), NA, X(1:NS), H(1:NS),
     &                 G(1:(NR+1)*NR/2), BETA(1:NR)
C-----------------------------------------------------------------------
      INTEGER NC, I, J, IJ, K, L, KL
C-----------------------------------------------------------------------
      NC = NS -NR
      DO I = 1, NR
        BETA(I) = 0.D0
        DO J = 1, NC
          BETA(I) = BETA(I) +NUIK((NC+I-1)*NS+J) *H(J) /NA
        ENDDO
        BETA(I)   = BETA(I) +NUIK((NC+I-1)*NS+NC+I) *H(NC+I) /NA
      ENDDO
      DO I = 1, NR
        DO J = I, NR
          IJ = ((I-1)*(2*NR-I)+2*J)/2 
          G(IJ) = 0.D0
          DO K = 1, NS-1
            DO L = K+1, NS
              KL = ((K-1)*(2*NS-K)+2*L)/2 
              G(IJ) = G(IJ) +X(K) *X(L) /BINIJ(KL) 
     &                *( NUIK((NC+I-1)*NS+K) /X(K) 
     &                  -NUIK((NC+I-1)*NS+L) /X(L)) 
     &                *( NUIK((NC+J-1)*NS+K) /X(K) 
     &                  -NUIK((NC+J-1)*NS+L) /X(L)) 
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE LAMBDAREASYST
C-----------------------------------------------------------------------

