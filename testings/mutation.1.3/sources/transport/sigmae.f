C-----------------------------------------------------------------------
      SUBROUTINE SIGMAE (WR1, LWR1, WR2, LWR2, X, TE, SIGMA, N)
C-----------------------------------------------------------------------
C     This subroutine computes the electrical conductivity
C     of the electron gas (N = Lagrange-Sonine polynomial order).
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2
      DOUBLE PRECISION WR1(1:LWR1), WR2(1:LWR2), X(1:NS), TE, SIGMA
C-----------------------------------------------------------------------
      DOUBLE PRECISION G(1:N*(N+1)/2), BETA(1:N), ALPHA(1:N), TEMP(1:N)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      BETA(1) = 3.D0 /X(1) *DSQRT(2.D0 *WR1(IUPI) *WR1(IUR) *TE
     &          /WR1(IMI))
      G(1) = WR2(IDQ00EE)

      SELECT CASE(N)
        CASE(1)
          ALPHA(1) = BETA(1) /G(1)

        CASE(2)
          G(2) = WR2(IDQ10EE)
          G(3) = WR2(IDQ11EE)
          ALPHA(1) = BETA(1) /(G(1) -G(2)**2 /G(3))

        CASE(3)
          G(2) = WR2(IDQ10EE)
          G(3) = WR2(IDQ20EE)
          G(4) = WR2(IDQ11EE)
          G(5) = WR2(IDQ12EE)
          G(6) = WR2(IDQ22EE)
          BETA(2) = 0.D0
          BETA(3) = 0.D0
          DO I = 1, N
             ALPHA(I) = BETA(I)
          ENDDO
          CALL EGSDEC (N, G, TEMP, IER)
          CALL EGSSOL (N, G, ALPHA)
        CASE DEFAULT
          WRITE(*,*) 'Lagrange-Sonine order not implemented'

      END SELECT

      SIGMA = X(1)**2 *WR1(IUE)**2 /(2.D0 *WR1(IUKB) *TE) *ALPHA(1)


      END SUBROUTINE SIGMAE
C-----------------------------------------------------------------------

