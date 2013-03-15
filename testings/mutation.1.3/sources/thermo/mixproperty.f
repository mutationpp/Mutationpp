C-----------------------------------------------------------------------
      SUBROUTINE MIXPROPERTY (WR1, LWR1, WI, LWI, TH, TE, TR, TV, P, X, 
     &                        PTY1, PTY2, PTY3, PTY4, PTY5, PTY6,
     &                        PROPERTY)
C-----------------------------------------------------------------------
C     This subroutine computes the mixture property per unit mole in 
C     thermal non-equilibrium. The subroutine PROPERTY is an external.
C     -PTY1: total mixture property (including PTY6)
C     -PTY2: translational mixture property
C     -PTY3: electronic mixture property
C     -PTY4: rotational mixture property
C     -PTY5: vibrational mixture property
C     -PTY6: formation mixture property 
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      INTEGER LWR1, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), TH, TE, TR, TV(1:NV), P, X(1:NS), 
     &                 PTY1, PTY2, PTY3, PTY4, PTY5, PTY6
C-----------------------------------------------------------------------
      DOUBLE PRECISION SPTY1(1:NS), SPTY2(1:NS), SPTY3(1:NS), 
     &                 SPTY4(1:NS), SPTY5(1:NS), SPTY6(1:NS)
C-----------------------------------------------------------------------
      CALL PROPERTY (WR1, LWR1, WI, LWI, TH, TE, TR, TV, P, SPTY1, 
     &               SPTY2, SPTY3, SPTY4, SPTY5, SPTY6)

      PTY1 = 0.D0; PTY2 = 0.D0; PTY3 = 0.D0
      PTY4 = 0.D0; PTY5 = 0.D0; PTY6 = 0.D0
      DO I = 1, NS
        PTY1 = PTY1 + X(I) *SPTY1(I)
        PTY2 = PTY2 + X(I) *SPTY2(I)
        PTY3 = PTY3 + X(I) *SPTY3(I)
        PTY4 = PTY4 + X(I) *SPTY4(I)
        PTY5 = PTY5 + X(I) *SPTY5(I)
        PTY6 = PTY6 + X(I) *SPTY6(I)
      ENDDO

      END SUBROUTINE MIXPROPERTY
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE MIXENTROPY (WR1, LWR1, WI, LWI, TH, TE, TR, TV, P, X, 
     &                       PTY1, PTY2, PTY3, PTY4, PTY5, PTY6)
C-----------------------------------------------------------------------
C     This subroutine computes the mixture property per unit mole in 
C     thermal non-equilibrium.
C     -PTY1: total mixture property (including PTY6)
C     -PTY2: translational mixture property
C     -PTY3: electronic mixture property
C     -PTY4: rotational mixture property
C     -PTY5: vibrational mixture property
C     -PTY6: entropy of mixing
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      INTEGER LWR1, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), TH, TE, TR, TV(1:NV), P, X(1:NS), 
     &                 PTY1, PTY2, PTY3, PTY4, PTY5, PTY6
C-----------------------------------------------------------------------
      DOUBLE PRECISION SPTY1(1:NS), SPTY2(1:NS), SPTY3(1:NS), 
     &                 SPTY4(1:NS), SPTY5(1:NS), SPTY6(1:NS)
C-----------------------------------------------------------------------
      CALL ENTROPY (WR1, LWR1, WI, LWI, TH, TE, TR, TV, P, SPTY1, 
     &              SPTY2, SPTY3, SPTY4, SPTY5, SPTY6)

      PTY1 = 0.D0; PTY2 = 0.D0; PTY3 = 0.D0
      PTY4 = 0.D0; PTY5 = 0.D0; PTY6 = 0.D0
      DO I = 1, NS
        PTY1 = PTY1 + X(I) *SPTY1(I)
        PTY2 = PTY2 + X(I) *SPTY2(I)
        PTY3 = PTY3 + X(I) *SPTY3(I)
        PTY4 = PTY4 + X(I) *SPTY4(I)
        PTY5 = PTY5 + X(I) *SPTY5(I)
        PTY6 = PTY6 - X(I) *WR1(IUR) *DLOG(X(I))
      ENDDO
      PTY1 = PTY1 +PTY6 

      END SUBROUTINE MIXENTROPY
C---------------------------------------------------------------------- 
