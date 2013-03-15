C-----------------------------------------------------------------------
      SUBROUTINE LAMBDAINT (WR1, LWR1, WR2, LWR2, CP, X, LAMBDA)
C-----------------------------------------------------------------------
C     This subroutine computes the internal thermal conductivity 
C     using Eucken's formula. Inelastic collisions are neglected.
C     The output is the electronic, rotational, vibrational, or total 
C     thermal conductivity according to the specific heat input. 
C     Free electrons are not accounted for.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2 
      DOUBLE PRECISION  WR1(1:LWR1), WR2(1:LWR2), X(1:NS), LAMBDA,
     &                  CP(1:NS)
C-----------------------------------------------------------------------
      INTEGER I
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      LAMBDA = 0.D0 
      DO I = NE+1, NS
        SUM = 0.D0
        DO J = NE+1, I-1 
          IJ = ((J-1)*(2*NS-J)+2*I)/2 
          SUM = SUM + X(J) /WR2(IBINIJ-1+IJ)
        ENDDO
        DO J = I, NS
          IJ = ((I-1)*(2*NS-I)+2*J)/2 
          SUM = SUM + X(J) /WR2(IBINIJ-1+IJ)
        ENDDO
        LAMBDA = LAMBDA + X(I) *CP(I) /WR1(IUNA) /SUM
      ENDDO

      END SUBROUTINE LAMBDAINT
C-----------------------------------------------------------------------
