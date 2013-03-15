C-----------------------------------------------------------------------
      SUBROUTINE EGSDEC (N, A, W, IER)
C-----------------------------------------------------------------------
C     From EGLIB
C     Symmetric positive definite LU decomposition, no pivoting.
C     A is stored in symmetric form
C     A(i,j) --> A(ij) with ij = N*(i-1) - i(i-1)/2 + j for j >= i
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*), W(N)
C-----------------------------------------------------------------------
      IER = 0
      DO  K = 1,N
         KM1 = K - 1
         KP1 = K + 1
         DO  J = 1, KM1
            jj = n*(j-1) - (j*(j-1))/2 + j
            kj = n*(j-1) - (j*(j-1))/2 + k
            W(J)= A(JJ)*A(KJ)
         ENDDO
         kk = n*(k-1) - (k*(k-1))/2 + k
         DO J = 1, KM1
            kj = n*(j-1) - (j*(j-1))/2 + k
            A(KK) = A(KK) - W(J)*A(KJ)
         ENDDO
         IF (A(KK) .EQ. 0.0D0) THEN
            WRITE(6, '(''SINGULAR MATRIX IN EGSDEC'')' )
            IER = K
            RETURN
         ENDIF
         FAC = 1.0D0/A(KK)
         DO J = 1, KM1
            DO I = KP1, N
               ik = n*(k-1) - (k*(k-1))/2 + i
               ij = n*(j-1) - (j*(j-1))/2 + i
               A(IK) = A(IK) - A(IJ)*W(J)
            ENDDO
         ENDDO
         DO I = KP1, N
            ik = n*(k-1) - (k*(k-1))/2 + i
            A(IK) = A(IK)/A(KK)
         ENDDO
      ENDDO
      END SUBROUTINE EGSDEC
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE EGSSOL (N, A, B)
C-----------------------------------------------------------------------
C     Backward substitution ? (EGSLIB)
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*), B(N)
C-----------------------------------------------------------------------
      NM1 = N - 1
      DO J = 1, NM1
         JP1 = J + 1
         FAC = -B(J)
         DO  K = JP1, N
            kj = n*(j-1) - (j*(j-1))/2 + k
            B(K) = B(K) + A(KJ)*FAC
         ENDDO
      ENDDO
      DO J = 1, N
         jj = n*(j-1) - (j*(j-1))/2 + j
         B(J) = B(J)/A(JJ)
      ENDDO
      DO JB = 1, NM1
         J = N + 1 - JB
         JM1 = J - 1
         FAC = -B(J)
         DO K = 1, JM1
            jk = n*(k-1) - (k*(k-1))/2 + j
            B(K) = B(K) + A(JK)*FAC
         ENDDO
      ENDDO
      END SUBROUTINE EGSSOL
C-----------------------------------------------------------------------
