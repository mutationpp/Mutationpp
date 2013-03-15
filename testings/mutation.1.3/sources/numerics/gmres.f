C-----------------------------------------------------------------------
      SUBROUTINE GMRES (N, G, P, X, B, TOLRES, M, I)
C-----------------------------------------------------------------------
C     GMRES with preconditioner P = diag(A).
C     Iterative process stopped when the residual tolerance or maximum 
C     number of iterations is reached. The initial guess is zero.
C-----------------------------------------------------------------------
      IMPLICIT NONE 
      INTEGER N, M , I
      DOUBLE PRECISION G(1:N,1:N), P(1:N), X(1:N), B(1:N), TOLRES
C-----------------------------------------------------------------------
      INTEGER I1, J, K, K1, L
      DOUBLE PRECISION V(1:N,1:M+1), H(1:M+1,1:M), C(1:M), S(1:M), 
     &                 R(1:M+1), Z(1:N), NORM, NORM0, T, GAM, EPS
      PARAMETER (EPS = 1.D-32)
C-----------------------------------------------------------------------
C     Initial guess is zero
      NORM0 = 0.D0
      DO J = 1, N
        X(J)   = 0.D0
        V(J,1) = B(J) 
        NORM0  = NORM0 + V(J,1)*V(J,1)
      ENDDO
      NORM0 = DSQRT(NORM0)
      
C     First term of Hessenberg system:
      IF (NORM0 /= 0.D0) THEN
        DO J = 1, N
          V(J,1) = V(J,1) /NORM0
        ENDDO
      ENDIF
      DO J = 1, M+1
        DO K = 1, M
          H(J,K) = 0.D0
        ENDDO
      ENDDO

C     Iterative loop
      NORM = NORM0; I = 0; R(1) = NORM
      DO WHILE ((I < M) .AND. (NORM > NORM0 *TOLRES))
        I = I+1; I1 = I+1
        
        DO J = 1, N
          Z(J) = V(J,I) /P(J)
        ENDDO
        DO J = 1, N
          V(J,I1) = 0.D0
          DO K = 1, N
           V(J,I1)= V(J,I1)+ G(J,K) *Z(K)
          ENDDO
        ENDDO

C     Modified Graham-Schmidt and Arnoldi step
        DO J = 1, I
          T = 0.D0
          DO K = 1, N
            T = T +V(K,J)*V(K,I1)
          ENDDO
          H(J,I) = T
          DO K = 1, N
            V(K,I1) = V(K,I1) -T *V(K,J)
          ENDDO
        ENDDO

        T = 0.D0
        DO J = 1, N
          T = T +V(J,I1) *V(J,I1)
        ENDDO
        T = DSQRT(T)
        H(I1,I) = T 
        IF ( T /= 0.D0 ) THEN
          DO J = 1, N
            V(J,I1) = V(J,I1) /T
          ENDDO
        ENDIF       

C     Factorization of H:
C     -Perform previous transformation on i-th column of H
        DO K = 2, I
          K1 = K-1
          T = H(K1,I)
          H(K1,I) =  C(K1) *T +S(K1) *H(K,I)
          H(K,I)  = -S(K1) *T +C(K1) *H(K,I)
        ENDDO
        GAM = DSQRT(H(I,I)**2 +H(I1,I)**2)
        IF (GAM == 0.D0) THEN
          GAM = EPS
        ENDIF

C     -Next plane rotation
        C(I)   =  H(I,I)  /GAM
        S(I)   =  H(I1,I) /GAM
        R(I1)  = -S(I) *R(I)
        R(I)   =  C(I) *R(I)

C     -Residual norm 
        H(I,I) =  C(I) *H(I,I) +S(I) *H(I1,I)
        NORM   =  DABS(R(I1))

      ENDDO

C     New solution:
C     -Solve for upper triangular system
        DO K = I, 1, -1
          T = 0.D0
          DO L = K+1, I
            T = T + H(K,L) *R(L)
          ENDDO
          DO L = I+1, K
            T = T + H(K,L) *R(L)
          ENDDO
          R(K) = (R(K) -T) /H(K,K)
        ENDDO

C     -Form linear combination of V(*,I)'s to get solution
        DO J = 1, N
          T = 0.D0
          DO K = 1, I
            T = T + V(J,K) *R(K)
          ENDDO
          Z(J) = T /P(J)
        ENDDO

        DO J = 1, N
          X(J) = X(J) +Z(J)
        ENDDO

      END SUBROUTINE GMRES
C-----------------------------------------------------------------------
