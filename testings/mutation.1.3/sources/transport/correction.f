C-----------------------------------------------------------------------
      SUBROUTINE CORRECTION (WR1, LWR1, WR2, LWR2, X, SONINE, FIJ)
C-----------------------------------------------------------------------
C     This subroutine computes the correction functions of the Stefan-
C     Maxwell equations of the heavy particle gas. 
C     A direct method is used to solve the linear system.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2, SONINE
      DOUBLE PRECISION WR1(1:LWR1), WR2(1:LWR2), X(1:NS),
     &                 FIJ(1:NS*(NS-1)/2)
C-----------------------------------------------------------------------
      DOUBLE PRECISION G(1:(NS-NE+1)*(NS-NE)/2),  
     &                 A(1:(NS-NE+1)*(NS-NE)/2), ALPHA(1:NS-NE),
     &                 TEMP(1:NS-NE), G01(1:NS-NE,1:NS-NE)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      SELECT CASE(SONINE)
        CASE(1)
          DO IS = 1, NS
            DO JS = IS+1, NS
              IJ = ((IS-1)*(2*NS-IS))/2+JS-IS
              FIJ(IJ) = 1.D0
            ENDDO
          ENDDO
        CASE(2)
          CALL CORRECTIONSYST (NS-NE, WR1(IMI+NE), WR2(IBINIJ+NE*NS),
     &                         WR2(IAIJ), WR2(IBIJ), WR2(ICIJ),
     &                         WR2(IETAI), WR1(IUNA), WR1(IUKB),
     &                         X(NE+1), G, G01)
          CALL EGSDEC (NS-NE, G, TEMP, IER)
          DO IS = NE+1, NS
            DO K = 1, NS-NE
              ALPHA(K) = G01(IS-NE,K)
            ENDDO
            DO K = 1, (NS-NE+1)*(NS-NE)/2
              A(K) = G(K)
            ENDDO
            CALL EGSSOL (NS-NE, A, ALPHA)
            DO JS = IS+1, NS
              IJ = ((IS-1)*(2*NS-IS))/2+JS-IS
              FIJ(IJ) = 0.D0
              DO K = 1, NS-NE
                FIJ(IJ) = FIJ(IJ) +G01(JS-NE,K) *ALPHA(K)
              ENDDO          
              FIJ(IJ) = 25.D0 *WR1(IUKB) /(4.D0 *X(IS) *X(JS))
     &                  *WR2(IBINIJ-1+(IS-1)*(2*NS-IS)/2+JS) *FIJ(IJ) 
              FIJ(IJ) = 1.D0 /(1.D0 +FIJ(IJ))
            ENDDO
          ENDDO
        CASE DEFAULT
          WRITE(*,*) 'Lagrange-Sonine order not implemented'
      END SELECT

      END SUBROUTINE CORRECTION 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE CORRECTIONSYST (N, MI, BINIJ, AIJ, BIJ, CIJ, ETAI, NA, 
     &                          KB, X, G, G01)
C-----------------------------------------------------------------------
C     This subroutine computes the system for the correction functions
C     of the Stefan-Maxwell equations of the heavy particle gas. 
C-----------------------------------------------------------------------
      IMPLICIT NONE 
      INTEGER N
      DOUBLE PRECISION MI(1:N), BINIJ(1:(N+1)*N/2), AIJ(1:(N+1)*N/2), 
     &                 BIJ(1:(N+1)*N/2), CIJ(1:(N+1)*N/2), ETAI(1:N), 
     &                 NA, KB, X(1:N), G(1:(N+1)*N/2), G01(1:N,1:N)
C-----------------------------------------------------------------------
      DOUBLE PRECISION FAC1, FAC2, MIIJ, MJIJ
      INTEGER I, J, II, IJ, JJ
C-----------------------------------------------------------------------
      FAC1 = 4.D0 /(15.D0 *KB)
      DO I = 1, N
        II      = ((I-1)*(2*N-I)+2*I)/2 
        G(II)   = X(I) * X(I) /ETAI(I) *MI(I) /NA *FAC1
      ENDDO
      DO I = 1, N
        DO J = I+1, N
          IJ = ((I-1)*(2*N-I)+2*J)/2 
          II = ((I-1)*(2*N-I)+2*I)/2 
          JJ = ((J-1)*(2*N-J)+2*J)/2 
          MIIJ = MI(I) /(MI(I) +MI(J)) 
          MJIJ = MI(J) /(MI(I) +MI(J))
          FAC1 = X(I) *X(J) /(BINIJ(IJ) *25.D0 *KB)
          G(IJ) = -FAC1 *MIIJ *MJIJ *(55.D0 -16.D0 *AIJ(IJ) 
     &            -12.D0 *BIJ(IJ))
          G(II) = G(II) + FAC1 *(30.D0 *MIIJ *MIIJ + 16.D0 *MIIJ *MJIJ 
     &            *AIJ(IJ) +(25.D0 -12.D0 *BIJ(IJ))*MJIJ*MJIJ)
          G(JJ) = G(JJ) + FAC1 *(30.D0 *MJIJ *MJIJ + 16.D0 *MIIJ *MJIJ 
     &            *AIJ(IJ) +(25.D0 -12.D0 *BIJ(IJ))*MIIJ*MIIJ) 
        ENDDO
      ENDDO
      DO I = 1, N
          G01(I,I) = 0.D0
          DO J = 1, I-1
            IJ   = ((J-1)*(2*N-J)+2*I)/2 
            FAC1 = X(I) *X(J) /(BINIJ(IJ) *25.D0 *KB) /(MI(I) +MI(J))
            FAC1 = FAC1 *(12.D0 *CIJ(IJ) -10.D0)
            G01(I,J) = FAC1 *MI(I)
            G01(I,I) = G01(I,I) -FAC1 *MI(J) 
          ENDDO
          DO J = I+1, N
            IJ = ((I-1)*(2*N-I)+2*J)/2 
            FAC1 = X(I) *X(J) /(BINIJ(IJ) *25.D0 *KB) /(MI(I) +MI(J))
            FAC1 = FAC1 *(12.D0 *CIJ(IJ) -10.D0)
            G01(I,J) = FAC1 *MI(I) 
            G01(I,I) = G01(I,I) -FAC1 *MI(J) 
          ENDDO
      ENDDO

      END SUBROUTINE CORRECTIONSYST
C-----------------------------------------------------------------------
