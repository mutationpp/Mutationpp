C-----------------------------------------------------------------------

C     Ionized mixture, system for E and V_i

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SMD (WR1, LWR1, WR2, LWR2, WI, LWI, X, TH, TE, ND, DF,
     &                FIJ, JDIF, EAMB)
C-----------------------------------------------------------------------
C     This subroutine computes the mass diffusion fluxes solution of the
C     Stefan-Maxwell equations in the case of a mixture with charged 
C     particles (ambipolar constraint) in thermal non-equilibrium. If 
C     considered, thermal diffusion is included in the driving forces.
C     A Gaussian elimination is used to solve for the non-symmetric 
C     indefinite linear system. The mass constraint is incorporated 
C     into the linear system.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), WR2(1:LWR2), X(1:NS), TH, TE, ND, 
     &                 DF(1:NS), FIJ(1:NS*(NS-1)/2), JDIF(1:NS), EAMB
C-----------------------------------------------------------------------
      DOUBLE PRECISION G(1:NS+1,1:NS+1), BETA(1:NS+1), INDX(1:NS+1), MM,
     &                 D
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
C     Molar mass
      MM = 0.D0
      DO I = 1, NS
        MM = MM +X(I) *WR1(IMI+I-1)
      ENDDO

C     System for V_i
      CALL SYSTSMD (NS, WR1(IMI), WR2(IBINIJ), WR1(IUE), WI(IQI), X, ND,
     &             DF, FIJ, MM, WR1(IUKB), TH, TE, G, BETA, SF)

C     Direct solution
      CALL LUDCMP(G, NS+1, NS+1, INDX, D)
      CALL LUBKSB(G, NS+1, NS+1, INDX, BETA)         

C     Mass diffusion fluxes
      DO I = 1, NS
       JDIF(I) = X(I) *BETA(I) *WR1(IMI+I-1) /WR1(IUNA) *ND
      ENDDO

C     Ambipolar field (scaling factor correction)
      EAMB = BETA(NS+1) /SF
                 
      END SUBROUTINE SMD
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SYSTSMD (NS, MI, BINIJ, E, QI, X, ND, DF, FIJ, MM, KB,
     &                   TH, TE, G, BETA, SF)
C-----------------------------------------------------------------------
C     This subroutine computes the system for the mass diffusion fluxes
C     using the Stefan-Maxwell equations in the case of a mixture 
C     with charged particles (ambipolar field constraint) in thermal 
C     non-equilibrium. 
C     If considered, thermal diffusion is included in the driving 
C     forces. The mass constraint is incorporated into the non-
C     symmetric linear system.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NS, QI(1:NS)
      DOUBLE PRECISION MI(1:NS), BINIJ(1:(NS+1)*NS/2), E, X(1:NS), ND,
     &                 DF(1:NS), FIJ(1:NS*(NS-1)/2), MM, KB, TH, TE,
     &                 G(1:NS+1,1:NS+1), BETA(1:NS+1), SF
C-----------------------------------------------------------------------
      INTEGER I, J, IJ, II, JJ
      DOUBLE PRECISION QT, DIJAV, M(1:NS), DIJ(1:NS,1:NS), KA(1:NS)
C-----------------------------------------------------------------------
      QT = 0.D0
      DO I = 1, NS
        BETA(I) =  -DF(I) 
        M(I)    = X(I) *MI(I) /MM
        QT = QT +X(I) *QI(I)
      ENDDO
      BETA(1) = BETA(1) *TH /TE 
      BETA(NS+1) = 0.D0

      SF = 0.D0
      DIJAV = 0.D0
      DO I = 1, NS
        KA(I) =  X(I) *(QI(I) -MI(I) /MM *QT) *E /(KB *TH)
        SF    =  SF +KA(I) *KA(I)
        DIJ(I,I)  = 1.D0
        G(I,I)    = 0.D0
        DO J = I+1, NS
          IJ = ((I-1)*(2*NS-I)+2*J)/2
          DIJ(I,J)   = BINIJ(IJ) /ND *FIJ(((I-1)*(2*NS-I))/2+J-I) 
          DIJ(J,I)   = DIJ(I,J)
          G(I,J)     = 0.D0
          G(J,I)     = 0.D0
          DIJAV      = DIJAV +DIJ(I,J) 
        ENDDO
      ENDDO             
      DIJAV = DIJAV /(NS*(NS-1)/2)
      SF = DSQRT(SF)

C     Electron contribution
      I = 1
        DO J = I+1, NS 
          G(I,J) = -X(I) *X(J) /DIJ(I,J)
          G(I,I) = G(I,I) -G(I,J)
          G(I,J) = G(I,J) *TE /TH
          G(J,I) = G(I,J)
          G(J,J) = G(J,J) -G(I,J) *TE /TH
        ENDDO 
        J = NS+1
          G(I,J) = -KA(I) *TH / (TE *SF)
Cthierry
C          G(I,J) = 0. 
Cthierry

C     Heavy particle contribution
      DO I = 2, NS
        DO J = I+1, NS
          G(I,J) = -X(I) *X(J) /DIJ(I,J)
          G(J,I) = G(I,J) 
          G(I,I) = G(I,I) -G(I,J)
          G(J,J) = G(J,J) -G(I,J)
        ENDDO
        J = NS+1
          G(I,J) = -KA(I) /SF
Cthierry
C          G(I,J) = 0. 
Cthierry
      ENDDO    

C     Mass constraint
      DO I = 1, NS
        DO J = 1, NS
          G(I,J) = G(I,J) +M(I) *M(J) /DIJAV
        ENDDO
      ENDDO    

C     Ambipolar constraint:
      I = NS+1
      DO J = 1, NS
          G(I,J) = -KA(J) /SF
Cthierry
C          G(I,J) = 0. 
Cthierry
      ENDDO    
      G(NS+1,NS+1) = 0.D0
Cthierry
C      G(NS+1,NS+1) = 1.D0
Cthierry

      END SUBROUTINE SYSTSMD
C----------------------------------------------------------------------- 
C-----------------------------------------------------------------------
      SUBROUTINE SMGMRES (WR1, LWR1, WR2, LWR2, WI, LWI, X, TH, TE, ND,
     &                    DF, FIJ, JDIF, EAMB)
C-----------------------------------------------------------------------
C     This subroutine computes the mass diffusion fluxes solution of the
C     Stefan-Maxwell equations in the case of a mixture with charged 
C     particles (ambipolar constraint) in thermal non-equilibrium. If 
C     considered, thermal diffusion is included in the driving forces.
C     The GMRES is used to solve for the non-symmetric indefinite linear
C     system. The mass constraint is incorporated into the linear 
C     system.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), WR2(1:LWR2), X(1:NS), TH, TE, ND, 
     &                 DF(1:NS), FIJ(1:NS*(NS-1)/2), JDIF(1:NS), EAMB
C-----------------------------------------------------------------------
      DOUBLE PRECISION G(1:NS+1,1:NS+1), BETA(1:NS+1), INDX(1:NS+1), MM,
     &                 D, M(1:NS), R(1:NS), MX, P(1:NS+1), 
     &                 ALPHA(1:NS+1), TOLRES
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
C     Molar mass
      MM = 0.D0
      DO I = 1, NS
        MM = MM +X(I) *WR1(IMI+I-1)
      ENDDO

C     System for V_i
      CALL SYSTSMK (NS, WR1(IMI), WR2(IBINIJ), WR1(IUE), WI(IQI), X, ND,
     &             DF, FIJ, MM, WR1(IUKB), TH, TE, G, BETA, SF, M, R)

C     GMRES iterative solution
      DO I = 1, NS
        P(I)= G(I,I)
      ENDDO
      P(NS+1) = 1.D0

      ITERMAX = 4; TOLRES = 1.D-2
      CALL GMRES (NS+1, G, P, ALPHA, BETA, TOLRES, ITERMAX, ITER)

C     Projection
      RM = 0.D0; MX = 0.D0
      DO I = 1, NS
        RM = RM +M(I) *R(I)
        MX = MX +M(I) *ALPHA(I)
      ENDDO
      PROJ = MX /RM
      DO I = 1, NS
       JDIF(I) = ALPHA(I) -PROJ *R(I)
      ENDDO

C     Mass diffusion fluxes
      DO I = 1, NS
       JDIF(I) = X(I) *JDIF(I) *WR1(IMI+I-1) /WR1(IUNA) *ND
      ENDDO

C     Ambipolar field (scaling factor correction)
      EAMB = ALPHA(NS+1) /SF
                 
      END SUBROUTINE SMGMRES
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SYSTSMK (NS, MI, BINIJ, E, QI, X, ND, DF, FIJ, MM, KB,
     &                    TH, TE, G, BETA, SF, M, R)
C-----------------------------------------------------------------------
C     This subroutine computes the system for the mass diffusion fluxes
C     using the Stefan-Maxwell equations in the case of a mixture 
C     with charged particles (ambipolar field constraint) in thermal 
C     non-equilibrium. 
C     If considered, thermal diffusion is included in the driving 
C     forces. The projector onto the mass constraint hyperplane along 
C     the nullspace of G is computed. 
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NS, QI(1:NS)
      DOUBLE PRECISION MI(1:NS), BINIJ(1:(NS+1)*NS/2), E, X(1:NS), ND,
     &                 DF(1:NS), FIJ(1:NS*(NS-1)/2), MM, KB, TH, TE,
     &                 G(1:NS+1,1:NS+1), BETA(1:NS+1), SF, M(1:NS), 
     &                 R(1:NS)
C-----------------------------------------------------------------------
      INTEGER I, J, IJ, II, JJ
      DOUBLE PRECISION QT, DIJ(1:NS,1:NS), KA(1:NS)
C-----------------------------------------------------------------------
      QT = 0.D0
      DO I = 1, NS
        BETA(I) =  -DF(I) 
        M(I)    = X(I) *MI(I) /MM
        R(I)    = 1.D0
        QT      = QT +X(I) *QI(I)
      ENDDO
      BETA(1)    = BETA(1) *TH /TE 
      R(1)       = R(1) *TH /TE 
      BETA(NS+1) = 0.D0

      SF = 0.D0
      DO I = 1, NS
        KA(I)     =  X(I) *(QI(I) -MI(I) /MM *QT) *E /(KB *TH)
        SF        =  SF +KA(I) *KA(I)
        DIJ(I,I)  = 1.D0
        G(I,I)    = 0.D0
        DO J = I+1, NS
          IJ = ((I-1)*(2*NS-I)+2*J)/2
          DIJ(I,J)   = BINIJ(IJ) /ND *FIJ(((I-1)*(2*NS-I))/2+J-I) 
          DIJ(J,I)   = DIJ(I,J)
          G(I,J)     = 0.D0
          G(J,I)     = 0.D0
        ENDDO
      ENDDO             
      SF = DSQRT(SF)

C     Electron contribution
      I = 1
        DO J = I+1, NS 
          G(I,J) = -X(I) *X(J) /DIJ(I,J)
          G(I,I) = G(I,I) -G(I,J)
          G(I,J) = G(I,J) *TE /TH
          G(J,I) = G(I,J)
          G(J,J) = G(J,J) -G(I,J) *TE /TH
        ENDDO 
        J = NS+1
          G(I,J) = -KA(I) *TH / (TE *SF)

C     Heavy particle contribution
      DO I = 2, NS
        DO J = I+1, NS
          G(I,J) = -X(I) *X(J) /DIJ(I,J)
          G(J,I) = G(I,J) 
          G(I,I) = G(I,I) -G(I,J)
          G(J,J) = G(J,J) -G(I,J)
        ENDDO
        J = NS+1
          G(I,J) = -KA(I) /SF
      ENDDO    

C     Ambipolar constraint:
      I = NS+1
      DO J = 1, NS
          G(I,J) = -KA(J) /SF
      ENDDO    
      G(NS+1,NS+1) = 0.D0

      END SUBROUTINE SYSTSMK
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SMDSPD (WR1, LWR1, WR2, LWR2, WI, LWI, X, TH, TE, ND,
     &                   DF, FIJ, JDIF, EAMB)
C-----------------------------------------------------------------------
C     This subroutine computes the mass diffusion fluxes solution of the
C     Stefan-Maxwell equations in the case of a mixture with charged 
C     particles (ambipolar constraint) in thermal non-equilibrium. If 
C     considered, thermal diffusion is included in the driving forces.
C     A variant of Choleki's method is used twice to solve for the
C     linear systems for the driving forces and the ambipolar field.
c     The ambipolar field is deduced from a projection onto the 
C     ambipolar constraint hyperplane.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), WR2(1:LWR2), X(1:NS), TH, TE, ND,
     &                 DF(1:NS), FIJ(1:NS*(NS-1)/2), JDIF(1:NS), EAMB
C-----------------------------------------------------------------------
      DOUBLE PRECISION G(1:(NS+1)*NS/2), BETA1(1:NS), BETA2(1:NS), MM, 
     &                 KA(1:NS), NUM, DEN, EPRIME, TEMP(1:NS)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
C     Molar mass
      MM = 0.D0
      DO I = 1, NS
        MM = MM +X(I) *WR1(IMI+I-1)
      ENDDO

C     System
      CALL SYSTSMDSPD (NS, WR1(IMI), WR2(IBINIJ), WI(IQI), X, ND, DF, 
     &                 FIJ, G, BETA1, BETA2, KA, MM, WR1(IUKB), TH, TE, 
     &                 SF, WR1(IUE))

C     LU decomposition
      CALL EGSDEC (NS, G, TEMP, IER)

C     Driving forces contribution
      CALL EGSSOL (NS, G, BETA1)

C     Ambipolar field contribution
      CALL EGSSOL (NS, G, BETA2)

C     Ambipolar projection
      NUM = 0.D0; DEN = 0.D0
      DO I = 1, NS
        NUM = NUM + KA(I) *BETA1(I)
        DEN = DEN + KA(I) *BETA2(I)
      ENDDO
      EPRIME = - NUM /DEN

C     Mass diffusion fluxes
      DO I = 1, NS
       JDIF(I) = (BETA1(I) +EPRIME *BETA2(I)) *X(I) *WR1(IMI+I-1)
     &           /WR1(IUNA) *ND
      ENDDO

      EAMB = EPRIME /SF

      END SUBROUTINE SMDSPD
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SYSTSMDSPD (N, MI, BINIJ, QI, X, ND, DF, FIJ, G, BETA1,
     &                       BETA2, KA, MM, KB, TH, TE, SF, E)
C-----------------------------------------------------------------------
C     This subroutine computes the system for the mass diffusion fluxes
C     using the Stefan-Maxwell equations in the case of a mixture
C     with charged particles (ambipolar field), for the DSP method. If 
C     considered, thermal diffusion is included in the driving forces.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N, QI(1:N)
      DOUBLE PRECISION MI(1:N), BINIJ(1:(N+1)*N/2), X(1:N), ND,
     &                 DF(1:N), FIJ(1:N*(N-1)/2), G(1:(N+1)*N/2), 
     &                 BETA1(1:N), BETA2(1:N), KA(1:N), MM, KB, TH, TE,
     &                 SF, E
C-----------------------------------------------------------------------
      INTEGER I, J, IJ, II, JJ
      DOUBLE PRECISION FAC, QT, DIJAV, M(1:N)
C-----------------------------------------------------------------------
      QT = 0.D0
      DO I = 1, N
        BETA1(I) = -DF(I)
        M(I)     = X(I) *MI(I) /MM
        QT       = QT +X(I) *QI(I)
        II      = ((I-1)*(2*N-I)+2*I)/2
        G(II)   = 0.D0
      ENDDO
      BETA1(1) = BETA1(1) *TH /TE

      SF = 0.D0
      DIJAV = 0.D0
      DO I = 1, N
        KA(I)    =  X(I) *(QI(I) -MI(I) /MM *QT) *E /(KB *TH)
        SF    =  SF +KA(I) *KA(I)
        DO J = I+1, N
          IJ = ((I-1)*(2*N-I)+2*J)/2
          DIJAV = DIJAV +BINIJ(IJ) /ND
        ENDDO
      ENDDO
      DIJAV = DIJAV /(N*(N-1)/2)
      SF = DSQRT(SF)
      DO I = 1, N
        KA(I)    = KA(I) /SF
        BETA2(I) = KA(I) 
      ENDDO
      BETA2(1) = BETA2(1) *TH /TE

      I = 1
        DO J = I+1, N
          IJ = ((I-1)*(2*N-I)+2*J)/2
          II = ((I-1)*(2*N-I)+2*I)/2
          JJ = ((J-1)*(2*N-J)+2*J)/2
          FAC   = X(I)*X(J) *ND /(BINIJ(IJ) *FIJ(((I-1)*(2*N-I))/2+J-I))
          G(IJ) = -FAC *TE /TH
          G(II) = G(II) + FAC
          G(JJ) = G(JJ) - G(IJ) *TE /TH
        ENDDO
      DO I = 2, N
        DO J = I+1, N
          IJ = ((I-1)*(2*N-I)+2*J)/2
          II = ((I-1)*(2*N-I)+2*I)/2
          JJ = ((J-1)*(2*N-J)+2*J)/2
          FAC   = X(I)*X(J) *ND /(BINIJ(IJ) *FIJ(((I-1)*(2*N-I))/2+J-I))
          G(IJ) = -FAC
          G(II) = G(II) + FAC
          G(JJ) = G(JJ) + FAC
        ENDDO
      ENDDO

C     Mass constraint
      DO I = 1, N
        DO J = I, N
          IJ = ((I-1)*(2*N-I)+2*J)/2
          G(IJ) = G(IJ) +M(I) *M(J) /DIJAV
        ENDDO
      ENDDO

      END SUBROUTINE SYSTSMDSPD
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SMDSPCG (WR1, LWR1, WR2, LWR2, WI, LWI, X, TH, TE, ND,
     &                    DF, FIJ, JDIF, EAMB)
C-----------------------------------------------------------------------
C     This subroutine computes the mass diffusion fluxes solution of the
C     Stefan-Maxwell equations in the case of a mixture with charged 
C     particles (ambipolar constraint) in thermal non-equilibrium. If 
C     considered, thermal diffusion is included in the driving forces.
C     The CG or MINRES is used twice to solve for the linear systems 
C     for the driving forces and the ambipolar field. The ambipolar 
C     field is deduced from a projection onto the ambipolar constraint 
C     hyperplane.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), WR2(1:LWR2), X(1:NS), TH, TE, ND,
     &                 DF(1:NS), FIJ(1:NS*(NS-1)/2), JDIF(1:NS), EAMB
C-----------------------------------------------------------------------
      DOUBLE PRECISION G(1:(NS+1)*NS/2), BETA1(1:NS), BETA2(1:NS), MM, 
     &                 KA(1:NS), NUM, DEN, EPRIME, TOLRES, R(1:NS),
     &                 M(1:NS), V1(1:NS), V2(1:NS), MX, P(1:NS)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
C     Molar mass
      MM = 0.D0
      DO I = 1, NS
        MM = MM +X(I) *WR1(IMI+I-1)
      ENDDO

C     System
      CALL SYSTSMDSPK (NS, WR1(IMI), WR2(IBINIJ), WI(IQI), X, ND, DF, 
     &                 FIJ, G, R, M, BETA1, BETA2, KA, MM, WR1(IUKB),
     &                 TH, TE, SF, WR1(IUE))

C     Krylov iterative solution
      DO I = 1, NS
        II = ((I-1)*(2*NS-I)+2*I)/2
        P(I) = G(II)
      ENDDO
      ITERMAX = 3; TOLRES1 = 6.D-3; TOLRES2 = 1.D-2


C     Driving forces contribution
      CALL CG (NS, G, P, V1, BETA1, TOLRES1, ITERMAX, ITER)
      RM = 0.D0; MX = 0.D0
      DO I = 1, NS
        RM = RM +M(I) *R(I)
        MX = MX +M(I) *V1(I)
      ENDDO
      PROJ = MX /RM
      DO I = 1, NS
       V1(I) = V1(I) -PROJ *R(I)
      ENDDO

C     Ambipolar field contribution
      CALL CG (NS, G, P, V2, BETA2, TOLRES2, ITERMAX, ITER)
      MX = 0.D0
      DO I = 1, NS
        MX = MX +M(I) *V2(I)
      ENDDO
      PROJ = MX /RM
      DO I = 1, NS
       V2(I) = V2(I) -PROJ *R(I)
      ENDDO

C     Ambipolar projection
      NUM = 0.D0; DEN = 0.D0
      DO I = 1, NS
        NUM = NUM + KA(I) *V1(I)
        DEN = DEN + KA(I) *V2(I)
      ENDDO
      EPRIME = - NUM /DEN

C     Mass diffusion fluxes
      DO I = 1, NS
       JDIF(I) = (V1(I) +EPRIME *V2(I)) *X(I) *WR1(IMI+I-1)
     &           /WR1(IUNA) *ND
      ENDDO

      EAMB = EPRIME /SF

      END SUBROUTINE SMDSPCG
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SYSTSMDSPK (N, MI, BINIJ, QI, X, ND, DF, FIJ, G, R, M,
     &                       BETA1, BETA2, KA, MM, KB, TH, TE, SF, E)
C-----------------------------------------------------------------------
C     This subroutine computes the system for the mass diffusion fluxes
C     using the Stefan-Maxwell equations in the case of a mixture
C     with charged particles (ambipolar field), for the DSP method. If 
C     considered, thermal diffusion is included in the driving forces.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N, QI(1:N)
      DOUBLE PRECISION MI(1:N), BINIJ(1:(N+1)*N/2), X(1:N), ND,
     &                 DF(1:N), FIJ(1:N*(N-1)/2), G(1:(N+1)*N/2), 
     &                 R(1:N), M(1:N), BETA1(1:N), BETA2(1:N), KA(1:N),
     &                 MM, KB, TH, TE, SF, E
C-----------------------------------------------------------------------
      INTEGER I, J, IJ, II, JJ
      DOUBLE PRECISION FAC, QT
C-----------------------------------------------------------------------
      QT = 0.D0
      DO I = 1, N
        BETA1(I) = -DF(I)
        M(I)     = X(I) *MI(I) /MM
        R(I)     = 1.D0
        QT       = QT +X(I) *QI(I)
        II       = ((I-1)*(2*N-I)+2*I)/2
        G(II)    = 0.D0
      ENDDO
      BETA1(1) = BETA1(1) *TH /TE
      R(1)     = R(1) *TH /TE

      SF = 0.D0
      DO I = 1, N
        KA(I)    =  X(I) *(QI(I) -MI(I) /MM *QT) *E /(KB *TH)
        SF    =  SF +KA(I) *KA(I)
        DO J = I+1, N
          IJ = ((I-1)*(2*N-I)+2*J)/2
        ENDDO
      ENDDO
      SF = DSQRT(SF) 
      DO I = 1, N
        KA(I)    = KA(I) /SF
        BETA2(I) = KA(I) 
      ENDDO
      BETA2(1) = BETA2(1) *TH /TE

      I = 1
        DO J = I+1, N
          IJ = ((I-1)*(2*N-I)+2*J)/2
          II = ((I-1)*(2*N-I)+2*I)/2
          JJ = ((J-1)*(2*N-J)+2*J)/2
          FAC   = X(I)*X(J) *ND /(BINIJ(IJ) *FIJ(((I-1)*(2*N-I))/2+J-I))
          G(IJ) = -FAC *TE /TH
          G(II) = G(II) + FAC
          G(JJ) = G(JJ) - G(IJ) *TE /TH
        ENDDO
      DO I = 2, N
        DO J = I+1, N
          IJ = ((I-1)*(2*N-I)+2*J)/2
          II = ((I-1)*(2*N-I)+2*I)/2
          JJ = ((J-1)*(2*N-J)+2*J)/2
          FAC   = X(I)*X(J) *ND /(BINIJ(IJ) *FIJ(((I-1)*(2*N-I))/2+J-I))
          G(IJ) = -FAC
          G(II) = G(II) + FAC
          G(JJ) = G(JJ) + FAC
        ENDDO
      ENDDO

      END SUBROUTINE SYSTSMDSPK
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SMSUT (WR1, LWR1, WR2, LWR2, WI, LWI, X, TH, TE, ND,
     &                  DF, FIJ, JDIF, EAMB)
C-----------------------------------------------------------------------
C     This subroutine computes the mass diffusion fluxes solution of the
C     Stefan-Maxwell equations in the case of a mixture with charged
C     particles (ambipolar constraint) in thermal non-equilibrium. If
C     considered, thermal diffusion is included in the driving forces.
C     The Sutton & Gnoffo algorithm is used. The generalization of this
C     algorithm to ionized mixtures is due to David Vanden Abeele and 
C     Paolo Barbante. 
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), WR2(1:LWR2), X(1:NS), TH, TE, ND,
     &                 DF(1:NS), FIJ(1:NS*(NS-1)/2), JDIF(1:NS), EAMB
C-----------------------------------------------------------------------
      DOUBLE PRECISION D(1:NS,1:NS), DEF(1:NS), NOLD(1:NS), NNEW(1:NS),
     &                 RES(1:NS), MM, TOLRES, Q
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
C     Binary diffusion coefficient. 
      DO I = 1, NS
        II     = ((I-1)*(2*NS-I)+2*I)/2
        D(I,I) = WR2(IBINIJ+II-1) /ND
        DO J = I+1, NS
          IJ     = ((I-1)*(2*NS-I)+2*J)/2
          D(I,J) = WR2(IBINIJ+IJ-1) /ND *FIJ(((I-1)*(2*NS-I))/2+J-I)
          D(J,I) = D(I,J)
        ENDDO
      ENDDO

C     Mixture molar mass
      MM = 0.D0
      DO I = 1, NS
        MM = MM +X(I) *WR1(IMI+I-1)
      ENDDO

C     Effective binary diffusion coefficient
      I = 1
        SUM = 0.D0
        DO J = 1, NS
          SUM = SUM +X(J) /D(I,J)
        ENDDO
        DEF(I) = 1.D0 /SUM
      DO I = 2, NS
        SUM = 0.D0
        J = 1
          SUM = SUM +X(J) /D(I,J) *(TE /TH)**2
        DO J = 2, NS 
          SUM = SUM +X(J) /D(I,J)
        ENDDO
        DEF(I) = 1.D0 /SUM
      ENDDO

C     Initial conditions, number density diffusion fluxes set to zero
      RES2INI = 0.D0
      I = 1
        NNEW(I) = 0.D0
        NOLD(I) = NNEW(I)
        RES(I)  = -ND *DEF(I) *DF(I) *TH /TE
        RES2INI = RES2INI +RES(I) *RES(I)
      DO I = 2, NS
        NNEW(I) = 0.D0
        NOLD(I) = NNEW(I)
        RES(I)  = -ND *DEF(I) *DF(I)
        RES2INI = RES2INI +RES(I) *RES(I)
      ENDDO
      RES2 = 1.D0

C     -Iterative process, number density diffusion fluxes
      TOLRES = 1.D-3; ITERMAX = 8 
      TOLRES2 = TOLRES *TOLRES
      ITER = 0
      DO WHILE ( (ITER< ITERMAX) .AND. (RES2 > TOLRES2) )
        ITER = ITER +1

C     -Driving force contribution to residual
        I = 1
          SUM = 0.D0
          DO J = 2, NS
            SUM = SUM +NOLD(J) /D(I,J)
          ENDDO
          SUM = SUM *TE /TH
          J = 1
            SUM = SUM +NOLD(J) /D(I,J)
          RES(I)  = -NOLD(I) -ND *DEF(I) *DF(I) *TH /TE 
     &              +X(I) *DEF(I) *SUM
        DO I = 2, NS
          SUM = 0.D0
          J = 1
            SUM = SUM +NOLD(J) /D(I,J) *TE /TH
          DO J = 2, NS 
            SUM = SUM +NOLD(J) /D(I,J)
          ENDDO
          RES(I)  = -NOLD(I) -ND *DEF(I) *DF(I) +X(I) *DEF(I) *SUM
        ENDDO

C     -Ambipolar driving force
        SUM1 = 0.D0; SUM2 = 0.D0
        DO I = 1, NS
          Q    = DBLE(WI(IQI+I-1))
          SUM1 = SUM1 + Q *(NOLD(I) +RES(I))
          SUM2 = SUM2 + X(I) *DEF(I) *Q *Q
        ENDDO
        SUM2 = SUM2 *ND
        DAMB = -SUM1 /SUM2

C     -Ambipolar driving force contribution to residual
        RES2 = 0.D0
        DO I = 1, NS
          RES(I) = RES(I) + DAMB * ND *X(I) *DEF(I) *DBLE(WI(IQI+I-1))
          RES2 = RES2 +RES(I) *RES(I)
          NNEW(I) = NOLD(I) +RES(I)
        ENDDO
        RES2 = RES2 /RES2INI

C     -Correction for mass constraint
        DO I = 1, NS
          SUM = 0.D0
          DO J = 1, NS
            SUM = SUM +WR1(IMI+J-1) *NNEW(J)
          ENDDO
          NNEW(I) = NNEW(I) -X(I) /MM *SUM
        ENDDO

C     -Update
        DO I = 1, NS
          NOLD(I) = NNEW(I)
        ENDDO

      ENDDO
C     Mass diffusion fluxes
      DO I = 1, NS
        JDIF(I) = WR1(IMI+I-1) /WR1(IUNA) *NNEW(I)
      ENDDO

C     Ambipolar field
      EAMB = DAMB *WR1(IUKB) *TH /WR1(IUE)

      END SUBROUTINE SMSUT
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SMKOL (WR1, LWR1, WR2, LWR2, WI, LWI, X, TH, TE, ND,
     &                  DF, FIJ, JDIF, EAMB)
C-----------------------------------------------------------------------
C     This subroutine computes the mass diffusion fluxes solution of the
C     Kolesnikov equations in the case of a mixture with charged
C     particles (ambipolar constraint) in thermal non-equilibrium. If
C     considered, thermal diffusion is included in the driving forces.
C     A direct method is used to solve for the linear system. The
C     mass constraint is incorporated into the non-symmetric linear
C     system.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), WR2(1:LWR2), X(1:NS), TH, TE, ND, 
     &                 DF(1:NS), FIJ(1:NS*(NS-1)/2), JDIF(1:NS), EAMB
C-----------------------------------------------------------------------
      DOUBLE PRECISION G(1:NS+1,1:NS+1), BETA(1:NS+1), INDX(1:NS+1), MM,
     &                 D
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
C     Molar mass
      MM = 0.D0
      DO I = 1, NS
        MM = MM +X(I) *WR1(IMI+I-1)
      ENDDO

C     System for V_i
      CALL SYSTSMKOL (NS, WR1(IMI), WR2(IBINIJ), WR1(IUE), WI(IQI), X,
     &                ND, DF, FIJ, MM, WR1(IUKB), TH, TE, G, BETA, SF)

C     Direct solution
      CALL LUDCMP(G, NS+1, NS+1, INDX, D)
      CALL LUBKSB(G, NS+1, NS+1, INDX, BETA)         

C     Mass diffusion fluxes
      DO I = 1, NS
       JDIF(I) = X(I) *BETA(I) *WR1(IMI+I-1) /WR1(IUNA) *ND
      ENDDO

C     Ambipolar field (scaling factor correction)
      EAMB = BETA(NS+1) /SF
                 
      END SUBROUTINE SMKOL
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SYSTSMKOL (NS, MI, BINIJ, E, QI, X, ND, DF, FIJ, MM, 
     &                      KB, TH, TE, G, BETA, SF)
C-----------------------------------------------------------------------
C     This subroutine computes the system for the mass diffusion fluxes
C     using the Kolesnikov equations in the case of a mixture 
C     with charged particles (ambipolar field constraint) in thermal 
C     non-equilibrium. 
C     If considered, thermal diffusion is included in the driving 
C     forces. The mass constraint is incorporated into the non-
C     symmetric linear system.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NS, QI(1:NS)
      DOUBLE PRECISION MI(1:NS), BINIJ(1:(NS+1)*NS/2), E, X(1:NS), ND,
     &                 DF(1:NS), FIJ(1:NS*(NS-1)/2), MM, KB, TH, TE,
     &                 G(1:NS+1,1:NS+1), BETA(1:NS+1), SF
C-----------------------------------------------------------------------
      INTEGER I, J, IJ, II, JJ
      DOUBLE PRECISION QT, DIJAV, M(1:NS), DIJ(1:NS,1:NS), KA(1:NS)
C-----------------------------------------------------------------------
      QT = 0.D0
      DO I = 1, NS
        BETA(I) =  -DF(I) 
        M(I)    = X(I) *MI(I) /MM
        QT = QT +X(I) *QI(I)
      ENDDO
      BETA(1) = BETA(1) *TH /TE 
      BETA(NS+1) = 0.D0

      SF = 0.D0
      DIJAV = 0.D0
      DO I = 1, NS
        KA(I) =  X(I) *(QI(I) -MI(I) /MM *QT) *E /(KB *TH)
        SF    =  SF +KA(I) *KA(I)
        DIJ(I,I)  = 1.D0
        G(I,I)    = 0.D0
        DO J = I+1, NS
          IJ = ((I-1)*(2*NS-I)+2*J)/2
          DIJ(I,J)   = BINIJ(IJ) /ND *FIJ(((I-1)*(2*NS-I))/2+J-I) 
          DIJ(J,I)   = DIJ(I,J)
          G(I,J)     = 0.D0
          G(J,I)     = 0.D0
          DIJAV      = DIJAV +DIJ(I,J) 
        ENDDO
      ENDDO             
      DIJAV = DIJAV /(NS*(NS-1)/2)
      SF = DSQRT(SF) 

C     Electron contribution
      I = 1
        DO J = I+1, NS 
          G(I,J) = 0.D0
          G(J,I) = -X(I) *X(J) /DIJ(I,J)
          G(I,I) = G(I,I) -G(J,I)
          G(J,I) = G(J,I) *TE /TH
        ENDDO 
        J = NS+1
          G(I,J) = -KA(I) *TH / (TE *SF)

C     Heavy particle contribution
      DO I = 2, NS
        DO J = I+1, NS
          G(I,J) = -X(I) *X(J) /DIJ(I,J)
          G(J,I) = G(I,J) 
          G(I,I) = G(I,I) -G(I,J)
          G(J,J) = G(J,J) -G(I,J)
        ENDDO
        J = NS+1
          G(I,J) = -KA(I) /SF
      ENDDO    

C     Mass constraint
      DO I = 1, NS
        DO J = 1, NS
          G(I,J) = G(I,J) +M(I) *M(J) /DIJAV
        ENDDO
      ENDDO    

C     Ambipolar constraint:
      I = NS+1
      DO J = 1, NS
          G(I,J) = -KA(J) /SF
      ENDDO    
      G(NS+1,NS+1) = 0.D0

      END SUBROUTINE SYSTSMKOL
C----------------------------------------------------------------------- 
C-----------------------------------------------------------------------
      SUBROUTINE SMFICK (WR1, LWR1, WR2, LWR2, WI, LWI, X, TH, TE, ND,
     &                   DF, FIJ, JDIF, JDIFR, EAMB)
C-----------------------------------------------------------------------
C     This subroutine computes the mass diffusion fluxes for Fick's law
C     with Ramshaw's correction.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), WR2(1:LWR2), X(1:NS), TH, TE, ND,
     &                 DF(1:NS), FIJ(1:NS*(NS-1)/2), JDIF(1:NS), EAMB,
     &                 JDIFR(1:NS)
C-----------------------------------------------------------------------
      INTEGER I, J
      DOUBLE PRECISION G(1:(NS-1)*NS/2), BETA(1:NS-1), MM, TEMP(1:NS-1),
     &                 NDIMI(1:NS), Y(1:NS), JPROJ, QT, KI(1:NS)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
C     Molar mass
      MM = 0.D0
      DO I = 1, NS
        MM = MM +X(I) *WR1(IMI+I-1)
      ENDDO

C     System
      CALL SYSTSMFICK (NS, WR1(IMI), WR2(IBINIJ), MM, X, Y, NDIMI) 

C     Ambipolar field
      IF (NE/=0) THEN      
        QT = 0.D0
        DO I = 1, NS
          QT = QT +X(I) *WI(IQI+I-1)
        ENDDO
        DO I = 1, NS
          KI(I) =  X(I) *(WI(IQI+I-1) -WR1(IMI+I-1) /MM *QT) *WR1(IUE) 
     &        /(WR1(IUKB) *TH)
        ENDDO
        EAMB = DF(1) /KI(1)
      ELSE
        DO I = 1, NS
          KI(I) =  0.D0
        ENDDO
        EAMB = 0.D0
      ENDIF
 
C     Multicomponent Fick law
      DO I = 1, NS
        JDIF(I) = -WR1(IMI+I-1) /WR1(IUNA) *(1.D0 -Y(I))/(1.D0 -X(I))
     &            *NDIMI(I) *(DF(I) -KI(I) *EAMB) 
      ENDDO

C     Mass diffusion fluxes
      JPROJ = 0.D0
      DO J = 1, NS
        JPROJ = JPROJ +JDIF(J) 
      ENDDO
      DO J = 1, NS
        JDIFR(J) = JDIF(J) -Y(J)*JPROJ
      ENDDO

      END SUBROUTINE SMFICK
C----------------------------------------------------------------------- 
C-----------------------------------------------------------------------
      SUBROUTINE DIJFICK (WR1, LWR1, WR2, LWR2, WI, LWI, X, ND, DIJ)
C-----------------------------------------------------------------------
C     This subroutine computes the diffusion coefficient, Ranshaw's 
C     approximation. 
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), WR2(1:LWR2), X(1:NS), ND,
     &                 DIJ(1:NS,1:NS)
C-----------------------------------------------------------------------
      INTEGER I, J
      DOUBLE PRECISION MM, Y(1:NS), NDIMI(1:NS)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
C     Molar mass
      MM = 0.D0
      DO I = 1, NS
        MM = MM +X(I) *WR1(IMI+I-1)
      ENDDO

C     System
      CALL SYSTSMFICKY (NS, WR1(IMI), WR2(IBINIJ), MM, X, Y, NDIMI) 

C     Multicomponent Diffusion coefficient 
      DO I = 1, NS
        DO J = 1, NS
          DIJ(I,J) = -Y(J) *NDIMI(J) /ND /X(J)
        ENDDO
        DIJ(I,I) = DIJ(I,I) + NDIMI(I) /ND /X(I) 
      ENDDO

      END SUBROUTINE DIJFICK
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SYSTSMFICKY (N, MI, BINIJ, MM, X, Y, NDIMI) 
C-----------------------------------------------------------------------
C     This subroutine computes the average diffusion coefficient
C     for Fick's law. 
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N 
      DOUBLE PRECISION MI(1:N), BINIJ(1:(N+1)*N/2), X(1:N), NDIMI(1:N),
     &                 MM, Y(1:N)
C-----------------------------------------------------------------------
      INTEGER I, J, IJ, II, JJ, IJBIN
      DOUBLE PRECISION SUM1, BINIJMAT(1:N,1:N) 
C-----------------------------------------------------------------------
      DO I = 1, N
          Y(I) = X(I) *MI(I) /MM
      ENDDO

      DO I = 1, N
        DO J = I+1, N
            IJ    = ((I-1)*(2*N-I)+2*J)/2
            BINIJMAT(I,J) = BINIJ(IJ)
            BINIJMAT(J,I) = BINIJ(IJ)
            BINIJMAT(I,I) = 0.D0 
        ENDDO
      ENDDO

      DO I = 1, N
        SUM1 = 0.D0
        DO J = 1, N
          IF (J/=I) THEN
            SUM1 = SUM1 +X(J) / BINIJMAT(I,J)
          ENDIF
        ENDDO
        NDIMI(I) = (1.D0 - Y(I)) / SUM1
      ENDDO

      END SUBROUTINE SYSTSMFICKY
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SYSTSMFICK (N, MI, BINIJ, MM, X, Y, NDIMI) 
C-----------------------------------------------------------------------
C     This subroutine computes the average diffusion coefficient
C     for Fick's law. 
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N 
      DOUBLE PRECISION MI(1:N), BINIJ(1:(N+1)*N/2), X(1:N), NDIMI(1:N),
     &                 MM, Y(1:N)
C-----------------------------------------------------------------------
      INTEGER I, J, IJ, II, JJ, IJBIN
      DOUBLE PRECISION SUM1, BINIJMAT(1:N,1:N) 
C-----------------------------------------------------------------------
      DO I = 1, N
          Y(I) = X(I) *MI(I) /MM
      ENDDO

      DO I = 1, N
        DO J = I+1, N
            IJ    = ((I-1)*(2*N-I)+2*J)/2
            BINIJMAT(I,J) = BINIJ(IJ)
            BINIJMAT(J,I) = BINIJ(IJ)
            BINIJMAT(I,I) = 0.D0 
        ENDDO
      ENDDO

      DO I = 1, N
        SUM1 = 0.D0
        DO J = 1, N
          IF (J/=I) THEN
            SUM1 = SUM1 +X(J) / BINIJMAT(I,J)
          ENDIF
        ENDDO
        NDIMI(I) = (1.D0 - X(I)) / SUM1
      ENDDO

      END SUBROUTINE SYSTSMFICK
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SMRAMD (WR1, LWR1, WR2, LWR2, WI, LWI, X, TH, TE, ND,
     &                   DF, FIJ, JDIF, EAMB)
C-----------------------------------------------------------------------
C     This subroutine computes the mass diffusion fluxes solution of the
C     Ramshaw equations in the case of a mixture with charged 
C     particles (ambipolar constraint) in thermal non-equilibrium. If 
C     considered, thermal diffusion is included in the driving forces.
C     A variant of Choleki's method is used to solve for the linear
C     system for the driving forces.
c     The ambipolar field is deduced from the electron driving force.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), WR2(1:LWR2), X(1:NS), TH, TE, ND,
     &                 DF(1:NS), FIJ(1:NS*(NS-1)/2), JDIF(1:NS), EAMB
C-----------------------------------------------------------------------
      DOUBLE PRECISION G(1:(NS-1)*NS/2), BETA(1:NS-1), MM, TEMP(1:NS-1)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
C     Molar mass
      MM = 0.D0
      DO I = 1, NS
        MM = MM +X(I) *WR1(IMI+I-1)
      ENDDO

C     System
      CALL SYSTSMRAMD (NS, WR1(IMI), WR2(IBINIJ), WI(IQI), X, ND, DF, 
     &                 FIJ, G, BETA, MM, WR1(IUKB), TH, WR1(IUE), EAMB)

C     LU decomposition
      CALL EGSDEC (NS-1, G, TEMP, IER)

C     Driving forces contribution
      CALL EGSSOL (NS-1, G, BETA)

      VE = 0.D0
      DO I = 2, NS
        VE = VE -X(I) *WI(IQI+I-1) *BETA(I-1)
      ENDDO
      VE = VE /(X(1) *WI(IQI))


C     Mass diffusion fluxes
      I = 1
        JDIF(I) = VE *X(I) *WR1(IMI+I-1) /WR1(IUNA) *ND
      DO I = 2, NS
        JDIF(I) = BETA(I-1) *X(I) *WR1(IMI+I-1) /WR1(IUNA) *ND
      ENDDO

      END SUBROUTINE SMRAMD
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SYSTSMRAMD (N, MI, BINIJ, QI, X, ND, DF, FIJ, G, BETA,
     &                       MM, KB, TH, E, EAMB)
C-----------------------------------------------------------------------
C     This subroutine computes the system for the mass diffusion fluxes
C     using the Ramshaw equations in the case of a mixture with charged
C     particles (ambipolar field). If considered, thermal diffusion is 
C     included in the driving forces.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N, QI(1:N)
      DOUBLE PRECISION MI(1:N), BINIJ(1:(N+1)*N/2), X(1:N), ND,
     &                 DF(1:N), FIJ(1:N*(N-1)/2), G(1:(N-1)*N/2), 
     &                 BETA(1:N-1), MM, KB, TH, E, EAMB
C-----------------------------------------------------------------------
      INTEGER I, J, IJ, II, JJ, IJBIN
      DOUBLE PRECISION FAC, QT, DIJAV, M(1:N-1), KA(1:N)
C-----------------------------------------------------------------------
      QT = 0.D0
      DO I = 1, N
        QT = QT +X(I) *QI(I)
      ENDDO

      DIJAV = 0.D0
      DO I = 1, N
        KA(I)   =  X(I) *(QI(I) -MI(I) /MM *QT) *E /(KB *TH)
        DO J = I+1, N
          IJ    = ((I-1)*(2*N-I)+2*J)/2
          DIJAV = DIJAV +BINIJ(IJ) /ND
        ENDDO
      ENDDO
      DIJAV = DIJAV /(N*(N-1)/2)

      EAMB = DF(1) /KA(1)
      DO I = 1, N-1
        II       = ((I-1)*(2*(N-1)-I)+2*I)/2
        G(II)    = 0.D0
        BETA(I)  = -DF(I+1) +KA(I+1) *EAMB
        M(I)     = X(I+1) /MM *(MI(I+1) -MI(1) *QI(I+1) /QI(1))
      ENDDO

      DO I = 1, N-1
        DO J = I+1, N-1
          IJBIN = (I*(2*N-I-1)+2*(J+1))/2
          IJ    = ((I-1)*(2*(N-1)-I)+2*J)/2
          II    = ((I-1)*(2*(N-1)-I)+2*I)/2
          JJ    = ((J-1)*(2*(N-1)-J)+2*J)/2
          FAC   = X(I+1)*X(J+1) *ND /(BINIJ(IJBIN) 
     &            *FIJ((I*(2*N-I-1))/2+J-I))
          G(IJ) = -FAC
          G(II) = G(II) + FAC
          G(JJ) = G(JJ) + FAC
        ENDDO
      ENDDO

C     Mass constraint
      DO I = 1, N-1
        DO J = I, N-1
          IJ = ((I-1)*(2*(N-1)-I)+2*J)/2
          G(IJ) = G(IJ) +M(I) *M(J) /DIJAV
        ENDDO
      ENDDO

      END SUBROUTINE SYSTSMRAMD
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SMRAMCG (WR1, LWR1, WR2, LWR2, WI, LWI, X, TH, TE, ND,
     &                    DF, FIJ, JDIF, EAMB)
C-----------------------------------------------------------------------
C     This subroutine computes the mass diffusion fluxes solution of the
C     Ramshaw equations in the case of a mixture with charged 
C     particles (ambipolar constraint) in thermal non-equilibrium. If 
C     considered, thermal diffusion is included in the driving forces.
C     The CG method is used to solve for the linear system for the 
C     driving forces.
c     The ambipolar field is deduced from the electron driving force.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), WR2(1:LWR2), X(1:NS), TH, TE, ND,
     &                 DF(1:NS), FIJ(1:NS*(NS-1)/2), JDIF(1:NS), EAMB
C-----------------------------------------------------------------------
      DOUBLE PRECISION G(1:(NS-1)*NS/2), BETA(1:NS-1), MM, 
     &                 MX, P(1:NS-1), V(1:NS-1), R(1:NS-1), M(1:NS-1),
     &                 TOLRES
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
C     Molar mass
      MM = 0.D0
      DO I = 1, NS
        MM = MM +X(I) *WR1(IMI+I-1)
      ENDDO

C     System
      CALL SYSTSMRAMCG (NS, WR1(IMI), WR2(IBINIJ), WI(IQI), X, ND, DF, 
     &                  FIJ, G, BETA, MM, WR1(IUKB), TH, WR1(IUE), EAMB,
     &                  R, M)

C     Krylov iterative solution
      DO I = 1, NS-1
        II = ((I-1)*(2*(NS-1)-I)+2*I)/2
        P(I) = G(II)
      ENDDO
      ITERMAX = 3; TOLRES = 2.D-2
      CALL CG (NS-1, G, P, V, BETA, TOLRES, ITERMAX, ITER)

C     Constraint
      RM = 0.D0; MX = 0.D0
      DO I = 1, NS-1
        RM = RM +M(I) *R(I)
        MX = MX +M(I) *V(I)
      ENDDO
      PROJ = MX /RM
      DO I = 1, NS-1
       V(I) = V(I) -PROJ *R(I)
      ENDDO


      VE = 0.D0
      DO I = 2, NS
        VE = VE -X(I) *WI(IQI+I-1) *V(I-1)
      ENDDO
      VE = VE /(X(1) *WI(IQI))


C     Mass diffusion fluxes
      I = 1
        JDIF(I) = VE *X(I) *WR1(IMI+I-1) /WR1(IUNA) *ND
      DO I = 2, NS
        JDIF(I) = V(I-1) *X(I) *WR1(IMI+I-1) /WR1(IUNA) *ND
      ENDDO

      END SUBROUTINE SMRAMCG
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SYSTSMRAMCG (N, MI, BINIJ, QI, X, ND, DF, FIJ, G, BETA,
     &                        MM, KB, TH, E, EAMB, R, M)
C-----------------------------------------------------------------------
C     This subroutine computes the system for the mass diffusion fluxes
C     using the Ramshaw equations in the case of a mixture with charged
C     particles (ambipolar field). If considered, thermal diffusion is 
C     included in the driving forces.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N, QI(1:N)
      DOUBLE PRECISION MI(1:N), BINIJ(1:(N+1)*N/2), X(1:N), ND,
     &                 DF(1:N), FIJ(1:N*(N-1)/2), G(1:(N-1)*N/2), 
     &                 BETA(1:N-1), MM, KB, TH, E, EAMB, R(1:N-1),
     &                 M(1:N-1)
C-----------------------------------------------------------------------
      INTEGER I, J, IJ, II, JJ, IJBIN
      DOUBLE PRECISION FAC, QT, KA(1:N)
C-----------------------------------------------------------------------
      QT = 0.D0
      DO I = 1, N
        QT = QT +X(I) *QI(I)
      ENDDO

      DO I = 1, N
        KA(I)   =  X(I) *(QI(I) -MI(I) /MM *QT) *E /(KB *TH)
      ENDDO

      EAMB = DF(1) /KA(1)
      DO I = 1, N-1
        II       = ((I-1)*(2*(N-1)-I)+2*I)/2
        G(II)    = 0.D0
        BETA(I)  = -DF(I+1) +KA(I+1) *EAMB
        M(I)     = X(I+1) /MM *(MI(I+1) -MI(1) *QI(I+1) /QI(1))
        R(I)     = 1.D0
      ENDDO

      DO I = 1, N-1
        DO J = I+1, N-1
          IJBIN = (I*(2*N-I-1)+2*(J+1))/2
          IJ    = ((I-1)*(2*(N-1)-I)+2*J)/2
          II    = ((I-1)*(2*(N-1)-I)+2*I)/2
          JJ    = ((J-1)*(2*(N-1)-J)+2*J)/2
          FAC   = X(I+1)*X(J+1) *ND /(BINIJ(IJBIN) 
     &            *FIJ((I*(2*N-I-1))/2+J-I))
          G(IJ) = -FAC
          G(II) = G(II) + FAC
          G(JJ) = G(JJ) + FAC
        ENDDO
      ENDDO

      END SUBROUTINE SYSTSMRAMCG
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SMRAMSUT (WR1, LWR1, WR2, LWR2, WI, LWI, X, TH, TE, ND,
     &                     DF1, FIJ, JDIF, EAMB)
C-----------------------------------------------------------------------
C     This subroutine computes the mass diffusion fluxes solution of the
C     Ramshaw equations in the case of a mixture with charged 
C     particles (ambipolar constraint) in thermal non-equilibrium. If 
C     considered, thermal diffusion is included in the driving forces.
C     The Sutton & Gnoffo algorithm is used. 
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), WR2(1:LWR2), X(1:NS), TH, TE, ND,
     &                 DF1(1:NS), FIJ(1:NS*(NS-1)/2), JDIF(1:NS), EAMB
C-----------------------------------------------------------------------
      DOUBLE PRECISION D(2:NS,2:NS), DEF(2:NS), NOLD(2:NS), NNEW(2:NS),
     &                 RES(2:NS), MM, TOLRES, QT, KA(1:NS), DF(2:NS)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
C     Ambipolar field
      MM = 0.D0; QT = 0.D0
      DO I = 1, NS
        MM = MM +X(I) *WR1(IMI+I-1)
        QT = QT +X(I) *WI(IQI+I-1)
      ENDDO
      DO I = 1, NS
        KA(I)   =  X(I) *(WI(IQI+I-1) -WR1(IMI+I-1) /MM *QT) *WR1(IUE) 
     &             /(WR1(IUKB) *TH)
      ENDDO
      EAMB = DF1(1) /KA(1)

C     Binary diffusion coefficient. 
      DO I = 2, NS
        DF(I)  = DF1(I) -KA(I) *EAMB
        II     = ((I-1)*(2*NS-I)+2*I)/2
        D(I,I) = WR2(IBINIJ+II-1) /ND
        DO J = I+1, NS
          IJ     = ((I-1)*(2*NS-I)+2*J)/2
          D(I,J) = WR2(IBINIJ+IJ-1) /ND *FIJ(((I-1)*(2*NS-I))/2+J-I)
          D(J,I) = D(I,J)
        ENDDO
      ENDDO

C     Effective binary diffusion coefficient
      DO I = 2, NS
        SUM = 0.D0
        DO J = 2, NS 
          SUM = SUM +X(J) /D(I,J)
        ENDDO
        DEF(I) = 1.D0 /SUM
      ENDDO

C     Initial conditions, number density diffusion fluxes set to zero
      RES2INI = 0.D0
      DO I = 2, NS
        NNEW(I) = 0.D0
        NOLD(I) = NNEW(I)
        RES(I)  = -ND *DEF(I) *DF(I)
        RES2INI = RES2INI +RES(I) *RES(I)
      ENDDO
      RES2 = 1.D0

C     -Iterative process, number density diffusion fluxes
      TOLRES = 1.D-3; ITERMAX = 8 
      TOLRES2 = TOLRES *TOLRES
      ITER = 0
      DO WHILE ( (ITER< ITERMAX) .AND. (RES2 > TOLRES2) )
        ITER = ITER +1

C     -Driving force contribution to residual
        RES2 = 0.D0
        DO I = 2, NS
          SUM = 0.D0
          DO J = 2, NS 
            SUM = SUM +NOLD(J) /D(I,J)
          ENDDO
          RES(I)  = -NOLD(I) -ND *DEF(I) *DF(I) +X(I) *DEF(I) *SUM
          RES2 = RES2 +RES(I) *RES(I)
          NNEW(I) = NOLD(I) +RES(I)
        ENDDO
        RES2 = RES2 /RES2INI

C     -Correction for mass constraint
        DO I = 2, NS
          SUM = 0.D0
          DO J = 2, NS
            SUM = SUM +WR1(IMI+J-1) *NNEW(J)
          ENDDO
          NNEW(I) = NNEW(I) -X(I) /MM *SUM
        ENDDO

C     -Update
        DO I = 2, NS
          NOLD(I) = NNEW(I)
        ENDDO

      ENDDO
 
C     Mass diffusion fluxes
      JDIF(1) = 0.D0
      DO I = 2, NS
        JDIF(I) = WR1(IMI+I-1) /WR1(IUNA) *NNEW(I)
        JDIF(1) = JDIF(1) -WI(IQI+I-1) *JDIF(I) /WR1(IMI+I-1)
      ENDDO
      JDIF(1) = JDIF(1) *WR1(IMI) /WI(IQI)

      END SUBROUTINE SMRAMSUT
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C     Neutral mixture, system for V_i

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SMNEUTD (WR1, LWR1, WR2, LWR2, X, ND, DF, FIJ, JDIF)
C-----------------------------------------------------------------------
C     This subroutine computes the mass diffusion fluxes solution of the
C     Stefan-Maxwell equations in the case of a neutral mixture. If
C     considered, thermal diffusion is included in the driving forces.
C     A variant of Choleki's method is used to solve for the
C     linear system for the driving forces.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2
      DOUBLE PRECISION WR1(1:LWR1), WR2(1:LWR2), X(1:NS), ND,
     &                 DF(1:NS), FIJ(1:NS*(NS-1)/2), JDIF(1:NS)
C-----------------------------------------------------------------------
      DOUBLE PRECISION G(1:(NS+1)*NS/2), BETA(1:NS), MM, TEMP(1:NS)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
C     Molar mass
      MM = 0.D0
      DO I = 1, NS
        MM = MM +X(I) *WR1(IMI+I-1)
      ENDDO

C     System
      CALL SYSTSMNEUTD (NS, WR1(IMI), WR2(IBINIJ), X, ND, DF, FIJ, G, 
     &                  BETA, MM)

C     LU decomposition
      CALL EGSDEC (NS, G, TEMP, IER)

C     Driving forces contribution
      CALL EGSSOL (NS, G, BETA)

C     Mass diffusion fluxes
      DO I = 1, NS
       JDIF(I) = BETA(I) *X(I) *WR1(IMI+I-1) /WR1(IUNA) *ND
      ENDDO

      END SUBROUTINE SMNEUTD
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SYSTSMNEUTD (N, MI, BINIJ, X, ND, DF, FIJ, G, BETA, MM)
C-----------------------------------------------------------------------
C     This subroutine computes the system for the mass diffusion fluxes
C     using the Stefan-Maxwell equations in the case of a neutral 
C     mixture. If considered, thermal diffusion is included in the 
C     driving forces.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION MI(1:N), BINIJ(1:(N+1)*N/2), X(1:N), ND,
     &                 DF(1:N), FIJ(1:N*(N-1)/2), G(1:(N+1)*N/2), 
     &                 BETA(1:N), MM
C-----------------------------------------------------------------------
      INTEGER I, J, IJ, II, JJ
      DOUBLE PRECISION FAC, DIJAV, M(1:N)
C-----------------------------------------------------------------------
      DO I = 1, N
        BETA(I) = -DF(I)
        M(I)     = X(I) *MI(I) /MM
        II      = ((I-1)*(2*N-I)+2*I)/2
        G(II)   = 0.D0
      ENDDO
      BETA(1) = BETA(1) 

      DIJAV = 0.D0
      DO I = 1, N
        DO J = I+1, N
          IJ = ((I-1)*(2*N-I)+2*J)/2
          DIJAV = DIJAV +BINIJ(IJ) /ND
        ENDDO
      ENDDO
      DIJAV = DIJAV /(N*(N-1)/2)

      DO I = 1, N
        DO J = I+1, N
          IJ = ((I-1)*(2*N-I)+2*J)/2
          II = ((I-1)*(2*N-I)+2*I)/2
          JJ = ((J-1)*(2*N-J)+2*J)/2
          FAC   = X(I)*X(J) *ND /(BINIJ(IJ) *FIJ(((I-1)*(2*N-I))/2+J-I))
          G(IJ) = -FAC
          G(II) = G(II) + FAC
          G(JJ) = G(JJ) + FAC
        ENDDO
      ENDDO

C     Mass constraint
      DO I = 1, N
        DO J = I, N
          IJ = ((I-1)*(2*N-I)+2*J)/2
          G(IJ) = G(IJ) +M(I) *M(J) /DIJAV
        ENDDO
      ENDDO

      END SUBROUTINE SYSTSMNEUTD
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SMNEUTCG (WR1, LWR1, WR2, LWR2, X, ND, DF, FIJ, JDIF)
C-----------------------------------------------------------------------
C     This subroutine computes the mass diffusion fluxes solution of the
C     Stefan-Maxwell equations in the case of a neutral mixture. If
C     considered, thermal diffusion is included in the driving forces.
C     The CG method is used to solve for the linear system for the 
C     driving forces.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2
      DOUBLE PRECISION WR1(1:LWR1), WR2(1:LWR2), X(1:NS), ND, DF(1:NS),
     &                 FIJ(1:NS*(NS-1)/2), JDIF(1:NS)
C-----------------------------------------------------------------------
      DOUBLE PRECISION G(1:(NS+1)*NS/2), BETA(1:NS), MM, TEMP(1:NS), 
     &                 R(1:NS), M(1:NS), P(1:NS), V(1:NS), MX, TOLRES
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
C     Molar mass
      MM = 0.D0
      DO I = 1, NS
        MM = MM +X(I) *WR1(IMI+I-1)
      ENDDO

C     System
      CALL SYSTSMNEUTCG (NS, WR1(IMI), WR2(IBINIJ), X, ND, DF, FIJ, G, 
     &                   BETA, MM, M, R)


C     Krylov iterative solution
      DO I = 1, NS
        II = ((I-1)*(2*NS-I)+2*I)/2
        P(I) = G(II)
      ENDDO
      ITERMAX = 3; TOLRES = 1.D-2
      CALL CG (NS, G, P, V, BETA, TOLRES, ITERMAX, ITER)

C     Constraint
      RM = 0.D0; MX = 0.D0
      DO I = 1, NS
        RM = RM +M(I) *R(I)
        MX = MX +M(I) *V(I)
      ENDDO
      PROJ = MX /RM
      DO I = 1, NS
       V(I) = V(I) -PROJ *R(I)
      ENDDO

C     Mass diffusion fluxes
      DO I = 1, NS
       JDIF(I) = V(I) *X(I) *WR1(IMI+I-1) /WR1(IUNA) *ND
      ENDDO

      END SUBROUTINE SMNEUTCG
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SYSTSMNEUTCG (N, MI, BINIJ, X, ND, DF, FIJ, G, BETA, 
     &                         MM, M, R)
C-----------------------------------------------------------------------
C     This subroutine computes the system for the mass diffusion fluxes
C     using the Stefan-Maxwell equations in the case of a neutral 
C     mixture. If considered, thermal diffusion is included in the 
C     driving forces.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION MI(1:N), BINIJ(1:(N+1)*N/2), X(1:N), ND,
     &                 DF(1:N), FIJ(1:N*(N-1)/2), G(1:(N+1)*N/2), 
     &                 BETA(1:N), MM, M(1:N), R(1:N)
C-----------------------------------------------------------------------
      INTEGER I, J, IJ, II, JJ
      DOUBLE PRECISION FAC
C-----------------------------------------------------------------------
      DO I = 1, N
        BETA(I) = -DF(I)
        R(I)    = 1.D0
        M(I)    = X(I) *MI(I) /MM
        II      = ((I-1)*(2*N-I)+2*I)/2
        G(II)   = 0.D0
      ENDDO
      BETA(1) = BETA(1) 

      DO I = 1, N
        DO J = I+1, N
          IJ = ((I-1)*(2*N-I)+2*J)/2
          II = ((I-1)*(2*N-I)+2*I)/2
          JJ = ((J-1)*(2*N-J)+2*J)/2
          FAC   = X(I)*X(J) *ND /(BINIJ(IJ) *FIJ(((I-1)*(2*N-I))/2+J-I))
          G(IJ) = -FAC
          G(II) = G(II) + FAC
          G(JJ) = G(JJ) + FAC
        ENDDO
      ENDDO

      END SUBROUTINE SYSTSMNEUTCG
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SMNEUTSUT (WR1, LWR1, WR2, LWR2, X, ND, DF, FIJ, JDIF)
C-----------------------------------------------------------------------
C     This subroutine computes the mass diffusion fluxes solution of the
C     Stefan-Maxwell equations in the case of a mixture without charged
C     particles (no ambipolar field). If considered, thermal diffusion
C     is included in the driving forces.
C     Sutton's algorithm is used to solve for the non-symmetric linear
C     system.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2
      DOUBLE PRECISION WR1(1:LWR1), WR2(1:LWR2), X(1:NS), ND, DF(1:NS),
     &                 FIJ(1:NS*(NS-1)/2), JDIF(1:NS)
C-----------------------------------------------------------------------
      DOUBLE PRECISION D(1:NS,1:NS), DEF(1:NS), NOLD(1:NS), NNEW(1:NS),
     &                 RES(1:NS), MM, TOLRES
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
C     Binary diffusion coefficient. Sutton's algorithm requires D_{ii}.
C     However, this coefficient should not play any role at convergence.
      DO I = 1, NS
        II     = ((I-1)*(2*NS-I)+2*I)/2
        D(I,I) = WR2(IBINIJ+II-1) /ND *1.5D0
        DO J = I+1, NS
          IJ     = ((I-1)*(2*NS-I)+2*J)/2
          D(I,J) = WR2(IBINIJ+IJ-1) /ND *FIJ(((I-1)*(2*NS-I))/2+J-I)
          D(J,I) = D(I,J)
        ENDDO
      ENDDO

C     Mixture molar mass
      MM = 0.D0
      DO I = 1, NS
        MM = MM +X(I) *WR1(IMI+I-1)
      ENDDO

C     Effective binary diffusion coefficient
      DO I = 1, NS
        SUM = 0.D0
        DO J = 1, NS
          SUM = SUM +X(J) /D(I,J)
        ENDDO
        DEF(I) = 1.D0 /SUM
      ENDDO

C     Initial conditions, number density diffusion fluxes set to zero
      RES2INI = 0.D0
      DO I = 1, NS
        NNEW(I) = 0.D0
        NOLD(I) = NNEW(I)
        RES(I)  = -ND *DEF(I) *DF(I)
        RES2INI = RES2INI +RES(I) *RES(I)
      ENDDO
      RES2 = 1.D0

C     -Iterative process, number density diffusion fluxes
      TOLRES = 1.D-2; ITERMAX = 8 
      TOLRES2 = TOLRES *TOLRES
      ITER = 0
      DO WHILE ( (ITER< ITERMAX) .AND. (RES2 > TOLRES2) )
        ITER = ITER +1

C     -Residual
        RES2 = 0.D0
        DO I = 1, NS
          SUM = 0.D0
          DO J = 1, NS
            SUM = SUM +NOLD(J) /D(I,J)
          ENDDO
          RES(I)  = -NOLD(I) -ND *DEF(I) *DF(I) +X(I) *DEF(I) *SUM
          RES2 = RES2 +RES(I) *RES(I)
          NNEW(I) = NOLD(I) +RES(I)
        ENDDO
        RES2 = RES2 /RES2INI

C     -Correction for mass constraint
        DO I = 1, NS
          SUM = 0.D0
          DO J = 1, NS
            SUM = SUM +WR1(IMI+J-1) *NNEW(J)
          ENDDO
          NNEW(I) = NNEW(I) -X(I) /MM *SUM
        ENDDO

C     -Update
        DO I = 1, NS
          NOLD(I) = NNEW(I)
        ENDDO

      ENDDO

C     Mass diffusion fluxes
      DO I = 1, NS
        JDIF(I) = WR1(IMI+I-1) /WR1(IUNA) *NNEW(I)
      ENDDO

      END SUBROUTINE SMNEUTSUT
C-----------------------------------------------------------------------

