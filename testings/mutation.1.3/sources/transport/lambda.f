C-----------------------------------------------------------------------
      SUBROUTINE LAMBDACHID (WR1, LWR1, WR2, LWR2, X, LAMBDA, CHI)
C-----------------------------------------------------------------------
C     This subroutine computes the translational thermal conductivity 
C     and the thermal diffusion ratios of the heavy particle gas. 
C     A direct method is used to solve the linear system.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2 
      DOUBLE PRECISION  WR1(1:LWR1), WR2(1:LWR2), X(1:NS), LAMBDA,
     &                  CHI(1:NS)
C-----------------------------------------------------------------------
      DOUBLE PRECISION G(1:(NS-NE+1)*(NS-NE)/2), BETA(1:NS-NE),
     &                 ALPHA(1:NS-NE), TEMP(1:NS-NE),
     &                 G01(1:NS-NE,1:NS-NE)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      CALL LAMBDACHISYST (NS-NE, WR1(IMI+NE), WR2(IBINIJ+NE*NS), 
     &                    WR2(IAIJ), WR2(IBIJ), WR2(ICIJ), WR2(IETAI), 
     &                    WR1(IUNA), WR1(IUKB), X(NE+1), G, BETA, G01)

      DO I = 1, NS-NE
        ALPHA(I) = BETA(I)
      ENDDO
      
      CALL EGSDEC (NS-NE, G, TEMP, IER)
      CALL EGSSOL (NS-NE, G, ALPHA)

      LAMBDA = 0.D0
      DO I = 1, NS-NE
        LAMBDA = LAMBDA +BETA(I) *ALPHA(I)
      ENDDO

      CHI(1) = 0.D0
C     Beware, the factor 5/2 is already in G01       
      CALL MATVEC(NS-NE, G01,ALPHA, CHI(NE+1))
      
      END SUBROUTINE LAMBDACHID
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE LAMBDACHICG (WR1, LWR1, WR2, LWR2, X, LAMBDA, CHI)
C-----------------------------------------------------------------------
C     This subroutine computes the translational thermal conductivity 
C     and the thermal diffusion ratios of the heavy particle gas. 
C     The conjugate gradient with preconditioning is used to solve the 
C     linear system.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2 
      DOUBLE PRECISION  WR1(1:LWR1), WR2(1:LWR2), X(1:NS), LAMBDA, 
     &                  CHI(1:NS)
C-----------------------------------------------------------------------
      DOUBLE PRECISION G(1:(NS-NE+1)*(NS-NE)/2), BETA(1:NS-NE),
     &                 ALPHA(1:NS-NE), TEMP(1:NS-NE), DMI(1:NS-NE), 
     &                 RN(1:NS-NE), G01(1:NS-NE,1:NS-NE), P(1:NS-NE)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      CALL LAMBDACHISYST (NS-NE, WR1(IMI+NE), WR2(IBINIJ+NE*NS), 
     &                    WR2(IAIJ), WR2(IBIJ), WR2(ICIJ), WR2(IETAI), 
     &                    WR1(IUNA), WR1(IUKB), X(NE+1), G, BETA, G01)

      DO I = 1, NS-NE
        RN(I) = BETA(I)
      ENDDO
 
      DO I = 1, NS-NE
        II = ((I-1)*(2*(NS-NE)-I)+2*I)/2
        P(I) = G(II)
      ENDDO

      ITERMAX = 3; TOLRES = 1.D-2 
      CALL CG (NS-NE, G, P, ALPHA, RN, TOLRES, ITERMAX, ITER)

      LAMBDA = 0.D0
      DO I = 1, NS-NE
        LAMBDA = LAMBDA +BETA(I) *ALPHA(I)
      ENDDO

      CHI(1) = 0.D0
      CALL MATVEC(NS-NE, G01, ALPHA, CHI(NE+1))

      END SUBROUTINE LAMBDACHICG
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE LAMBDACHISYST (N, MI, BINIJ, AIJ, BIJ, CIJ, ETAI, NA, 
     &                          KB, X, G, BETA, G01)
C-----------------------------------------------------------------------
C     This subroutine computes the system for the translational thermal 
C     conductivity and the thermal diffusion ratios of the heavy 
C     particle gas. 
C-----------------------------------------------------------------------
      IMPLICIT NONE 
      INTEGER N
      DOUBLE PRECISION MI(1:N), BINIJ(1:(N+1)*N/2), AIJ(1:(N+1)*N/2), 
     &                 BIJ(1:(N+1)*N/2), CIJ(1:(N+1)*N/2), ETAI(1:N), 
     &                 NA, KB, X(1:N), G(1:(N+1)*N/2), BETA(1:N),
     &                 G01(1:N,1:N)
C-----------------------------------------------------------------------
      DOUBLE PRECISION FAC1, FAC2, MIIJ, MJIJ
      INTEGER I, J, II, IJ, JJ
C-----------------------------------------------------------------------
      FAC1 = 4.D0 /(15.D0 *KB)
      DO I = 1, N
        BETA(I) =  X(I)
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
            FAC1 = FAC1 *5.D0 /2.D0 *(12.D0 *CIJ(IJ) -10.D0)
            G01(I,J) = FAC1 *MI(I)
            G01(I,I) = G01(I,I) -FAC1 *MI(J) 
          ENDDO
          DO J = I+1, N
            IJ = ((I-1)*(2*N-I)+2*J)/2 
            FAC1 = X(I) *X(J) /(BINIJ(IJ) *25.D0 *KB) /(MI(I) +MI(J))
            FAC1 = FAC1 *5.D0 /2.D0 *(12.D0 *CIJ(IJ) -10.D0)
            G01(I,J) = FAC1 *MI(I) 
            G01(I,I) = G01(I,I) -FAC1 *MI(J) 
          ENDDO
      ENDDO

      END SUBROUTINE LAMBDACHISYST
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE LAMBDAYOS (WR1, LWR1, WR2, LWR2, X, LAMBDA)
C-----------------------------------------------------------------------
C     This subroutine computes the translational thermal conductivity 
C     of the heavy particle gas. Yos' mixture rule is used.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2 
      DOUBLE PRECISION  WR1(1:LWR1), WR2(1:LWR2), X(1:NS), LAMBDA
C-----------------------------------------------------------------------
      INTEGER N
      DOUBLE PRECISION GIJ(1:(NS-NE+1)*(NS-NE)/2), GII(1:NS-NE)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      N = NS-NE
      CALL SYSTLAMBDAYOS (N, WR1(IMI+NE), WR2(IBINIJ+NE*NS), WR2(IAIJ),
     &                    WR2(IBIJ), WR1(IUKB), X(NE+1), GIJ, GII)


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

      LAMBDA = SUM /(1.D0 -GAV *SUM)

      END SUBROUTINE LAMBDAYOS
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SYSTLAMBDAYOS (N, MI, BINIJ, AIJ, BIJ, KB, X, GIJ, GII)
C-----------------------------------------------------------------------
C     This subroutine computes the system for the translational thermal 
C     conductivity of the heavy particle gas using Yos' rule. 
C-----------------------------------------------------------------------
      IMPLICIT NONE 
      INTEGER N
      DOUBLE PRECISION MI(1:N), BINIJ(1:(N+1)*N/2), AIJ(1:(N+1)*N/2), 
     &                 KB, X(1:N), GIJ(1:(N+1)*N/2), GII(1:N),
     &                 BIJ(1:(N+1)*N/2)
C-----------------------------------------------------------------------
      DOUBLE PRECISION B(1:N,1:N), MIJ, FAC
      INTEGER I, J, IJ
C-----------------------------------------------------------------------
      DO I = 1, N
        DO J = I, N
          IJ = ((I-1)*(2*N-I)+2*J)/2
          MIJ = MI(I) *MI(J) /(MI(I) +MI(J))**2
          GIJ(IJ) = MIJ /(75.D0 *KB *BINIJ(IJ))
     &            *(165.D0 -36.D0 *BIJ(IJ) -48.D0 *AIJ(IJ))
        ENDDO
      ENDDO

      DO I = 1, N
        DO J = 1, N 
          IF (J < I) THEN
            IJ = ((J-1)*(2*N-J)+2*I)/2
          ELSE
            IJ = ((I-1)*(2*N-I)+2*J)/2
          ENDIF
          MIJ = MI(I) *MI(J) /(MI(I) +MI(J))**2
          B(I,J) = MIJ /BINIJ(IJ) *(9.6D0 *AIJ(IJ) +(MI(I)-MI(J))
     &             *(9.D0/MI(J) +(-7.5D0 +3.6D0 *BIJ(IJ)) 
     &             /MI(I)))
        ENDDO
      ENDDO

      FAC = 2.D0 /(15.D0 *KB)
      DO I = 1, N
        GII(I) = 0.D0
        DO J = 1, N 
          GII(I) = GII(I) + X(J) *B(I,J) *FAC
        ENDDO
      ENDDO

      END SUBROUTINE SYSTLAMBDAYOS
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE LAMBDAWILKE (WR1, LWR1, WR2, LWR2, X, LAMBDA)
C-----------------------------------------------------------------------
C     This subroutine computes the translational thermal conductivity 
C     of the heavy particle gas. Wilke's mixture rule is used.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2 
      DOUBLE PRECISION  WR1(1:LWR1), WR2(1:LWR2), X(1:NS), LAMBDA
C-----------------------------------------------------------------------
      DOUBLE PRECISION KI(1:NS-NE), GI(1:NS-NE)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      N = NS-NE
      DO I = 1, N
        KI(I) = 3.75D0 *WR1(IUR) /WR1(IMI+NE+I-1) *WR2(IETAI+I-1)
      ENDDO
      CALL SYSTWILKE (N, WR1(IMI+NE), X(NE+1), GI, KI)

      LAMBDA = 0.D0      
      DO I = 1, N
        LAMBDA = LAMBDA + X(NE+I) *KI(I) /GI(I)
      ENDDO

      END SUBROUTINE LAMBDAWILKE
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE LAMBDAD (WR1, LWR1, WR2, LWR2, X, LAMBDA)
C-----------------------------------------------------------------------
C     This subroutine computes the translational thermal conductivity 
C     of the heavy particle gas. 
C     A direct method is used to solve the linear system.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2 
      DOUBLE PRECISION  WR1(1:LWR1), WR2(1:LWR2), X(1:NS), LAMBDA
C-----------------------------------------------------------------------
      DOUBLE PRECISION G(1:(NS-NE+1)*(NS-NE)/2), BETA(1:NS-NE),
     &                 ALPHA(1:NS-NE), TEMP(1:NS-NE)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      CALL LAMBDASYST (NS-NE, WR1(IMI+NE), WR2(IBINIJ+NE*NS), WR2(IAIJ),
     &                 WR2(IBIJ), WR2(ICIJ), WR2(IETAI), WR1(IUNA), 
     &                 WR1(IUKB), X(NE+1), G, BETA)

      DO I = 1, NS-NE
        ALPHA(I) = BETA(I)
      ENDDO
      
      CALL EGSDEC (NS-NE, G, TEMP, IER)
      CALL EGSSOL (NS-NE, G, ALPHA)

      LAMBDA = 0.D0
      DO I = 1, NS-NE
        LAMBDA = LAMBDA +BETA(I) *ALPHA(I)
      ENDDO

      END SUBROUTINE LAMBDAD
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE LAMBDACG (WR1, LWR1, WR2, LWR2, X, LAMBDA)
C-----------------------------------------------------------------------
C     This subroutine computes the translational thermal conductivity 
C     of the heavy particle gas. 
C     The conjugate gradient with preconditioning is used to solve the 
C     linear system.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2 
      DOUBLE PRECISION  WR1(1:LWR1), WR2(1:LWR2), X(1:NS), LAMBDA 
C-----------------------------------------------------------------------
      DOUBLE PRECISION G(1:(NS-NE+1)*(NS-NE)/2), BETA(1:NS-NE),
     &                 ALPHA(1:NS-NE), TEMP(1:NS-NE),DMI(1:NS-NE), 
     &                 ZN(1:NS-NE), RN(1:NS-NE)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      CALL LAMBDASYST (NS-NE, WR1(IMI+NE), WR2(IBINIJ+NE*NS), WR2(IAIJ),
     &                 WR2(IBIJ), WR2(ICIJ), WR2(IETAI), WR1(IUNA), 
     &                 WR1(IUKB), X(NE+1), G, BETA)

      DO I = 1, NS-NE
        RN(I) = BETA(I)
      ENDDO
 
      ITERMAX = 2 
      CALL EGSCG1 (NS-NE, G, DMI, ALPHA, ZN, RN, TEMP, ITERMAX)    

      LAMBDA = 0.D0
      DO I = 1, NS-NE
        LAMBDA = LAMBDA +BETA(I) *ALPHA(I)
      ENDDO

      END SUBROUTINE LAMBDACG
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE LAMBDASYST (N, MI, BINIJ, AIJ, BIJ, CIJ, ETAI, NA, KB, 
     &                       X, G, BETA)
C-----------------------------------------------------------------------
C     This subroutine computes the system for the translational thermal 
C     conductivity of the heavy particle gas. 
C-----------------------------------------------------------------------
      IMPLICIT NONE 
      INTEGER N
      DOUBLE PRECISION MI(1:N), BINIJ(1:(N+1)*N/2), AIJ(1:(N+1)*N/2), 
     &                 BIJ(1:(N+1)*N/2), CIJ(1:(N+1)*N/2), ETAI(1:N), 
     &                 NA, KB, X(1:N), G(1:(N+1)*N/2), BETA(1:N)
C-----------------------------------------------------------------------
      DOUBLE PRECISION FAC1, FAC2, MIIJ, MJIJ
      INTEGER I, J, II, IJ, JJ
C-----------------------------------------------------------------------
      FAC1 = 4.D0 /(15.D0 *KB)
      DO I = 1, N
        BETA(I) =  X(I)
        II = ((I-1)*(2*N-I)+2*I)/2 
        G(II) = X(I) * X(I) /ETAI(I) *MI(I) /NA *FAC1
      ENDDO
      DO I = 1, N
        DO J = I+1, N
          IJ = ((I-1)*(2*N-I)+2*J)/2 
          II = ((I-1)*(2*N-I)+2*I)/2 
          JJ = ((J-1)*(2*N-J)+2*J)/2 
          MIIJ = MI(I) /(MI(I) +MI(J)) 
          MJIJ = MI(J) /(MI(I) +MI(J))
          FAC2 = X(I) *X(J) /(BINIJ(IJ) *25.D0 *KB)
          G(IJ) = -FAC2 *MIIJ *MJIJ *(55.D0 -16.D0 *AIJ(IJ) 
     &            -12.D0 *BIJ(IJ))
          G(II) = G(II) + FAC2 *(30.D0 *MIIJ *MIIJ + 16.D0 *MIIJ *MJIJ 
     &            *AIJ(IJ) +(25.D0 -12.D0 *BIJ(IJ))*MJIJ*MJIJ)
          G(JJ) = G(JJ) + FAC2 *(30.D0 *MJIJ *MJIJ + 16.D0 *MIIJ *MJIJ 
     &            *AIJ(IJ) +(25.D0 -12.D0 *BIJ(IJ))*MIIJ*MIIJ) 
        ENDDO
      ENDDO

      END SUBROUTINE LAMBDASYST
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      SUBROUTINE KCONDUCTIVITY (WR1, LWR1, WR2, LWR2, WI, LWI, P, TH, X,
     &                          XTOL, XN, EPS, LAMBDATOT)
C-----------------------------------------------------------------------
C     This subroutine computes the thermal conductivity for a LTE
C     mixture. Barodiffusion is neglected.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), WR2(1:LWR2), P, TH, X(1:NS), 
     &                 XTOL(1:NS), XN(1:NC), EPS, LAMBDATOT
      INTEGER IV, IS
      DOUBLE PRECISION EPSP1, THP, TV(1:NV+1), TVP(1:NV+1),
     &                 EM1(1:NS), EM2(1:NS), EM3(1:NS), XP(1:NS),
     &                 EM4(1:NS), EM5(1:NS), EM6(1:NS), DF(1:NS),
     &                 EM1P(1:NS), EM2P(1:NS), EM3P(1:NS), 
     &                 EM4P(1:NS), EM5P(1:NS), EM6P(1:NS), JDIF(1:NS),
     &                 CPR, CPV, CPE, CPINT(1:NS), LAMBDAINTH, 
     &                 LAMBDATH, LAMBDATE, CHIH(1:NS), CHIE(1:NS),
     &                 FIJ(1:NS*(NS-1)/2), EAMB, ND, LAMBDAR
C-----------------------------------------------------------------------
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      EPSP1 = EPS +1.D0
      THP   = TH *EPSP1
      IF (NV==0) THEN
        TV(1)      = 1.D0
        TVP(1)     = TV(1)*EPSP1
      ELSE
        DO IV = 1, NV
          TV(IV)      = TH
          TVP(IV)     = TV(IV)*EPSP1 
        ENDDO
      ENDIF

C     -Eucken corrections
      CALL ENERGYMASS (WR1, LWR1, WI, LWI, TH, TH, TH, TV, P, EM1, EM2, 
     &                 EM3, EM4, EM5, EM6)
      CALL ENERGYMASS (WR1, LWR1, WI, LWI, THP, THP, THP, TVP, P, EM1P,
     &                 EM2P, EM3P, EM4P, EM5P, EM6P)
      DO IS = 1, NS
        CPE       = (EM3P(IS) -EM3(IS)) / (EPS *TH)
        CPR       = (EM4P(IS) -EM4(IS)) / (EPS *TH)
        CPV       = (EM5P(IS) -EM5(IS)) / (EPS *TH)
        CPINT(IS) = WR1(IMI+IS-1) *(CPE +CPR +CPV)
      ENDDO
      CALL LAMBDAINT (WR1, LWR1, WR2, LWR2, CPINT, XTOL, LAMBDAINTH)

C     -Heavy-particle translational thermal conductivity and thermal
C      diffusion ratios
      CALL LAMBDACHICG (WR1, LWR1, WR2, LWR2, XTOL, LAMBDATH, CHIH)
C                           CHIH(1:NS) with CHIH(1) = 0.D0 IF NE = 1 
C     -Electron translational thermal conductivity and thermal
C      diffusion ratios
      IF (NE /= 0) THEN
        CALL TDIFE (WR1, LWR1, WR2, LWR2, XTOL, TH, CHIE, 3)
        CALL LAMBDAE (WR1, LWR1, WR2, LWR2, XTOL, TH, LAMBDATE, 3)
      ELSE
        LAMBDATE = 0.D0
        DO IS = 1, NS
          CHIE(IS) = 0.D0
        ENDDO
      ENDIF

C     -The driving forces are concentration gradients due to unity 
C      temperature gradient (constant pressure).
      CALL COMPOSITION (WR1, LWR1, WI, LWI, THP, P, XN, X, XP)
      DO IS = 1, NS
        DF(IS) = (XP(IS) -X(IS)) /(TH *EPS) +(CHIE(IS) +CHIH(IS)) /TH 
      ENDDO 

C     -First-order Stefan-Maxwell equations for heavy particles
      CALL CORRECTION (WR1, LWR1, WR2, LWR2, XTOL, 1, FIJ)

C     -Stefan-Maxwell equations (Ramshaw): diffusion fluxes and 
C      ambipolar electric field 
      CALL NUMBERD (WR1, LWR1, P, TH, TH, X, ND)
      IF (NE /= 0) THEN
        CALL SMRAMCG (WR1, LWR1, WR2, LWR2, WI, LWI, XTOL, TH, TH, ND,
     &                DF, FIJ, JDIF, EAMB)
      ELSE
        CALL SMNEUTCG (WR1, LWR1, WR2, LWR2, XTOL, ND, DF, FIJ, JDIF) 
        EAMB = 0.D0
      ENDIF

C     -Reactive thermal conductivity based on the Stefan-Maxwell 
C      equations 
      LAMBDAR = 0.D0
      DO IS = 1, NS
        LAMBDAR  =  LAMBDAR -(EM1(IS) +2.D0/3.D0 *EM2(IS))* JDIF(IS)
      ENDDO

C     -Total thermal conductivity
      LAMBDATOT = LAMBDAINTH +LAMBDATH +LAMBDATE +LAMBDAR


      END SUBROUTINE KCONDUCTIVITY
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
      SUBROUTINE KCONDUCTIVITYNEQ (WR1, LWR1, WR2, LWR2, WI, LWI, P, TH,
     &                             X, XTOL, XN, EPS, LAMBDATOT)
C-----------------------------------------------------------------------
C     This subroutine computes the thermal conductivity for a mixture
C     in thermal equilibrium and chemical nonequilibrium. 
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), WR2(1:LWR2), P, TH, X(1:NS), 
     &                 XTOL(1:NS), XN(1:NC), EPS, LAMBDATOT
      INTEGER IV, IS
      DOUBLE PRECISION EPSP1, THP, TV(1:NV+1), TVP(1:NV+1),
     &                 EM1(1:NS), EM2(1:NS), EM3(1:NS), XP(1:NS),
     &                 EM4(1:NS), EM5(1:NS), EM6(1:NS), DF(1:NS),
     &                 EM1P(1:NS), EM2P(1:NS), EM3P(1:NS), 
     &                 EM4P(1:NS), EM5P(1:NS), EM6P(1:NS), JDIF(1:NS),
     &                 CPR, CPV, CPE, CPINT(1:NS), LAMBDAINTH, 
     &                 LAMBDATH, LAMBDATE, CHIH(1:NS), CHIE(1:NS),
     &                 FIJ(1:NS*(NS-1)/2), EAMB, ND, LAMBDAR
C-----------------------------------------------------------------------
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      EPSP1 = EPS +1.D0
      THP   = TH *EPSP1
      IF (NV==0) THEN
        TV(1)      = 1.D0
        TVP(1)     = TV(1)*EPSP1
      ELSE
        DO IV = 1, NV
          TV(IV)      = TH
          TVP(IV)     = TV(IV)*EPSP1 
        ENDDO
      ENDIF

C     -Eucken corrections
      CALL ENERGYMASS (WR1, LWR1, WI, LWI, TH, TH, TH, TV, P, EM1, EM2, 
     &                 EM3, EM4, EM5, EM6)
      CALL ENERGYMASS (WR1, LWR1, WI, LWI, THP, THP, THP, TVP, P, EM1P,
     &                 EM2P, EM3P, EM4P, EM5P, EM6P)
      DO IS = 1, NS
        CPE       = (EM3P(IS) -EM3(IS)) / (EPS *TH)
        CPR       = (EM4P(IS) -EM4(IS)) / (EPS *TH)
        CPV       = (EM5P(IS) -EM5(IS)) / (EPS *TH)
        CPINT(IS) = WR1(IMI+IS-1) *(CPE +CPR +CPV)
      ENDDO
      CALL LAMBDAINT (WR1, LWR1, WR2, LWR2, CPINT, XTOL, LAMBDAINTH)

C     -Heavy-particle translational thermal conductivity and thermal
C      diffusion ratios
      CALL LAMBDACHICG (WR1, LWR1, WR2, LWR2, XTOL, LAMBDATH, CHIH)
C                           CHIH(1:NS) with CHIH(1) = 0.D0 IF NE = 1 
C     -Electron translational thermal conductivity and thermal
C      diffusion ratios
      IF (NE /= 0) THEN
        CALL LAMBDAE (WR1, LWR1, WR2, LWR2, XTOL, TH, LAMBDATE, 3)
      ELSE
        LAMBDATE = 0.D0
      ENDIF

C     -Total thermal conductivity
      LAMBDATOT = LAMBDAINTH +LAMBDATH +LAMBDATE


      END SUBROUTINE KCONDUCTIVITYNEQ
C-----------------------------------------------------------------------

