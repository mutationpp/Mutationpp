C-----------------------------------------------------------------------
      SUBROUTINE COLLISION (WR1, LWR1, WR2, LWR2, TH, TE, ND, X) 
C-----------------------------------------------------------------------
C     This subroutine computes the collision integral ratios, modified 
C     binary diffusion coefficients (multiplied by the number density n) 
C     and species viscosity stored in the work array WR2.
C     This subroutine must be called by the user after every
C     (translational) temperature and composition (if electrons) update.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR1, LWR2
      DOUBLE PRECISION WR1(1:LWR1), WR2(1:LWR2), TH, TE, X(1:NS), ND
C-----------------------------------------------------------------------
      DOUBLE PRECISION MI, MJ, C(1:5), SQRT2,
     &                 LNTH1, LNTH2, LNTH3, LNTE1, LNTE2, LNTE3,
     &                 LNTSH1, LNTSH2, LNTSH3, LNTSH4,
     &                 LNTSE1, LNTSE2, LNTSE3, LNTSE4,
     &                 Q11EI(1:(NN+NI)), Q12EI(1:(NN+NI)), 
     &                 Q13EI(1:(NN+NI)), Q14EI(1:(NN+NI)), 
     &                 Q15EI(1:(NN+NI))
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      FAC1  = 3.D0 /16.D0 *DSQRT(2.D0 *WR1(IUPI) *WR1(IUR))
      SQRT2 = DSQRT(2.D0)
      LNTH1 = DLOG(TH)    ; LNTE1 = DLOG(TE)
      LNTH2 = LNTH1 *LNTH1; LNTE2 = LNTE1 *LNTE1
      LNTH3 = LNTH2 *LNTH1; LNTE3 = LNTE2 *LNTE1

C     1. Charge-charge interactions
      IF (NE /= 0) THEN
        FAC2 = FAC1 *DSQRT(TE / WR1(IMI)) 

C     -Average closest impact parameter (e-ion and e-e interactions)
        FAC0 = WR1(IUE)*WR1(IUE)/(8.D0*WR1(IUPI)*WR1(IUE0)*WR1(IUKB))
        BE   = FAC0 /TE

C     -Average closest impact parameter (ion-ion collision)
        BH   = FAC0 /TH

C     -Debye shielding distance (e and ion contribution)
        DS = DSQRT(WR1(IUE0) *WR1(IUKB) *TE 
     &            /(2.D0 *ND *X(1) *WR1(IUE) *WR1(IUE)))
        DS = MIN(DS, 10000.D0*(BE +BH))

C     -Non-dimensional temperature for charge-charge interactions
        TSE   = MAX(DS /(2.D0 *BE), 0.1D0)
        TSH   = MAX(DS /(2.D0 *BH), 0.1D0)
        LNTSE1 = DLOG(TSE)
        LNTSE2 = LNTSE1 *LNTSE1
        LNTSE3 = LNTSE2 *LNTSE1
        LNTSE4 = LNTSE3 *LNTSE1
        LNTSH1 = DLOG(TSH)
        LNTSH2 = LNTSH1 *LNTSH1
        LNTSH3 = LNTSH2 *LNTSH1
        LNTSH4 = LNTSH3 *LNTSH1
        EFAC  = WR1(IUPI) *DS *DS /(TSE *TSE)
        HFAC  = WR1(IUPI) *DS *DS /(TSH *TSH)

C     -Attractive interaction 
        C(1)=-1.67070199D-4; C(2)=4.87662165D-3; C(3)=-6.34929831D-2;
        C(4)=5.49816993D-1; C(5)=-7.91157696D-1;
        AEQ11 = EFAC *DEXP(C(1)*LNTSE4+C(2)*LNTSE3+C(3)*LNTSE2
     &          +C(4)*LNTSE1+C(5))
        C(1)=2.15877964D-5; C(2)=-6.06756574D-4; C(3)=7.22979820D-3;
        C(4)=-4.82251886D-2; C(5)=5.21304137D-1;
        AEQ12 = AEQ11*MAX(C(1)*LNTSE4+C(2)*LNTSE3+C(3)*LNTSE2
     &          +C(4)*LNTSE1+C(5),1.D0/3.D0)
        C(1)=1.59680456D-5; C(2)=-5.13956309D-4; C(3)=7.74980901D-3;
        C(4)=-6.71047378D-2; C(5)=1.32373242D+0;
        AEQ13 = 1.25D0 *AEQ12 -0.25D0 *AEQ11 *MAX(C(1)*LNTSE4
     &          +C(2)*LNTSE3+C(3)*LNTSE2+C(4)*LNTSE1+C(5),1.D0)
        C(1)=-1.85869672D-4; C(2)=4.45192523D-3; C(3)=-4.77048273D-2;
        C(4)=3.89672330D-1; C(5)=-2.31134796D+0;
        AEQ14 = EFAC *DEXP(C(1)*LNTSE4+C(2)*LNTSE3+C(3)*LNTSE2
     &          +C(4)*LNTSE1+C(5))
        C(1)=-1.91341716D-4; C(2)=4.39527014D-3; C(3)=-4.53877276D-2;
        C(4)=3.69204436D-1; C(5)=-2.62434507D+0;
        AEQ15 = EFAC *DEXP(C(1)*LNTSE4+C(2)*LNTSE3+C(3)*LNTSE2
     &          +C(4)*LNTSE1+C(5))

C     -Repulsive interaction 
        C(1)=-3.25157374D-4; C(2)=9.56207818D-3; C(3)=-1.17684421D-1;
        C(4)=8.41093450D-1; C(5)=-1.39975085D+0;
        REQ11 = EFAC *DEXP(C(1)*LNTSE4+C(2)*LNTSE3+C(3)*LNTSE2
     &          +C(4)*LNTSE1+C(5))
        C(1)=3.41942D-3;   C(2)=-8.01279D-2;C(3)=7.5847D-1;
        C(4)=-1.3513D0
        RHQ11 = HFAC *DEXP(C(1)*LNTSH3+C(2)*LNTSH2+C(3)*LNTSH1+C(4))
        C(1)=-3.96530664D-4; C(2)=1.08344320D-2; C(3)=-1.22277658D-1;
        C(4)=8.03918135D-1; C(5)=-1.10681226D+0;
        REQ22 = EFAC *DEXP(C(1)*LNTSE4+C(2)*LNTSE3+C(3)*LNTSE2
     &          +C(4)*LNTSE1+C(5))
        RHQ22 = HFAC *DEXP(C(1)*LNTSH4+C(2)*LNTSH3+C(3)*LNTSH2
     &          +C(4)*LNTSH1+C(5))
        C(1)=2.14301076D-5; C(2)=-7.35578975D-4; C(3)=9.78964744D-3;
        C(4)=-6.37201837D-2; C(5)=7.01123442D-1;
        REQ23 = REQ22 *MAX(C(1)*LNTSE4+C(2)*LNTSE3+C(3)*LNTSE2
     &          +C(4)*LNTSE1+C(5), 0.5D0)
        C(1)=-3.97192162D-4; C(2)=1.00233388D-2; C(3)=-1.03318360D-1;
        C(4)=6.47397675D-1; C(5)=-1.75549845D+0;
        REQ24 = EFAC *DEXP(C(1)*LNTSE4+C(2)*LNTSE3+C(3)*LNTSE2
     &          +C(4)*LNTSE1+C(5))
        RHA   = RHQ22/RHQ11
        C(1)=3.1204D-6;  C(2)=-0.0002626D0; C(3)=0.007332D0;
        C(4)=-0.0869D0;  C(5)=1.4292D0
        RHB   = MAX(C(1)*LNTSH4+C(2)*LNTSH3+C(3)*LNTSH2
     &          +C(4)*LNTSH1+C(5),1.D0)
        C(1)=4.8891D-7;  C(2)=-0.0001389D0; C(3)=0.004867D0;
        C(4)=-0.05934D0; C(5)=0.6008D0
        RHC   = MAX(C(1)*LNTSH4+C(2)*LNTSH3+C(3)*LNTSH2
     &          +C(4)*LNTSH1+C(5),1.D0/3.D0)

C     2. Electron gas
C     2.1 Electron-electron interactions
        WR2(IQEE22) = REQ22
        Q22EE = REQ22 
        Q23EE = REQ23
        Q24EE = REQ24
        WR2(IBINIJ) = FAC2 /REQ11

C     2.2 Electron-neutral interactions
        DO J = 1, NN
          Q11EI(J) = 1.D-20*DEXP( LNTE3*WR1(IFITQ11EI+4*J-4)
     &                           +LNTE2*WR1(IFITQ11EI+4*J-3)
     &                           +LNTE1*WR1(IFITQ11EI+4*J-2)
     &                                 +WR1(IFITQ11EI+4*J-1))
C-----------------------------------------------------------------------
C       To be modified (B* = 0. to avoid problems in e lambda_e, \xi=3): 
          Q12EI(J)= Q11EI(J)
          Q13EI(J)= 1.25D0 *Q12EI(J) -0.25D0 *Q11EI(J)
     &              *1.D-20*DEXP( LNTE3*WR1(IFITQ13EI+4*J-4)
C     &              *DEXP( LNTE3*WR1(IFITQ13EI+4*J-4)
     &                    +LNTE2*WR1(IFITQ13EI+4*J-3)
     &                    +LNTE1*WR1(IFITQ13EI+4*J-2)
     &                          +WR1(IFITQ13EI+4*J-1))
          Q14EI(J)= Q13EI(J)
          Q15EI(J)= Q13EI(J) 
C-----------------------------------------------------------------------
          WR2(IBINIJ+NE+J-1) = FAC2 /Q11EI(J)
          IF (60.D0*Q12EI(J)-48.D0*Q13EI(J)>=25.D0*Q11EI(J)) THEN
            WRITE(*,*) 'Incorrect electron-heavy interaction'
          ENDIF 
        ENDDO

C     2.3 Electron-ion interactions
        DO J = NN +1, NN +NI
          Q11EI(J) = AEQ11
          Q12EI(J) = AEQ12
          Q13EI(J) = AEQ13
          Q14EI(J) = AEQ14
          Q15EI(J) = AEQ15
          WR2(IBINIJ+NE+J-1) = FAC2 /AEQ11  
        ENDDO
       
C     2.4 Devoto collision integrals (divided by ne *n)
        WR2(IDQ00EE) = 0.D0; WR2(IDQ10EE) = 0.D0; WR2(IDQ20EE) = 0.D0
        WR2(IDQ11EE) = 0.D0; WR2(IDQ12EE) = 0.D0; WR2(IDQ22EE) = 0.D0
        DO I = 1, NN +NI
          FAC3 = 8.D0 *X(I+1) 
          WR2(IDQ00EE)     = WR2(IDQ00EE) +FAC3 *Q11EI(I)
          WR2(IDQ10EI+I-1) = -FAC3 /2.D0 *(5.D0*Q11EI(I)-6.D0*Q12EI(I))
          WR2(IDQ10EE)     = WR2(IDQ10EE) -WR2(IDQ10EI+I-1)
          WR2(IDQ20EI+I-1) = -FAC3 /8.D0 *
     &                      (35.D0*Q11EI(I)-84.D0*Q12EI(I)+48*Q13EI(I))
          WR2(IDQ20EE) = WR2(IDQ20EE) -WR2(IDQ20EI+I-1)
          WR2(IDQ11EE) = WR2(IDQ11EE) +FAC3 /4.D0 *(25.D0*Q11EI(I) 
     &    -60.D0*Q12EI(I) +48.D0*Q13EI(I))
          WR2(IDQ12EE) = WR2(IDQ12EE) +FAC3 /16.D0 *(175.D0*Q11EI(I) 
     &    -630.D0*Q12EI(I) +912.D0*Q13EI(I) -480.D0*Q14EI(I))
          WR2(IDQ22EE) = WR2(IDQ22EE) +FAC3 /64.D0 *(1225.D0*Q11EI(I) 
     &    -5880.D0*Q12EI(I) +12768.D0*Q13EI(I)-13440.D0*Q14EI(I) 
     &    +5760.D0*Q15EI(I))
        ENDDO
        FAC3 = 8.D0 *SQRT2 *X(1)
        WR2(IDQ11EE) = WR2(IDQ11EE) +FAC3 *Q22EE
        WR2(IDQ12EE) = WR2(IDQ12EE) +FAC3/4.D0*(7.D0*Q22EE -8.D0*Q23EE)
        WR2(IDQ22EE) = WR2(IDQ22EE) +FAC3/16.D0*(77.D0*Q22EE 
     &  -112.D0*Q23EE +80.D0*Q24EE)
      ENDIF

C     3. Heavy particle gas
C     3.1. Neutral-neutral and neutral-ion interactions      
      DO I = 1, NN
        J = I 
          MI = WR1(IMI+NE+I-1)
          FAC2 = FAC1 *DSQRT(2.D0 *TH /MI) 
          IJ = ((I-1)*(2*(NS-NE)-I)+2*J)/2
          Q11IJ = 1.D-20*DEXP( LNTH3*WR1(IFITQ11IJ+4*IJ-4)
     &                        +LNTH2*WR1(IFITQ11IJ+4*IJ-3)
     &                        +LNTH1*WR1(IFITQ11IJ+4*IJ-2)
     &                              +WR1(IFITQ11IJ+4*IJ-1))
C-----------------------------------------------------------------------
C       To be modified:
C          AIJ   = DEXP(FIT(LNTH1, WR1(IFITAIJ+4*IJ-4) ,4))
          AIJ = 1.D-20*DEXP( LNTH3*WR1(IFITAIJ+4*IJ-4)
     &                        +LNTH2*WR1(IFITAIJ+4*IJ-3)
     &                        +LNTH1*WR1(IFITAIJ+4*IJ-2)
     &                              +WR1(IFITAIJ+4*IJ-1)) /Q11IJ
C-----------------------------------------------------------------------
          WR2(IAIJ+IJ-1) = AIJ
          WR2(IBIJ+IJ-1) = DEXP( LNTH3*WR1(IFITBIJ+4*IJ-4)
     &                          +LNTH2*WR1(IFITBIJ+4*IJ-3)
     &                          +LNTH1*WR1(IFITBIJ+4*IJ-2)
     &                                +WR1(IFITBIJ+4*IJ-1)) 
          WR2(ICIJ+IJ-1) = DEXP( LNTH3*WR1(IFITCIJ+4*IJ-4)
     &                          +LNTH2*WR1(IFITCIJ+4*IJ-3)
     &                          +LNTH1*WR1(IFITCIJ+4*IJ-2)
     &                                +WR1(IFITCIJ+4*IJ-1))
          WR2(IETAI+I-1) = 5.D0/16.D0*DSQRT(WR1(IUPI) *WR1(IUKB) *MI
     &                     /WR1(IUNA) *TH)/(Q11IJ *AIJ)
          WR2(IBINIJ+NE*NS+IJ-1)= FAC2 /Q11IJ
C          IF (WR2(IBIJ+IJ-1)>=25.D0/12.D0) THEN 
C            WRITE(*,*) 'Incorrect heavy-heavy interaction'
C          ENDIF
        DO J = I+1, NN+NI
          MJ = WR1(IMI+NE+J-1)
          FAC2 = FAC1 *DSQRT(TH *(1.D0 /MI +1.D0 /MJ)) 
          IJ = ((I-1)*(2*(NS-NE)-I)+2*J)/2
          Q11IJ = 1.D-20 *DEXP( LNTH3*WR1(IFITQ11IJ+4*IJ-4)
     &                         +LNTH2*WR1(IFITQ11IJ+4*IJ-3)
     &                         +LNTH1*WR1(IFITQ11IJ+4*IJ-2)
     &                               +WR1(IFITQ11IJ+4*IJ-1)) 
C-----------------------------------------------------------------------
C       To be modified:
C          WR2(IAIJ+IJ-1) = DEXP(FIT(LNTH1, WR1(IFITAIJ+4*IJ-4) ,4))
          WR2(IAIJ+IJ-1) = 1.D-20 *DEXP( LNTH3*WR1(IFITAIJ+4*IJ-4)
     &                                  +LNTH2*WR1(IFITAIJ+4*IJ-3)
     &                                  +LNTH1*WR1(IFITAIJ+4*IJ-2)
     &                                  +WR1(IFITAIJ+4*IJ-1)) /Q11IJ 
C-----------------------------------------------------------------------
          WR2(IBIJ+IJ-1) = DEXP( LNTH3*WR1(IFITBIJ+4*IJ-4)
     &                          +LNTH2*WR1(IFITBIJ+4*IJ-3)
     &                          +LNTH1*WR1(IFITBIJ+4*IJ-2)
     &                                +WR1(IFITBIJ+4*IJ-1))
          WR2(ICIJ+IJ-1) = DEXP( LNTH3*WR1(IFITCIJ+4*IJ-4)
     &                          +LNTH2*WR1(IFITCIJ+4*IJ-3)
     &                          +LNTH1*WR1(IFITCIJ+4*IJ-2)
     &                                +WR1(IFITCIJ+4*IJ-1))
          WR2(IBINIJ+NE*NS+IJ-1)= FAC2 /Q11IJ
C          IF (WR2(IBIJ+IJ-1)>=25.D0/12.D0) THEN 
C            WRITE(*,*) 'Incorrect heavy-heavy interaction'
C          ENDIF
        ENDDO 
      ENDDO


C     3.2 Ion-Ion interactions (if NE /= 0)
      DO I = NN +1, NN+NI
        J = I
          MI = WR1(IMI+NE+I-1)
          FAC2 = FAC1 *DSQRT( 2.D0 *TH /MI ) 
          IJ = ((I-1)*(2*(NS-NE)-I)+2*J)/2
          WR2(IAIJ+IJ-1)  = RHA
          WR2(IBIJ+IJ-1)  = RHB
          WR2(ICIJ+IJ-1)  = RHC
          WR2(IETAI+I-1) = 5.D0/16.D0*DSQRT(WR1(IUPI) *WR1(IUKB) *MI
     &                     /WR1(IUNA) *TH)/(RHQ11 *RHA)
          WR2(IBINIJ+NE*NS+IJ-1)= FAC2 /RHQ11
C          IF (WR2(IBIJ+IJ-1)>=25.D0/12.D0) THEN 
C            WRITE(*,*) 'Incorrect heavy-heavy interaction'
C          ENDIF
        DO J = I+1, NN+NI
          MJ = WR1(IMI+NE+J-1)
          FAC2 = FAC1 *DSQRT(TH *(1.D0/MI +1.D0/MJ)) 
          IJ = ((I-1)*(2*(NS-NE)-I)+2*J)/2
          WR2(IAIJ+IJ-1) = RHA
          WR2(IBIJ+IJ-1) = RHB
          WR2(ICIJ+IJ-1) = RHC
          WR2(IBINIJ+NE*NS+IJ-1)= FAC2 /RHQ11
C          IF (WR2(IBIJ+IJ-1)>=25.D0/12.D0) THEN 
C            WRITE(*,*) 'Incorrect heavy-heavy interaction'
C          ENDIF
        ENDDO 
      ENDDO

      END SUBROUTINE COLLISION 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE MEANFP (WR1, LWR1, WR2, LWR2, T, ND, X, MFP)
C-----------------------------------------------------------------------
C     This subroutine computes the mean free path. 
C     The average cross section is considered and not the transport 
C     collision integral. Then, pi sigma^2 Omega* is used.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LWR2
      DOUBLE PRECISION WR1(1:LWR1), WR2(1:LWR2), X(1:NS), T, ND, MFP,
     &                 Q(1:NS,1:NS), MI, MJ
C-----------------------------------------------------------------------
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      DO I = 1, NS
        MI = WR1(IMI+I-1)
        DO J = I, NS
          MJ = WR1(IMI+J-1)
          FAC = SQRT(2.D0*WR1(IUPI)*WR1(IUR)*T*(MI+MJ)/(MI*MJ))
          IJ = ((I-1)*(2*NS-I)+2*J)/2
          Q(I,J)   = 3.D0 *FAC /(16.D0 *WR2(IBINIJ+IJ-1))
          IF (I.NE.J) THEN
            Q(J,I)   = Q(I,J)
          ENDIF
        ENDDO
      ENDDO

      MFP = 0.D0
      DO I = 1, NS
        CS = 0.D0
        DO J = 1, NS
          CS = CS + X(J) *Q(I,J)
        ENDDO
        MFP = MFP + X(I) *CS
      ENDDO
      MFP = 1.D0 /(MFP *ND)

      END SUBROUTINE 
C-----------------------------------------------------------------------
