
C     MUTATION library version 1.0
C     von Karman Institute, Rhode-Saint-Genese (Belgium), 09-30-2002

C     MUTATION library version 1.1
C     Ecole Centrale Paris, Chatenay-Malabry (France), 05-31-2005

C     MUTATION library version 1.2
C     Stanford University, Palo Alto (USA), 12-31-2006

C     Thierry E. Magin

C-----------------------------------------------------------------------
      PROGRAM MUTATION
C-----------------------------------------------------------------------
C     This program computes the equilibrium composition and MUlti-
C     component Thermodynamic And Transport properties of an IONized 
C     mixture.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      INTEGER LWR1, LWR2, LWR3, LWR4, LWI, LWC, LWR1MAX, LWR2MAX, 
     &        LWR3MAX, LWR4MAX, LWIMAX, LWCMAX, NMAX, OUT,NMIXMAX,
     &        NSPCMAX 
      PARAMETER (LWR1MAX = 100000, LWR2MAX = 20000, LWR3MAX = 1000,
     &           LWR4MAX = 1000000, LWIMAX = 10000, LWCMAX = 500,
     &           NMAX = 400, OUT = 50, NMIXMAX= 100, NSPCMAX = 40)
      INTEGER WI(1:LWIMAX), IS, J, IT, NT, IP, NP, RULE, IC1,
     &        MIX(1:NMIXMAX), NMIX, SPC(1:NSPCMAX), NSPC, 
     &        IHSONINE, IESONINE, IMOD
      DOUBLE PRECISION WR1(1:LWR1MAX), WR2(1:LWR2MAX), WR3(1:LWR3MAX),
     &                 WR4(1:LWR4MAX),
     &                 TMIN, TMAX, DELTAT, PMIN, PMAX, DELTAP, P, TH,
     &                 TE, TR, TV(1:NMAX), ND, RHO, TOL, EPS, EPSP1,
     &                 XN(1:NMAX), XINI(1:NMAX), X(1:NMAX), XP(1:NMAX),
     &                 XTOL(1:NMAX), YINI(1:NMAX), Y(1:NMAX), MM,
     &                 MIXH1, MIXH2, MIXH3, MIXH4, MIXH5, MIXH6,
     &                 MIXH1P, MIXH2P, MIXH3P, MIXH4P, MIXH5P, MIXH6P,
     &                 MIXS1, MIXS2, MIXS3, MIXS4, MIXS5, MIXS6, 
     &                 S1(1:NMAX), S2(1:NMAX), S3(1:NMAX), S4(1:NMAX), 
     &                 S5(1:NMAX), S6(1:NMAX), MIXCP, MIXFCP, 
     &                 DF(1:NMAX), HM1(1:NMAX), HM2(1:NMAX), 
     &                 HM3(1:NMAX), HM4(1:NMAX), HM5(1:NMAX), 
     &                 HM6(1:NMAX), H1(1:NMAX), H2(1:NMAX), H3(1:NMAX),
     &                 H4(1:NMAX), H5(1:NMAX), H6(1:NMAX), H1P(1:NMAX), 
     &                 H2P(1:NMAX), H3P(1:NMAX), H4P(1:NMAX), 
     &                 H5P(1:NMAX), NDP, RHOP, H6P(1:NMAX), CPR(1:NMAX),
     &                 CPV(1:NMAX), CPE(1:NMAX), CPINT(1:NMAX), U, 
     &                 ONE(1:NMAX), MMP, GAMMAF, SOUNDF, GAMMAE, DRHODP,
     &                 SOUNDE, MACH, GAMMAI, TT, PT, PTB, PTSV, PTIS
      DOUBLE PRECISION MIXMCP, ETA, LAMBDATH, CHIH(1:NMAX), LAMBDATE, 
     &                 LAMBDATD, LAMBDATDH, LAMBDATDE, LAMBDATOT,
     &                 CHIE(1:NMAX), LAMBDAINTE, LAMBDAINTR, LAMBDAINTV,
     &                 LAMBDAINTH, LAMBDAREA, LAMBDASM, LAMBDAR, 
     &                 LAMBDAK, JDIF(1:NMAX), JDIFD(1:NMAX), LAMBDAFI, 
     &                 JDIFR(1:NMAX), JDIFK(1:NMAX), JDIFF(1:NMAX),
     &                 JDIFFR(1:NMAX), LAMBDAFR,
     &                 FIJ(1:NMAX*(NMAX-1)/2), SIGMA, SIGMAI, 
     &                 EAMB, EAMBR, EAMBF, EAMBK, DUM, LAMBDAI(1:NMAX), 
     &                 ETAI(1:NMAX), MFP, ARRAY(1:NMIXMAX+NMAX*NSPCMAX),
     &                 THEPSP1, TEEPSP1, TREPSP1, TVEPSP1(1:NMAX) 
      CHARACTER WC(1:LWCMAX)
      CHARACTER(10) MIXTURE, REACTION, TRANSFER
      CHARACTER(100) PATH      
      EXTERNAL ENTHALPY, MINRES, CG
C-----------------------------------------------------------------------
      DUM = 0.D0

C     Presentation screen
      CALL PRESENTATION (1)

C     Read inputfile
      PATH = '..'
      LPATH = LCHAR(PATH)
      CALL READINPUT (RULE, PMIN, PMAX, DELTAP, TMIN, TMAX, DELTAT, TOL,
     &                PATH, MIXTURE, U, NMIXMAX, NMIX, MIX, NSPCMAX, 
     &                NSPC, SPC, IHSONINE, IESONINE, IMOD)

      NP = 1; NT = 1
      IF (DELTAP /= 0.D0) NP = 1 +INT((PMAX -PMIN) /DELTAP)
      IF (DELTAT /= 0.D0) NT = 1 +INT((TMAX -TMIN) /DELTAT)

C     LTE, no finite rate chemistry      
      REACTION = 'empty'
      TRANSFER = 'empty'

C     Epsilon for finite difference (specific heat)
      EPS = 1.D-2; EPSP1 = 1.D0 +EPS

C     Pseudo-dynamic allocation
      CALL LENGTH (PATH, MIXTURE, REACTION, TRANSFER, LWR1, LWR2, LWR3, 
     &             LWR4, LWI, LWC, NS, NE, NC, NREA, NV, NMAX, NVIB,
     &             NELEQ)
      IF (     ( LWR1 > LWR1MAX ) .OR. ( LWR2 > LWR2MAX ) 
     &    .OR. ( LWR3 > LWR3MAX ) .OR. ( LWR4 > LWR4MAX )
     &    .OR. ( LWI  > LWIMAX  ) .OR. ( LWC  > LWCMAX  ) 
     &    .OR. ( NS   > NMAX    ) .OR. ( NMIX > NMIXMAX )
     &    .OR. ( NSPC > NSPCMAX ) .OR. ( NV > NMAX      ) ) THEN
         WRITE(*,*) LWR1, LWR2, LWR3, LWR4, LWI, LWC, NS, NMIX, NSPC, NV
         WRITE(*,*) 'Error: work arrays out of space.'
         PAUSE
      ENDIF 

C     Initialize data independent on temperature and pressure
      CALL INITIALIZE (PATH, MIXTURE, REACTION, TRANSFER, WR1, LWR1, 
     &                 WR3, LWR3, WR4, LWR4, WI, LWI, WC, LWC, IMOD)

C     Initial composition and nuclear fraction ratios
      DO IS = 1, NS
        X(IS) = 1.D0; Y(IS) = 1.D0
      ENDDO

      CALL NUCLEAR (WR1, LWR1, XN)

C     Open output file
      OPEN(UNIT=OUT, FILE=PATH(1:LPATH)//'/output/output.dat',
     &     STATUS='UNKNOWN')

C     Pressure
      DO IP = 1, NP

      P = PMIN + (IP -1) *(DELTAP)
      WRITE(*,*) '      -Pressure: ',P,' [Pa]'

C     Temperature 
      DO IT = 1, NT

      TH      = TMIN + (IT -1) *DELTAT 
C      TH      = TMIN + (NT-IT) *DELTAT 
      THEPSP1 = TH*EPSP1
      WRITE(*,*) '       -Temperature: ',TH,' [K]'
      TE = TH; TEEPSP1 = TE*EPSP1
      TR = TH; TREPSP1 = TR*EPSP1
      IF (NV==0) THEN
        TV(1)      = 1.D0
        TVEPSP1(1) = TV(1)*EPSP1
      ELSE
        DO J = 1, NV
          TV(J)      = TH
          TVEPSP1(J) = TV(J)*EPSP1 
        ENDDO
      ENDIF

C     Equilibrium composition
        DO IS = 1, NS
          XINI(IS) = X(IS); YINI(IS) = Y(IS)
        ENDDO
        CALL COMPOSITION (WR1, LWR1, WI, LWI, TH, P, XN, XINI, X)
        CALL COMPOSITION (WR1, LWR1, WI, LWI, THEPSP1, P, XN, X, XP)

        CALL NUMBERD (WR1, LWR1, P, TH, TE, X, ND)
        CALL NUMBERD (WR1, LWR1, P, THEPSP1, TEEPSP1, XP, NDP)
        CALL DENSITY (WR1, LWR1, X, ND, RHO)
        CALL DENSITY (WR1, LWR1, XP, NDP, RHOP)

        CALL MOLARMASS (WR1, LWR1, RHO, ND, MM)
        CALL MOLARMASS (WR1, LWR1, RHOP, NDP, MMP)

        CALL MASSCOMPOSITION (WR1, LWR1, WI, LWI, TH, RHO, XN, YINI, Y)

C     Thermodynamic properties
C     -Species enthalpy per unit mass
        CALL ENTHALPYMASS (WR1, LWR1, WI, LWI, TH, TE, TR, TV, P, HM1, 
     &                     HM2, HM3, HM4, HM5, HM6)

C     -Mixture enthalpy per unit mole
        CALL MIXPROPERTY (WR1, LWR1, WI, LWI, TH, TE, TR, TV, DUM, X,
     &                    MIXH1, MIXH2, MIXH3, MIXH4, MIXH5, MIXH6, 
     &                    ENTHALPY)

C     -Mixture equilibrium specific heat per unit mole
        CALL MIXPROPERTY (WR1, LWR1, WI, LWI,THEPSP1, TEEPSP1, TREPSP1, 
     &                    TVEPSP1, DUM, XP, MIXH1P, MIXH2P, MIXH3P, 
     &                    MIXH4P, MIXH5P, MIXH6P, ENTHALPY)
        MIXCP = (MIXH1P -MIXH1) /(TH *EPS)

C     -Mixture equilibrium specific heat per unit mass
        MIXMCP = (MIXH1P/MMP -MIXH1/MM) /(TH *EPS)

C     -Mixture frozen specific heat per unit mole
        CALL MIXPROPERTY (WR1, LWR1, WI, LWI, THEPSP1, TEEPSP1, TREPSP1,
     &                    TVEPSP1, DUM, X, MIXH1P, MIXH2P, MIXH3P, 
     &                    MIXH4P, MIXH5P, MIXH6P, ENTHALPY)
        MIXFCP =(MIXH1P -MIXH1) /(TH *EPS)

C     -Mixture entropy per unit mole
        CALL MIXENTROPY (WR1, LWR1, WI, LWI, TH, TE, TR, TV, P, X,
     &                   MIXS1, MIXS2, MIXS3, MIXS4, MIXS5, MIXS6)

C     -Species specific entropy per unit mole
        CALL ENTROPY (WR1, LWR1, WI, LWI, TH, TE, TR, TV, P, S1, 
     &                S2, S3, S4, S5, S6)

C     -Species specific heat per unit mole
        CALL ENTHALPY (WR1, LWR1, WI, LWI, TH, TE, TR, TV, P, H1, H2,
     &                 H3, H4, H5, H6)

        CALL ENTHALPY (WR1, LWR1, WI, LWI, THEPSP1, TEEPSP1, TREPSP1,
     &                 TVEPSP1, P, H1P, H2P, H3P, H4P, H5P, H6P)
        DO IS = 1, NS
          CPE(IS)   = (H3P(IS) -H3(IS)) / (EPS *TE)
          CPR(IS)   = (H4P(IS) -H4(IS)) / (EPS *TR)
          CPV(IS)   = (H5P(IS) -H5(IS)) / (EPS *TV(1))
          CPINT(IS) = CPE(IS) +CPR(IS) +CPV(IS)
        ENDDO     

C     -Frozen speed of sound
        CALL FROZENGAMMA (WR1, LWR1, WI, LWI, TH, P, X, EPS, GAMMAF)
        SOUNDF = DSQRT(GAMMAF *P /RHO)

C     -Equilibrium speed of sound
        CALL EQUIGAMMA (WR1, LWR1, WI, LWI, TH, P, RHO, XN, X, EPS, 
     &                  GAMMAE, DRHODP)
        SOUNDE = DSQRT(GAMMAE /DRHODP)
        GAMMAI = RHO /P *GAMMAE /DRHODP
        MACH   = U /SOUNDE

C     -Total pressure and temperature
        CALL BERNOULLI (P, TH, RHO, U, PTB)
        CALL SAINTVENANT (P, MACH, GAMMAE, PTSV)
        CALL SAINTVENANT (P, MACH, GAMMAI, PTIS)
        CALL TOTAL (WR1, LWR1, WI, LWI, TH, P, U, XN, X, EPS, TT, 
     &              PT)

        IF (IMOD==1) THEN
C     Transport properties
C     BEWARE! A small number is added to the mole fractions to compute 
C             the transport properties, respecting the mass 
C             conservation. Don't used the modified mole fractions as 
C             initial guess to compute the new composition, otherwise
C             the Newton convergence will be killed.
          CALL COMPOTOL (X, TOL, XTOL)
C     Update of kinetic data
          CALL COLLISION (WR1, LWR1, WR2, LWR2, TH, TE, ND, X)
          CALL MEANFP (WR1, LWR1, WR2, LWR2, TH, ND, X, MFP)

C     A. Heavy particle mixture properties
          DO IS = 1, NS-NE
            ETAI(IS+NE)    = WR2(IETAI+IS-1)
            LAMBDAI(IS+NE) = 3.75D0 *WR1(IUR) /WR1(IMI+NE+IS-1) 
     &                      *WR2(IETAI+IS-1)
          ENDDO
C     Eucken corrections for internal energy (no inelastic collisions)
          CALL CORRECTION (WR1, LWR1, WR2, LWR2, XTOL, IHSONINE, FIJ)
          CALL LAMBDAINT (WR1, LWR1, WR2, LWR2, CPE, XTOL, LAMBDAINTE)
          CALL LAMBDAINT (WR1, LWR1, WR2, LWR2, CPR, XTOL, LAMBDAINTR)
          CALL LAMBDAINT (WR1, LWR1, WR2, LWR2, CPV, XTOL, LAMBDAINTV)
          CALL LAMBDAINT (WR1, LWR1, WR2, LWR2, CPINT, XTOL, LAMBDAINTH)

          SELECT CASE(RULE)
C     A.1. Direct method
            CASE(1,5,6,7)
              CALL ETAD (WR1, LWR1, WR2, LWR2, XTOL, ETA)
 
              CALL LAMBDACHID (WR1, LWR1, WR2, LWR2, XTOL, LAMBDATH,
     &                         CHIH)
              CALL LAMBDAREAD (WI, LWI, WR1, LWR1, WR2, LWR2, H1, TH,
     &                         XTOL, LAMBDAREA)

C     A.2. Preconditioned conjugate gradient
            CASE(2)
              CALL ETACG (WR1, LWR1, WR2, LWR2, XTOL, ETA)

              CALL LAMBDACHICG (WR1, LWR1, WR2, LWR2, XTOL, LAMBDATH, 
     &                         CHIH)
              CALL LAMBDAREACG (WI, LWI, WR1, LWR1, WR2, LWR2, H1, TH, 
     &                          XTOL, LAMBDAREA)

C     A.3. Yos' rule
            CASE(3)
              CALL ETAYOS (WR1, LWR1, WR2, LWR2, XTOL, ETA)

              CALL LAMBDACHID (WR1, LWR1, WR2, LWR2, XTOL, LAMBDATH, 
     &                         CHIH)
              CALL LAMBDAYOS (WR1, LWR1, WR2, LWR2, XTOL, LAMBDATH)
              CALL LAMBDAREAYOS (WI, LWI, WR1, LWR1, WR2, LWR2, H1, TH,
     &                           XTOL, LAMBDAREA)

C     A.4. Wilke's rule
            CASE(4)
              CALL ETAWILKE (WR1, LWR1, WR2, LWR2, XTOL, ETA)

              CALL LAMBDACHID (WR1, LWR1, WR2, LWR2, XTOL, LAMBDATH, 
     &                         CHIH)
              CALL LAMBDAWILKE (WR1, LWR1, WR2, LWR2, XTOL, LAMBDATH)
              CALL LAMBDAREAD (WI, LWI, WR1, LWR1, WR2, LWR2, H1, TH,
     &                         XTOL, LAMBDAREA)

          END SELECT

C     B. Electron gas properties
          IF (NE /= 0) THEN
            CALL PHIE (WR1, LWR1, WR2, LWR2, XTOL, TE, FIJ, IESONINE)
            CALL TDIFE (WR1, LWR1, WR2, LWR2, XTOL, TE, CHIE, 3)
            CALL LAMBDAE (WR1, LWR1, WR2, LWR2, X, TE, LAMBDATE, 3)
            CALL SIGMAE (WR1, LWR1, WR2, LWR2, X, TE, SIGMA, 2)
            ETAI(1)    = 5.D0 /16.D0 *DSQRT(WR1(IUPI) *WR1(IUKB) 
     &                   *WR1(IMI) /WR1(IUNA) *TE) /WR2(IQEE22)
            LAMBDAI(1) = 3.75D0 *WR1(IUR) /WR1(IMI+1) *ETAI(1)
          ELSE
            LAMBDATE = 0.D0
            SIGMA = 0.D0
            DO IS = 1, NS
              CHIE(IS) = 0.D0
            ENDDO
          ENDIF

C     C. Mixed heavy and electron gas properties
C     -Electrical conductivity (ion and electron contributions)
          DO IS = 1, NS*(NS-1)/2
            IF (FIJ(IS)<=-1) THEN
              WRITE(*,*) 'Incorrect correction function for SM'
            ENDIF 
          ENDDO
          DO IS = 1, NS
            DF(IS) = -(X(IS) *WR1(IUE) *WI(IQI+IS-1)) /(WR1(IUKB) *TH) 
          ENDDO
          CALL SMNEUTD (WR1, LWR1, WR2, LWR2, XTOL, ND, DF, FIJ, JDIF)
          SIGMAI = 0.D0
          DO IS = NE+1, NS
            SIGMAI = SIGMAI +JDIF(IS) *WR1(IUE) *WI(IQI+IS-1)
     &               /WR1(IMI+IS-1) *WR1(IUNA)
          ENDDO         

C     -Stefan-Maxwell equations (1rst order Sonine expansion), the 
C      driving forces are concentration gradients due to unity 
C      temperature gradient (constant pressure).
          DO IS = 1, NS
C            DF(IS) = (XP(IS) -X(IS)) /(TH *EPS) +(CHIE(IS) +CHIH(IS)) 
C     &               /TH
C TEMPORARY Uncomment the following line to neglect thermal diffusion as 
C           driving force
           DF(IS) = (XP(IS) -X(IS)) /(TH *EPS) 
C TEMPORARY
          ENDDO        
C     Direct method used for thermal diffusion conductivity
          IF (NE /= 0) THEN
            CALL SMDSPD (WR1, LWR1, WR2, LWR2, WI, LWI, XTOL, TH, TE,
     &                   ND, DF, FIJ, JDIFD, EAMB)
            CALL SMD (WR1, LWR1, WR2, LWR2, WI, LWI, XTOL, TH, TE, ND, 
     &                DF, FIJ, JDIFD, EAMB)
            CALL SMRAMD (WR1, LWR1, WR2, LWR2, WI, LWI, XTOL, TH, TE, 
     &                   ND, DF, FIJ, JDIFR, EAMBR)
            CALL SMKOL (WR1, LWR1, WR2, LWR2, WI, LWI, XTOL, TH, TE, ND,
     &                  DF, FIJ, JDIFK, EAMBK)
            CALL SMFICK (WR1, LWR1, WR2, LWR2, WI, LWI, XTOL, TH, TE, 
     &                   ND, DF, FIJ, JDIFF, JDIFFR, EAMBF)
          ELSE
            CALL SMNEUTD (WR1, LWR1, WR2, LWR2, XTOL, ND, DF, FIJ, 
     &                    JDIFD)
            EAMB  = 0.D0; EAMBR = 0.D0; EAMBK = 0.D0
          ENDIF

          SELECT CASE(RULE)

C     C.1. Direct method
            CASE(1,3,4)
            DO IS = 1, NS
              JDIF(IS) = JDIFD(IS) 
            ENDDO        

C     C.2. Preconditioned Krylov method
            CASE(2)
            IF (NE /= 0) THEN
              CALL SMRAMCG (WR1, LWR1, WR2, LWR2, WI, LWI, XTOL, TH, TE,
     &                      ND, DF, FIJ, JDIFR, EAMBR)
              CALL SMGMRES (WR1, LWR1, WR2, LWR2, WI, LWI, XTOL, TH, TE,
     &                      ND, DF, FIJ, JDIF, EAMB)
              CALL SMDSPCG (WR1, LWR1, WR2, LWR2, WI, LWI, XTOL, TH, TE,
     &                      ND, DF, FIJ, JDIF, EAMB)

            ELSE
              CALL SMNEUTCG (WR1, LWR1, WR2, LWR2, XTOL, ND, DF, FIJ, 
     &                       JDIF)
            ENDIF
C     C.3. Sutton's algorithm 
            CASE(5)
            IF (NE /= 0) THEN
              CALL SMRAMSUT (WR1, LWR1, WR2, LWR2, WI, LWI, XTOL, TH, 
     &                       TE, ND, DF, FIJ, JDIFR, EAMBR)
              CALL SMSUT (WR1, LWR1, WR2, LWR2, WI, LWI, XTOL, TH, TE,
     &                    ND, DF, FIJ, JDIF, EAMB)
            ELSE
              CALL SMNEUTSUT (WR1, LWR1, WR2, LWR2, XTOL, ND, DF, FIJ, 
     &                        JDIF)
            ENDIF
          END SELECT
 
          LAMBDASM  = 0.D0; LAMBDAR = 0.D0; LAMBDATDH = 0.D0
          LAMBDATDE = 0.D0; LAMBDAK = 0.D0; LAMBDAFI  = 0.D0
          LAMBDAFR  = 0.D0 
          DO IS = 1, NS
            LAMBDASM  =  LAMBDASM -HM1(IS)* JDIF(IS)
            LAMBDAR   =  LAMBDAR  -HM1(IS)* JDIFR(IS)
            LAMBDAK   =  LAMBDAK  -HM1(IS)* JDIFK(IS)
            LAMBDAFI  =  LAMBDAFI -HM1(IS)* JDIFF(IS)
            LAMBDAFR  =  LAMBDAFR -HM1(IS)* JDIFFR(IS)
            LAMBDATDH =  LAMBDATDH -CHIH(IS)   
     &                   /(WR1(IMI+IS-1) *XTOL(IS)) *JDIFD(IS)
            LAMBDATDE =  LAMBDATDE -CHIE(IS)
     &                   /(WR1(IMI+IS-1) *XTOL(IS)) *JDIFD(IS)
          ENDDO
          LAMBDATDH = LAMBDATDH *P *WR1(IUNA) /ND
          LAMBDATDE = LAMBDATDE *P *WR1(IUNA) /ND
          LAMBDATD  = LAMBDATDH +LAMBDATDE
          LAMBDATOT = LAMBDATD +LAMBDAINTH +LAMBDATH +LAMBDATE +LAMBDASM
        ENDIF

C     -Species specific heat per unit mole
        CALL FILLARRAY (NMIX, MIX, NS, NSPC, SPC, ARRAY,
     &                 TH, P, RHO, ND, MM, TT, PT, PTB, PTSV, PTIS,
     &                 MIXH1, MIXH2, MIXH3, MIXH4, MIXH5, MIXH6,
     &                 MIXS1, MIXS2, MIXS3, MIXS4, MIXS5, MIXS6, 
     &                 MIXCP, MIXFCP, CPE, CPR, CPV, CPINT,
     &                 GAMMAF, GAMMAE, GAMMAI, DRHODP,
     &                 SOUNDF, SOUNDE, MACH, MIXMCP, 
     &                 ETA, LAMBDATH, LAMBDATE,
     &                 LAMBDAINTE, LAMBDAINTR, LAMBDAINTV, LAMBDAINTH, 
     &                 LAMBDAREA, LAMBDASM, LAMBDAR, LAMBDAFI, LAMBDAFR,
     &                 LAMBDAK, 
     &                 LAMBDATDH, LAMBDATDE, LAMBDATD, LAMBDATOT,
     &                 ETAI, LAMBDAI, SIGMA, SIGMAI, MFP, 
     &                 X, Y, CHIH, CHIE, JDIF, JDIFR, JDIFK, DF, EAMB, 
     &                 EAMBR, EAMBK, H1, H2, H3, H4, H5, H6, S1, S2, 
     &                 S3, S4, S5, WR1(IMI))

        WRITE (OUT,2010) (ARRAY(I),I=1,NMIX+NS*NSPC)

      ENDDO
      ENDDO

C     Output format
2010  FORMAT (100E14.6)

C     Close output file
      CLOSE(OUT)

C     Presentation screen
      CALL PRESENTATION (2)

      END PROGRAM MUTATION
C-----------------------------------------------------------------------
