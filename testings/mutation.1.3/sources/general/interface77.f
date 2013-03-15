
C     Interface for the mutation library version 1.3
C     Stanford University, Palo Alto (USA), 02-18-2009

C     Thierry E. Magin

C-----------------------------------------------------------------------
      PROGRAM INTERFACE77 
C-----------------------------------------------------------------------
C     This program is an interface in FORTRAN G77 for the mutation 
C     library
C-----------------------------------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------------------------------
C     Internal variables
      INTEGER LWR1MAX, LWR2MAX, LWR3MAX, LWR4MAX, LWIMAX, LWCMAX, NMAX,
     &        LCHAR 
      PARAMETER (LWR1MAX = 4000, LWR2MAX = 1000, LWR3MAX = 1000,
     &           LWR4MAX = 5000, LWIMAX = 1000, LWCMAX = 100, NMAX = 35)
      INTEGER IS, JS, IT, NT, IP, NP, IV, IVIB, LPATH
      DOUBLE PRECISION TMIN, TMAX, DELTAT, PMIN, PMAX, DELTAP, U,
     &                 EPS, EPSP1, RHOI2(1:NMAX), RHO2
C     Variables passed as arguments to the library:
C     -Mandatory routines for thermo and chemistry (LENGTH, INITIALIZE)
      INTEGER LWR1, LWR3, LWR4, LWI, LWC
      DOUBLE PRECISION WR1(1:LWR1MAX), WR3(1:LWR3MAX), WR4(1:LWR4MAX) 
      INTEGER WI(1:LWIMAX)
      CHARACTER WC(1:LWCMAX)
      INTEGER NS, NE, NC, NREA, NV, NVIB, NELEQ, IMOD
      CHARACTER(10) MIXTURE, REACTION, TRANSFER
      CHARACTER(100) PATH 
C     -Mandatory routines for transport (COLLISION) requiring as well 
C      the mandatory thermo and chemistry routines (LENGTH, INITIALIZE)
      INTEGER LWR2
      DOUBLE PRECISION WR2(1:LWR2MAX), TH, TE, ND, X(1:NMAX)
C     -Optional routines
      INTEGER  IHSONINE, IESONINE
      DOUBLE PRECISION P, TR, TV(1:NMAX), THP, TEP, TRP,
     &                 TVP(1:NMAX), MI(1:NMAX), XP(1:NMAX),
     &                 XINI(1:NMAX), YINI(1:NMAX), XN(1:NMAX),
     &                 MM, MMP, HM1(1:NMAX), HM2(1:NMAX),  HM3(1:NMAX), 
     &                 HM4(1:NMAX), HM5(1:NMAX), HM6(1:NMAX), 
     &                 RHOI(1:NMAX), P2, T2, T3, X2(1:NMAX), Y(1:NMAX),
     &                 EM1(1:NMAX), EM2(1:NMAX), EM3(1:NMAX), 
     &                 EM4(1:NMAX), EM5(1:NMAX), EM6(1:NMAX), 
     &                 EM1P(1:NMAX), EM2P(1:NMAX), EM3P(1:NMAX), 
     &                 EM4P(1:NMAX), EM5P(1:NMAX), EM6P(1:NMAX), 
     &                 TOLX, TOLY, TVEC(1:NMAX), OMEGA(1:NMAX),
     &                 RHOE, TINI, RHO, DUM, XTOL(1:NMAX), ETA, 
     &                 Y3(1:NMAX), Y4(1:NMAX), XN2(1:NMAX), LAMBDATE,
     &                 CPR, CPV, CPE, CPINT(1:NMAX), EAMB, DF(1:NMAX),
     &                 CHIE(1:NMAX), LAMBDATH, CHIH(1:NMAX),
     &                 FIJ(1:NMAX*(NMAX-1)/2), JDIF(1:NMAX), 
     &                 H1(1:NMAX), H2(1:NMAX), H3(1:NMAX),
     &                 H4(1:NMAX), H5(1:NMAX), H6(1:NMAX), LAMBDATOT,
     &                 LAMBDAR, LAMBDAINTH, LAMBDATOT2, Y5(1:NMAX), 
     &                 LAMBDATOT3, DIJ(1:NMAX,1:NMAX), JDIF2(1:NMAX), 
     &                 LAMBDAR2, XN1(1:NMAX)
C-----------------------------------------------------------------------
C     Tolerance for mole fractions (avoid singularity in the transport 
C     properties)
      TOLX = 1.D-16

C     Tolerance for mass fractions (avoid singularity in the Arrhenius 
C     rates)
      TOLY = 1.D-20

C     Epsilon for finite difference (specific heat, etc.)
      EPS = 1.D-2; EPSP1 = 1.D0 +EPS

C     Dummy variable
      DUM = 0.D0

C     Velocity for Mach number, etc. 
C     (quantity usually provided by the CFD code)
      U = 0.D0

      PMIN = 101325.D0; PMAX = 101325.D0; DELTAP = 0.D0
      TMIN = 500.D0;   TMAX = 15000.D0;  DELTAT = 500.D0
C      PMIN = 5491500.D0; PMAX = PMIN; DELTAP = 0.D0
C      TMIN = 3529.03D0;   TMAX = TMIN;  DELTAT = 0.D0
      NP = 1; NT = 1
      IF (DELTAP /= 0.D0) NP = 1 +INT((PMAX -PMIN) /DELTAP)
      IF (DELTAT /= 0.D0) NT = 1 +INT((TMAX -TMIN) /DELTAT)
C-----------------------------------------------------------------------
C     Type of routines used by the CFD code
C      flag IMOD = 0: only thermodynamic and chemistry routines
C                = 1: thermodynamic, chemistry and transport routines
      IMOD = 1 
C-----------------------------------------------------------------------
C     Files and path
C      -mixture file located in directory /data/mixture
      MIXTURE = 'air5'
C      -reaction file located in directory /data/chemistry/gasreact
C       for nonreacting mixture, use 'empty' 
      REACTION = 'air5'
C      -energy transfer file located in directory /data/transfer
      TRANSFER = 'empty'
C      -path for the library location
      PATH = '..'
      LPATH = LCHAR(PATH)
C-----------------------------------------------------------------------
C     Pseudo-dynamic allocation (mandatory routine)
      CALL LENGTH (PATH, MIXTURE, REACTION, TRANSFER, LWR1, LWR2, LWR3, 
     &             LWR4, LWI, LWC, NS, NE, NC, NREA, NV, NMAX, NVIB,
     &             NELEQ)
       write(*,*) NS
C      => output variables:
C     -LWR1:  length of WR1 vector
C     -LWR2:  length of WR2 vector
C     -LWR3:  length of WR3 vector
C     -LWR4:  length of WR4 vector
C     -LWI:   length of WI vector
C     -LWC:   length of WC vector
C     -NS:    number of species
C     -NE:    number of electrons (0/1)
C     -NC:    number of elements (including charge for ionized mixtures)
C     -NREA:  number of chemical reactions
C     -NV:    number of vibrational modes (cumulated over all molecules)
C     -NVIB:  number of distinct vibrational temperatures (0/1/2...)
C     -NELEQ: number of distinct electron/electronic temperatures (0/1)
C     Attention: variable NMAX will disappear!
C-----------------------------------------------------------------------
      IF (     ( LWR1 > LWR1MAX ) .OR. ( LWR2 > LWR2MAX ) 
     &    .OR. ( LWR3 > LWR3MAX ) .OR. ( LWR4 > LWR4MAX )
     &    .OR. ( LWI  > LWIMAX  ) .OR. ( LWC  > LWCMAX  ) 
     &    .OR. ( NS   > NMAX    ) .OR. ( NV > NMAX      ) ) THEN
         WRITE(*,*) LWR1, LWR2, LWR3, LWR4, LWI, LWC, NS,  NV
         WRITE(*,*) 'Error: work arrays out of space.'
         PAUSE
      ENDIF 
C-----------------------------------------------------------------------
C     Initialize data independent on temperature and pressure
C     The routine INITIALIZE is called only once (mandatory)
      CALL INITIALIZE (PATH, MIXTURE, REACTION, TRANSFER, WR1, LWR1, 
     &                 WR3, LWR3, WR4, LWR4, WI, LWI, WC, LWC, IMOD)
C      => output variables:
C     -WR1
C     -WR3
C     -WR4
C     -WI
C     -WC
C-----------------------------------------------------------------------
C     Initial composition fraction ratios (a good guess can be obtained 
C     by the composition at the previous iteration) 
C     -X(1:NS):  mole fractions
C     -Y(1:NS):  mass fractions
      DO IS = 1, NS
        X(IS) = 1.D0/NS; Y(IS) = 1.D0/NS
      ENDDO
C-----------------------------------------------------------------------
C     Element fractions for equilibrium composition calculations
C     read from data/chemistry/ceq/file.ceq
      CALL NUCLEAR (WR1, LWR1, XN)
C      => output variables:
C     -XN(1:NC): element fractions
C-----------------------------------------------------------------------
C     Open output file
      OPEN(UNIT=99, FILE=PATH(1:LPATH)//'/output/output.dat',
     &     STATUS='UNKNOWN')

C     Pressure
      DO IP = 1, NP

      P = PMIN + (IP -1) *(DELTAP)
      WRITE(*,*) '      -Pressure: ',P,' [Pa]'

C     Temperature
      DO IT = 1, NT

      TH      = TMIN + (IT -1) *DELTAT
C      TH      = TMIN + (NT-IT) *DELTAT
      THP = TH*EPSP1
      WRITE(*,*) '       -Temperature: ',TH,' [K]'
      TE = TH; TEP = TE*EPSP1
      TR = TH; TRP = TR*EPSP1
      IF (NV==0) THEN
        TV(1)      = 1.D0
        TVP(1)     = TV(1)*EPSP1
      ELSE
        DO IV = 1, NV
          TV(IV)      = TH
          TVP(IV)     = TV(IV)*EPSP1 
        ENDDO
      ENDIF
      TVEC(1)      = TH
      DO IVIB = 1, NVIB
        TVEC(1+IVIB)  = TV(1)
      ENDDO
      TVEC(NVIB+2) = TE
C-----------------------------------------------------------------------
C     Equilibrium composition (in thermal equilibrium)
C     -Initial composition (a good guess, important for the convergence
C      of the Newton method, can be obtained by the composition at the
C      previous iteration)
        DO IS = 1, NS
          XINI(IS) = X(IS); YINI(IS) = Y(IS)
        ENDDO
C     -The gas composition can also be directly provided by the CFD code
C     for chemical nonequilibrium simulations
C-----------------------------------------------------------------------
C     A. Mole fractions as a function of p and T
        CALL COMPOSITION (WR1, LWR1, WI, LWI, TH, P, XN, XINI, X)
C      => output variable: X(1:NS) mole fractions 
        CALL COMPOSITION (WR1, LWR1, WI, LWI, THP, P, XN, X, XP)
C     Number density of the mixture
        CALL NUMBERD (WR1, LWR1, P, TH, TE, X, ND)
C      => output variable: ND
C     Mass density for the mixture
        CALL DENSITY (WR1, LWR1, X, ND, RHO)
C      => output variable: RHO
C     Molar mass of the mixture
        CALL MOLARMASS (WR1, LWR1, RHO, ND, MM)
C      => output variable: MM 
C     Particle mass of the species 
        CALL MOLE2MASS (WR1, LWR1, ND, X, RHOI)
C      => output variable: RHOI(1:NS)
        CALL MOLE2MASSFRAC (WR1, LWR1, X, Y5)
C      => output variable: Y5(1:NS)
C-----------------------------------------------------------------------
C     B. Mass fractions as a function of rho and T
        CALL MASSCOMPOSITION (WR1, LWR1, WI, LWI, TH, RHO, XN, YINI, Y)
C      => output variable: Y(1:NS) mass fractions 
C     Pressure
        CALL MASS2PRESSURE (WR1, LWR1, RHO, TH, TE, Y, P2)
C      => output variable: P2 pressure 
C     Mole fractions 
        CALL MASS2MOLE (WR1, LWR1, Y, X2)
C      => output variable: X2(1:NS) 
C     Element mole fractions 
        CALL MOLE2ELEMENT (WI, LWI, X2, XN2)
C      => output variable: XN2(1:NC)

        RHO2 = 0.D0
        DO IS = 1, NS
          RHOI2(IS) = Y(IS) *RHO
          RHO2 = RHO2 + RHOI2(IS)
        ENDDO
C-----------------------------------------------------------------------
C     Species mass
        CALL SETMASS (WR1, LWR1, MI)
C      => output variable: MI(1:NS) species molar mass
C-----------------------------------------------------------------------
C     Species enthalpy per unit mass
        CALL ENTHALPYMASS (WR1, LWR1, WI, LWI, TH, TE, TR, TV, P, HM1, 
     &                     HM2, HM3, HM4, HM5, HM6)
C      => output variables:
C     -HM1(1:NS): total species enthalpy 
C                 (including the formation enthalpy)
C     -HM2(1:NS): translational species enthalpy
C     -HM3(1:NS): electronic species enthalpy
C     -HM4(1:NS): rotational species enthalpy
C     -HM5(1:NS): vibrational species enthalpy
C     -HM6(1:NS): formation species enthalpy
C-----------------------------------------------------------------------
C     Species energy per unit mass
        CALL ENERGYMASS (WR1, LWR1, WI, LWI, TH, TE, TR, TV, P, EM1, 
     &                   EM2, EM3, EM4, EM5, EM6)
C      => output variables:
C     -EM1(1:NS): total species energy 
C                 (including the formation enthalpy)
C     -EM2(1:NS): translational species energy
C     -EM3(1:NS): electronic species energy
C     -EM4(1:NS): rotational species energy
C     -EM5(1:NS): vibrational species energy
C     -EM6(1:NS): formation species enthalpy
C-----------------------------------------------------------------------
C     Specific heat per unit mole
        CALL ENERGYMASS (WR1, LWR1, WI, LWI, THP, TEP, TRP, TVP, P, 
     &                   EM1P, EM2P, EM3P, EM4P, EM5P, EM6P)
        DO IS = 1, NS
          CPE       = (EM3P(IS) -EM3(IS)) / (EPS *TH)
          CPR       = (EM4P(IS) -EM4(IS)) / (EPS *TH)
          CPV       = (EM5P(IS) -EM5(IS)) / (EPS *TH)
          CPINT(IS) = MI(IS) *(CPE +CPR +CPV)
        ENDDO

      RHOE = 0.D0
      DO IS = 1, NS
        RHOE     = RHOE +RHOI(IS) *EM1(IS)
      ENDDO
C-----------------------------------------------------------------------
C     Temperature in thermal equilibrium models as a function of rho*e 
C     and rho_i
C     -Temperature guess based on a linear approximation (vibration and 
C      electronic energies neglected). Temperature at the previous time
C      iteration should be used as guess in a CFD code
      CALL TGUESS (WR1, LWR1, WI, LWI, RHOE, RHOI, TINI)
C      => output variable: TINI 
C     -Newton method
      CALL TCNEQNEWTON (WR1, LWR1, WI, LWI, RHOE, RHOI, TINI, T2)
C      => output variable: T2
C-----------------------------------------------------------------------
C     Temperature and composition in thermal equilibrium models as a
C     function of rho*e and rho
C     -Temperature and composition guess such that YINI(TINI,RHO)
C      Values at the previous time iteration should be used as guess in 
C      a CFD code
       TINI = TH *1.2D0
       DO IS = 1, NS
         YINI(IS) = 1.D0
       ENDDO
      CALL MASSCOMPOSITION (WR1, LWR1, WI, LWI, TINI, RHO, XN, YINI, Y3)
       DO IS = 1, NS
         YINI(IS) = Y3(IS) 
       ENDDO

C     -Newton method
      CALL TCEQNEWTON (WR1, LWR1, WI, LWI, RHOE, RHO, XN, YINI, TINI, 
     &                 T3, Y4)
C      => output variables: T3
C                           Y4(1:NS)
C-----------------------------------------------------------------------
C     Chemical production rates expressed in [mole m^-3 s^-1].
      CALL ARRHENIUS (WR1, LWR1, WR3, LWR3, WI, LWI, Y5, TOLY, P, TVEC, 
     &                RHO, OMEGA)
C      => output variable: OMEGA(1:NS)
C-----------------------------------------------------------------------
C     Initialize data dependent on temperature and pressure
C     The routine COLLISION is called after every pressure, temperature,
C     mole fractions updates
C     Beware: use X and not XTOL (to obtain accurate electron transport 
C     properties)
      CALL COLLISION (WR1, LWR1, WR2, LWR2, TH, TE, ND, X)
C      => output variable: WR2
C-----------------------------------------------------------------------
C     BEWARE! A small number is added to the mole fractions to compute
C             the transport properties, respecting the mass
C             conservation. Don't used the modified mole fractions as
C             initial guess to compute the new composition, otherwise
C             the Newton convergence will be killed.
      CALL COMPOTOL (X, TOLX, XTOL)
C      => output variable: XTOL(1:NS) 
C-----------------------------------------------------------------------
C     Viscosity (conjugate gradient)
      CALL ETACG (WR1, LWR1, WR2, LWR2, XTOL, ETA)
C      => output variable: ETA 
C-----------------------------------------------------------------------
C     Thermal conductivity and diffusion velocities
C     -Eucken corrections
      CALL LAMBDAINT (WR1, LWR1, WR2, LWR2, CPINT, XTOL, LAMBDAINTH)
C      => output variable: LAMBDAINTH 
C     -Heavy-particle translational thermal conductivity and thermal
C      diffusion ratios
      CALL LAMBDACHICG (WR1, LWR1, WR2, LWR2, XTOL, LAMBDATH, CHIH)
C      => output variables: LAMBDATH 
C                           CHIH(1:NS) with CHIH(1) = 0.D0 IF NE = 1 
C     -Electron translational thermal conductivity and thermal
C      diffusion ratios
      IF (NE /= 0) THEN
        CALL TDIFE (WR1, LWR1, WR2, LWR2, XTOL, TE, CHIE, 3)
C      => output variable: CHIE(1:NS) 
        CALL LAMBDAE (WR1, LWR1, WR2, LWR2, XTOL, TE, LAMBDATE, 3)
C      => output variable: LAMBDATE 
      ELSE
        LAMBDATE = 0.D0
        DO IS = 1, NS
          CHIE(IS) = 0.D0
        ENDDO
      ENDIF
C-----------------------------------------------------------------------
C     -The driving forces are concentration gradients due to unity 
C      temperature gradient (constant pressure).
      DO IS = 1, NS
        DF(IS) = (XP(IS) -X(IS)) /(TH *EPS) +(CHIE(IS) +CHIH(IS)) /TH 
      ENDDO 
C     -First-order Stefan-Maxwell equations for heavy particles
      CALL CORRECTION (WR1, LWR1, WR2, LWR2, XTOL, 1, FIJ)
C      => output variable: FIJ(1:NS*(NS-1)/2)
C     -Stefan-Maxwell equations (Ramshaw): diffusion fluxes and 
C      ambipolar electric field 
      IF (NE /= 0) THEN
        CALL SMRAMCG (WR1, LWR1, WR2, LWR2, WI, LWI, XTOL, TH, TE, ND,
     &            DF, FIJ, JDIF, EAMB)
C      => output variables: JDIF(1:NS)
C                           EAMB
      ELSE
        CALL SMNEUTCG (WR1, LWR1, WR2, LWR2, XTOL, ND, DF, FIJ, JDIF) 
C      => output variables: JDIF(1:NS)
        EAMB = 0.D0
      ENDIF
      CALL DIJFICK (WR1, LWR1, WR2, LWR2, WI, LWI, XTOL, ND, DIJ)
      DO IS = 1, NS
        JDIF2(IS) = 0.D0
        DO JS = 1, NS
          JDIF2(IS) = JDIF2(IS) -RHOI(IS) *DIJ(IS,JS) *DF(JS)
        ENDDO
      ENDDO
C-----------------------------------------------------------------------
C     -Reactive thermal conductivity based on the Stefan-Maxwell 
C      equations (valid in both chemical equilibrium)
      LAMBDAR = 0.D0; LAMBDAR2 = 0.D0
      DO IS = 1, NS
        LAMBDAR   =  LAMBDAR  -(EM1(IS) +2.D0/3.D0 *EM2(IS))* JDIF(IS)
        LAMBDAR2  =  LAMBDAR2 -(EM1(IS) +2.D0/3.D0 *EM2(IS))* JDIF2(IS)
      ENDDO
C-----------------------------------------------------------------------
C     -Total thermal conductivity 
C      *in LTE, the following routine can be used
      CALL KCONDUCTIVITY (WR1, LWR1, WR2, LWR2, WI, LWI, P, TH, X, XTOL,
     &                    XN, EPS, LAMBDATOT)
C      => output variables: LAMBDATOT 
C      *Otherwise, add the various contributions
      LAMBDATOT2 = LAMBDAINTH +LAMBDATH +LAMBDATE +LAMBDAR
      LAMBDATOT2 = LAMBDAINTH +LAMBDATH +LAMBDATE 
C      *in chemical nonequilibrium & thermal equilibrium, the following 
C       routine can be used
      CALL KCONDUCTIVITYNEQ (WR1, LWR1, WR2, LWR2, WI, LWI, P, TH, X,
     &                       XTOL, XN, EPS, LAMBDATOT3)
C      => output variables: LAMBDATOT3 
      CALL XTOXN (WI, LWI, X, XN1) 
      WRITE(*,*) (XN1(IS), IS=1,NC)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     Output data
C        WRITE(99,2010)  TH, (OMEGA(IS), IS = 1,NS) 
C        WRITE(99,2010)  TH, (Y(IS)-Y4(IS), IS = 1,NS) 
        WRITE(99,2010)  TH, LAMBDAR, LAMBDAR2 
2010  FORMAT (100E14.6)

      ENDDO
      ENDDO

      CLOSE(99)

      END PROGRAM INTERFACE77 
C-----------------------------------------------------------------------

C--------------------------- -------------------------------------------
      SUBROUTINE XTOXN (WI, LWI, X, XN)
C-----------------------------------------------------------------------
C     This subroutine computes the elemental mole fraction based on the 
C     species mole fractions
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER LWI
      INTEGER WI(1:LWI)
      DOUBLE PRECISION X(1:NS), XN(1:NC)
C-----------------------------------------------------------------------
      INTEGER I,J
      DOUBLE PRECISION SUM, NTOT
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      DO J = 1, NS
      WRITE(*,*) (WI(INUIK+(J-1)*NS+I-1), I=1,NS)
      ENDDO
      NTOT = 0.D0
      DO I = 1, NC
        XN(I) = 0.D0
        DO J = 1, NS
          XN(I) = XN(I) +WI(INUIK+(I-1)*NS+J-1) *X(J)
        ENDDO
        NTOT = NTOT +XN(I)
      ENDDO
      DO I = 1, NC
        XN(I) = XN(I) / NTOT
      ENDDO

C-----------------------------------------------------------------------
      END SUBROUTINE XTOXN 
C--------------------------- -------------------------------------------
