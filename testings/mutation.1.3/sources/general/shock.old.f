
C     SHOCKING code


C     Ecole Centrale Paris, Chatenay-Malabry (France), 05-31-2005

C     Stanford University, Palo Alto (USA), 12-31-2006

C     Thierry E. Magin


C     von Karman Institute, Rhode-Saint-Genese (Belgium), 12-31-2006

C     Marco Panesi

C-----------------------------------------------------------------------
      PROGRAM SHOCKING
C-----------------------------------------------------------------------
C     This program computes the flowfield in a shock tube. 
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE '../general/memory.cmn'
      INCLUDE '../general/memoryshock.cmn'
      INCLUDE '../general/memoryradiation.cmn'
C-----------------------------------------------------------------------
      INTEGER LRW, LIW, N, NB, IWRITE, C
      PARAMETER (LRW= 900000, LIW = 1200)
      INTEGER  NEQ, IS, IMOD, ITERMAX, ITER, INDX(1:3), I, IV, ICUT, MF,
     &         ITOL, ITASK, ISTATE, IOPT, IWORK(LRW), LCHAR, IC1, IC2,
     &         LSPECIES
      DOUBLE PRECISION X(1:NMAX), P1, T1, US,
     &                 P2, T2, U2, RESMIN, RES, RESINI, R0(1:3), 
     &                 J(1:3,1:3), D, MM, RHO1, POSITION, Z(NMAX+3), 
     &                 ZP(NMAX+3), ND, RHO, TINI, TOUT, TMAX,
     &                 DT, RTOL, ATOL, RWORK(LRW), ALPHA, XTOL
      CHARACTER(10) MIXTURE, REACTION, TRANSFER      
      CHARACTER(100) PATH
      EXTERNAL FUN, DUMMY 
C-----------------------------------------------------------------------
      PATH = '..'

      CALL READSHOCK (PATH, MIXTURE, REACTION, TRANSFER, P1, T1, US)

C     Link with the MUTATION library
C     ------------------------------
      IMOD = 0

C     Pseudo-dynamic allocation
      CALL LENGTH (PATH, MIXTURE, REACTION, TRANSFER, LWR1, LWR2, LWR3, 
     &             LWR4, LWI, LWC, NS, NE, NC, NREA, NV, NMAX, NVIB,
     &             NELEQ)
      IF (     ( LWR1 > LWR1MAX ) .OR. ( LWR2 > LWR2MAX ) 
     &    .OR. ( LWR3 > LWR3MAX ) .OR. ( LWR4 > LWR4MAX )
     &    .OR. ( LWI  > LWIMAX  ) .OR. ( LWC  > LWCMAX  )
     &    .OR. ( NS   > NMAX    ) .OR. ( NV   > NMAX    ) ) THEN
         WRITE(*,*) LWR1, LWR2, LWR3, LWR4, LWI, LWC, NS, NV
         WRITE(*,*) 'Error: work arrays out of space.'
         PAUSE
      ENDIF

C     Number of temperatures
      NT = 1 +NVIB +NELEQ

      WRITE(*,*)
      WRITE(*,*) 'Shock-tube'
      WRITE(*,*) '=========='
      CALL PRESENTATION (3)
      IF (NT==1) THEN
        WRITE(*,*) '    in thermal equilibrium.' 
      ELSEIF (NT>1) THEN
        WRITE(*,*) '    in thermal nonequilibrium.' 
      ENDIF 
      WRITE(*,*)

C     Initialize data independent on temperature and pressure
      CALL INITIALIZE (PATH, MIXTURE, REACTION, TRANSFER, WR1, LWR1, 
     &                 WR3, LWR3, WR4, LWR4, WI, LWI, WC, LWC, IMOD)
      IF (NVIB >= 1) THEN
        IC1 = 0
        WRITE(*,*) 
        DO IV = 1, NVIB     
          WRITE(*,*) '    -Energy equation', IV,' for the vibration of:'
          DO I = 1, WI(INVIBSPEI+IV-1)
            IC1 = IC1 +1
            IC2 = 0
            DO IS = 1, WI(IVIBSPEI+IC1-1)-1
              IC2 = IC2 +WI(ILNAMEI+IS-1)
            ENDDO
            LSPECIES = WI(ILNAMEI+WI(IVIBSPEI+IC1-1)-1)
            WRITE(*,*) '      ', (WC(INAMEI+IC2+IS-1), IS=1,LSPECIES)
          ENDDO 
        ENDDO
        WRITE(*,*)
      ENDIF

C     Variable passed in common in both memory.cmn and memoryshock.cmn
C     ----------------------------------------------------------------
      NNS = NS; NNE = NE; NNVIB= NVIB; NNELEQ = NELEQ
      IC1 = 0
      DO IV = 1, NVIB
        NNVIBSPE(IV) = WI(INVIBSPEI+IV-1)
        DO I = 1, WI(INVIBSPEI+IV-1)
          IC1 = IC1 +1
          IIVIBSPE(IC1) = WI(IVIBSPEI+IC1-1)
        ENDDO
      ENDDO

C     Initial composition
C     -------------------
      XTOL = 1.D-32
      CALL READCOMPO (PATH, WI, LWI, WC, LWC, XTOL, X)

C     Rankine-Hugoniot's jump relations 
C     ---------------------------------
C     Initialization (cold gas)
      MM = 0.D0
      DO IS= 1, NS
        MM  =  MM + X(IS) *WR1(IMI+IS-1)
      ENDDO

      CALL RHCOLD (WR1(IUR)/MM, US, P1, T1, U2, P2, T2)
      WRITE(*,*) '*Post-shock conditions (cold)'
      WRITE(*,*) '-----------------------------'
      WRITE(*,*) ' -P2:                  ', P2/101325.,'[atm]'
      WRITE(*,*) ' -P2:                  ', P2,'[Pa]'
      WRITE(*,*) ' -T2:                  ', T2,'[K]'
      WRITE(*,*) ' -US-U2:               ', US-U2, '[m/s]' 
      WRITE(*,*)
    
      P2 = 54000.
      T2 = 3900.
      U2 = 3000.-434. 
 
      CALL SYSTEMRH (WR1, LWR1, WI, LWI,  US, P1, T1, U2, P2, T2, X, 
     &               R0, J, NT)  

      RES = R0(1) *R0(1) +R0(2) *R0(2) +R0(3) *R0(3)
      RESINI = RES
      RES = 1.D0

C     Test (Newton's method, residual and maximum number of iterations)
      RESMIN = 1.D-10; ITERMAX = 100000; ITER = 0
C     Underrelaxation parameter ALPHA \in ]0, 1]
      ALPHA = 1.D0
      RESMIN = RESMIN *RESMIN

      DO WHILE ( (ITER < ITERMAX) .AND.  (RES > RESMIN) )
        CALL LUDCMP (J, 3, 3, INDX, D)
        CALL LUBKSB (J, 3, 3, INDX, R0)

        U2 = U2 +ALPHA *R0(1)
        P2 = P2 +ALPHA *R0(2)
        T2 = T2 +ALPHA *R0(3)

        CALL SYSTEMRH (WR1, LWR1, WI, LWI,  US, P1, T1, U2, P2, T2, X, 
     &                 R0, J, NT)  
        RES  = (R0(1) *R0(1) +R0(2) *R0(2) +R0(3) *R0(3)) /RESINI
        ITER = ITER +1                 
      ENDDO

      WRITE(*,*) '*Post-shock conditions (high temperatures)'
      WRITE(*,*) '-------------------------------------------'
      WRITE(*,*) ' -Number of Newton iterations:', ITER
      WRITE(*,*) ' -Residual:            ', DSQRT(RES)
      WRITE(*,*) ' -P2:                  ', P2/101325.,'[atm]'
      WRITE(*,*) ' -P2:                  ', P2,'[Pa]'
      WRITE(*,*) ' -T2:                  ', T2,'[K]'
      WRITE(*,*) ' -US-U2:               ', US-U2, '[m/s]'

C     One-dimensional Euler equations after the shock
C     -----------------------------------------------
C     Parameters for LSODE
C     ITOL = 1 or 2 according as atol (below) is a scalar or array.
      ITOL  = 1
C     RTOL = relative tolerance parameter (scalar).
      RTOL  = 1.D-6
C     ATOL = absolute tolerance parameter (scalar or array).
C          the estimated local error in y(i) will be controlled so as
C          to be roughly less (in magnitude) than
C          ewt(i) = rtol*abs(y(i)) + atol     if itol = 1, or
C          ewt(i) = rtol*abs(y(i)) + atol(i)  if itol = 2.
C          thus the local error test passes if, in each component,
C          either the absolute error is less than atol (or atol(i)),
C          or the relative error is less than rtol.
C          use rtol = 0.0 for pure absolute error control, and
C          use atol = 0.0 (or atol(i) = 0.0) for pure relative error
C          control.  caution... actual (global) errors may exceed these
C          local tolerances, so choose them conservatively.
      ATOL  = 1.D-15
C     ITASK = 1 for normal computation of output values of y at t = 
C          tout.
      ITASK = 1
C     ISTATE = integer flag (input and output). Set ISTATE = 1
      ISTATE= 1
C     IOPT = 0 to indicate no optional inputs used.
      IOPT  = 0
C     MF = method flag.  standard values are:
C          10 for nonstiff (adams) method, no jacobian used.
C          21 for stiff (bdf) method, user-supplied full jacobian.
C          22 for stiff method, internally generated full jacobian.
C          24 for stiff method, user-supplied banded jacobian.
C          25 for stiff method, internally generated banded jacobian.
      MF    = 22
C     NEQ = number of first order ode-s.
      NEQ = NS +1 +NT

C     Tolerance for mole fractions in Arrhenius rates
      YTOL = 1.D-20

C     Initial point
C     No electrons such that the perfect gas equation in thermal 
C     equilibrium is still valid.
      RHO1 = P1*MM /(WR1(IUR) *T1)
      MDOT = RHO1 *US
      DO I = 1, NS
        Z(I) = X(I) *WR1(IMI+I-1) /MM
      ENDDO
      Z(NS+1) = US -U2
      Z(NS+2) = T2
      DO I = 1, NVIB+NELEQ
        Z(NS+2+I) = T1
      ENDDO
C     Loop      
      OPEN(UNIT=INOUT1,FILE='../output/shockflow.dat',STATUS='UNKNOWN')
      TMAX  = 2.D-1
Cthierry
C      TMAX  = 20.D-1
Cthierry
      TINI  = 0.D0 
      TOUT  = TINI
      WRITE(*,*) 
      WRITE(*,*)'*Flow field'
      WRITE(*,*)'-----------'
      WRITE(*,*)' -First position: x = 0. [m]'
      K = 0
      CALL WRITESOL(TOUT, NEQ, Z)
      N = 5; DT   = 10.D0**(-N)
C     Write solution every TOUT = 10^(-3)
      NB = 10**(N-4)
Cthierry
C      NB = 10**(N-3)
Cthierry
      I = 1; TOUT = TOUT +DT
      DOWHILE (TOUT < TMAX)
        CALL LSODE (FUN, NEQ, Z, TINI, TOUT, ITOL, RTOL, ATOL, ITASK,
     &              ISTATE, IOPT, RWORK, LRW, IWORK, LIW, DUMMY, MF)
        IF (I-INT(I/NB)*NB == 0) THEN          
          CALL WRITESOL(TOUT, NEQ, Z)
        ENDIF
        I = I +1; TOUT = TOUT +DT     
        IF (ISTATE <= -2) THEN
          WRITE(*,*) 'ISTATE = ', ISTATE
        STOP
        ENDIF
      ENDDO 
      CLOSE(INOUT1)
      WRITE(*,*)' -Final position: x =', TMAX, ' [m]'
      WRITE(*,*)' -Corresponding time: t =', TIMEK(K), ' [s]'
      WRITE(*,*)' -Number of cells in the mesh: K =', K
      WRITE(*,*)
      IF (K>KMAX) THEN
         WRITE(*,*) 'K', K, ':', KMAX
         WRITE(*,*) 'Error: work arrays out of space.'
      ENDIF

      CALL PRESENTATION (2)

      END PROGRAM SHOCKING
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE DUMMY()
      END SUBROUTINE DUMMY
C----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE RHCOLD (R, US, P1, T1, U2, P2, T2)
C-----------------------------------------------------------------------
C     This subroutine computes the variables after the shock based on 
C     Rankine-Hugoniot relations in the case of a cold gas (specific 
C     heat ration = 1.4).
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION R, US, P1, T1, U2, P2, T2
C-----------------------------------------------------------------------
      DOUBLE PRECISION MS, G, GP1, C1, MS2, A, B, C
C-----------------------------------------------------------------------
      G   = 1.4D0
      GP1 = G +1.D0

      C1  = SQRT(G * R *T1)
      MS  = US / C1
      MS2 = MS *MS
      T2 = T1 *(2.D0 *G *MS2 -G +1.D0)*(G -1.D0 +2.D0 /MS2) /(GP1 *GP1)

      U2 = C1 *2.D0 /GP1 *(MS - 1.D0 /MS)
      P2 = P1 *(2.D0 * G *MS2 -G +1.D0) /GP1

      END SUBROUTINE RHCOLD
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SYSTEMRH (WR1, LWR1, WI, LWI,  US, P1, T1, U2, P2, T2,
     &                     X, R0, J, NT)  
C-----------------------------------------------------------------------
C     This subroutine computes the residual and the Jacobian of the 
C     Newton  method for the Rankine-Hugoniot relations in thermal 
C     equilibrium and nonequilibrium situations. Electronic levels of 
C     atoms and molecules are not taken into account. Moreover, it is 
C     assumed that no free electrons are present in the mixture 
C     (OK behind the shock in a shock tube) such that the perfect gas 
C     law in thermal equilibrium is still valid. 
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------

      INTEGER LWR1, LWI, WI(1:LWI), NT
      DOUBLE PRECISION WR1(LWR1), US, P1, T1, U2, P2, T2, X(1:NS), 
     &                 R0(1:3), J(1:3, 1:3)
C-----------------------------------------------------------------------
      INTEGER I
      DOUBLE PRECISION PHI1, RHO1, RHO2, E1, E2, MM, CP2, EPS,
     &                 EPSP1V, EPSP1, PHI12, E11(1:NS), E12(1:NS), 
     &                 E13(1:NS), E14(1:NS), E15(1:NS), E16(1:NS), 
     &                 E21(1:NS), E22(1:NS), E23(1:NS), E24(1:NS), 
     &                 E25(1:NS), E26(1:NS), E21P(1:NS), E22P(1:NS), 
     &                 E23P(1:NS), E24P(1:NS), E25P(1:NS), E26P(1:NS),
     &                 T1V(1:NV), T2V(1:NV), T2E, T2VEPSP1(1:NV)
C-----------------------------------------------------------------------
      EPS = 1.D-2; EPSP1 = 1.D0 +EPS
C     Thermal equilibrium: the rotational and  vibrational energy 
C     populations are Boltzmann distributions at the translational
C     temperature T2.
      IF (NT==1) THEN
        DO I = 1, NV
          T1V(I)  = T1
          T2V(I)  = T2
          T2E     = T2
        ENDDO
        EPSP1V = EPSP1
C     Thermal nonequilibrium: the rotational energy population is a
C     Boltzmann distribution at the translational temperature T2, 
C     the vibrational population is frozen at the temperature T1. 
      ELSE
        DO I = 1, NV
          T1V(I)  = T1
          T2V(I)  = T1
          T2E     = T1
        ENDDO
        EPSP1V = 1.D0
      ENDIF
      DO I = 1, NV
        T2VEPSP1(I) = T2V(I)*EPSP1V
      ENDDO

      CALL ENERGY (WR1, LWR1, WI, LWI, T1, T1, T1, T1V, P1, E11, E12, 
     &             E13, E14, E15, E16)
      CALL ENERGY (WR1, LWR1, WI, LWI, T2, T2E, T2, T2V, P2, E21, E22, 
     &             E23, E24, E25, E26)
      CALL ENERGY (WR1, LWR1, WI, LWI, EPSP1*T2, EPSP1V*T2E, EPSP1*T2,
     &             T2VEPSP1, P2, E21P, E22P, E23P, E24P, E25P, E26P)
      MM = 0.D0; E1 = 0.D0; E2 = 0.D0 ; CP2 = 0.D0
      DO I = 1, NS
        MM  =  MM + X(I) *WR1(IMI+I-1)
Cthierry
C        E1  =  E1 + X(I) *(E11(I) -E13(I))
C        E2  =  E2 + X(I) *(E21(I) -E23(I))
C        CP2 = CP2 + X(I) *(E21P(I) -E21(I) -E23P(I)+E23(I)) /(T2 *EPS)
        E1  =  E1 + X(I) *E11(I)
        E2  =  E2 + X(I) *E21(I)
        CP2 = CP2 + X(I) *(E21P(I) -E21(I)) /(T2 *EPS)
Cthierry
      ENDDO

C     Energy per unit mole => energy per unit mass
      E1  =  E1 /MM
      E2  =  E2 /MM
      CP2 = CP2 /MM

      RHO1  = MM *P1 /(WR1(IUR) *T1)
      RHO2  = MM *P2 /(WR1(IUR) *T2)
      PHI1  = -RHO1 *US
      PHI12 = PHI1 *PHI1

C     Residual
C     --------
      R0(1) = -(U2 +(P2 -P1) /PHI1)
      R0(2) = -(1.D0 /RHO2 -1.D0 /RHO1 +(P2 -P1) /PHI12)
      R0(3) = -(E2 -E1 -(P2*P2 -P1*P1) /(2.D0 *PHI12))

C     Jacobian
C     --------
      J(1,1) = 1.D0
      J(1,2) = 1.D0 /PHI1
      J(1,3) = 0.D0
      J(2,1) = 0.D0  
      J(2,2) = -1.D0 /(RHO2*T2)  +1.D0 /PHI12
      J(2,3) = 1.D0 /(RHO2*T2)
      J(3,1) = 0.D0
      J(3,2) = -P2 /PHI12
      J(3,3) = CP2

      END SUBROUTINE SYSTEMRH
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE FUN (NEQ, X, Z, ZP)  
C-----------------------------------------------------------------------
C     This subroutine computes the source terms of the ODEs (1D Euler 
C     equations). 
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE '../general/memoryshock.cmn'
C-----------------------------------------------------------------------
      INTEGER NEQ
      DOUBLE PRECISION X, Z(1:NEQ), ZP(1:NEQ)
C-----------------------------------------------------------------------
      INTEGER I, ICC, IV
      DOUBLE PRECISION Y(1:NNS), U, TH, TE, CP(1:NNS,1:2), P, FACE,
     &                 H(1:NNS), OMEGA(1:NNS), DET, OMEGASUM, HSUM, 
     &                 MASS(1:NNS), RHO, R, YSUM, CPRTSUM, CPVSUM, A, B,
     &                 C, D, E, F, OMEGAVT(1:NNVIB+1), TVEC(1:NT+1), 
     &                 TVIB(1:NNVIB+1), FAC(1:NNVIB), CPVIB(1:NNVIB), 
     &                 CPVIBSUM, OMEGAVV(1:NNVIB+1), OMEGAVE(1:NNVIB+1),
     &                 OMEGAET, OMEGAEV
C-----------------------------------------------------------------------
      DO I = 1, NNS
        Y(I) = Z(I)
      ENDDO
      U  = Z(NNS+1)
      TH = Z(NNS+2)
      IF (NNVIB == 0) THEN
        TVIB(1) = TH
      ELSE
        DO I = 1, NNVIB
          TVIB(I) = Z(NNS+2+I)
        ENDDO
      ENDIF
      IF (NNELEQ == 0) THEN
        TE = TVIB(1)
      ELSEIF (NNELEQ == 1) THEN
        TE = Z(NNS+NNVIB+3)
      ENDIF
      TVEC(1)      = TH
      DO I = 1, NNVIB
        TVEC(1+I)  = TVIB(I)
      ENDDO
      TVEC(NNVIB+2) = TE

      RHO = MDOT /U
      CALL SOURCE (LWR1, WR1, LWR3, WR3, LWR4, WR4, WI, LWI, Y, YTOL, U,
     &             TVEC, RHO, R, CP, P, H, OMEGA, OMEGAVT, OMEGAVE,
     &             OMEGAVV, OMEGAET, MASS)  
      HSUM = 0.D0
      DO I = 1, NNS
        ZP(I)    = OMEGA(I) *MASS(I) /MDOT
        HSUM     = HSUM + H(I) *OMEGA(I) *MASS(I)
      ENDDO
      OMEGASUM = 0.D0; YSUM = 0.D0 ; CPRTSUM = 0.D0; CPVSUM = 0.D0
      DO I = 2, NNS
        OMEGASUM = OMEGASUM +OMEGA(I)
        YSUM     = YSUM    +Y(I) /MASS(I)
        CPRTSUM  = CPRTSUM +Y(I) *CP(I,1)
        CPVSUM   = CPVSUM  +Y(I) *CP(I,2)
      ENDDO
Cthierry
      IF (NNE==0) THEN
        I = 1
          CPRTSUM  = CPRTSUM +Y(I) *CP(I,1)
          CPVSUM   = CPVSUM  +Y(I) *CP(I,2)
      ENDIF
Cthierry

      OMEGASUM = OMEGASUM +OMEGA(1) *TE /TH
      A = MDOT *U /(R* TH) -RHO *(YSUM +Y(1) /MASS(1) *TE /TH)
      B = MDOT /TH *YSUM
      C = -OMEGASUM
      D = MDOT *U
      E = MDOT *CPRTSUM
      F = -HSUM
      IF (NNVIB==0) THEN
        B = B +MDOT *Y(1) /(TH *MASS(1))
        E = E +MDOT *(Y(1) *CP(1,1) +CPVSUM)
      ELSE
        ICC = 0
        DO IV = 1, NNVIB
          CPVIB(IV) = 0.D0
          DO I = 1, NNVIBSPE(IV)
            ICC       = ICC +1
            CPVIB(IV) = CPVIB(IV) +Y(IIVIBSPE(ICC)) *CP(IIVIBSPE(ICC),2)
          ENDDO
          FAC(IV) = (OMEGAVT(IV) +OMEGAVE(IV) +OMEGAVV(IV)) /CPVIB(IV)
          CPVSUM  = CPVSUM -CPVIB(IV)
          F = F -CPVIB(IV) *FAC(IV)
        ENDDO
        IF (NNELEQ==0) THEN
          FACE = FAC(1)
        ELSEIF (NNELEQ==1) THEN
          OMEGAEV = 0.D0
          DO IV = 1, NNVIB
            OMEGAEV = OMEGAEV -OMEGAVE(IV)
          ENDDO
          FACE = (OMEGAEV +OMEGAET) /(Y(1) *0.6D0 *CP(1,1))
        ENDIF
        C = C -Y(1) /(TH * MASS(1)) *FACE
        F = F -Y(1) *CP(1,1) *FACE -CPVSUM *FAC(1)
      ENDIF
      DET = A *E - B *D
      ZP(NNS+1) = (C *E -B *F) /DET
      ZP(NNS+2) = (A *F -C *D) /DET

      DO I = 1, NNVIB
        ZP(NNS+2+I) = FAC(I) /MDOT
      ENDDO
      IF (NNELEQ==1) THEN
        ZP(NNS+NNVIB+3) = FACE /MDOT
      ENDIF
    
      END SUBROUTINE FUN
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SOURCE (LWR1, WR1, LWR3, WR3, LWR4, WR4, WI, LWI, Y,
     &                   YTOL, U, TVEC, RHO, R, CP, P, H, OMEGA, 
     &                   OMEGAVT, OMEGAVE, OMEGAVV, OMEGAET, MASS) 
C-----------------------------------------------------------------------
C     This subroutine computes quantities necessary to evaluate the 
C     source terms of the 1D Euler equations. The electronic levels are 
C     not accounted for.
C     -Specific heat (per unit mass) CP(I,K)
C       *I:     species.
C       *K = 1: translation (at TH for heavy particles or TE for
C               electrons) + rotation (at TH).
C       *K = 2: vibration (at TV(IV)).
C     -Enthalpy H(I): translation (at TH or TE), rotation (at TH), 
C      vibration (at TV(IV), and formation.
C-----------------------------------------------------------------------
      IMPLICIT NONE 
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      INTEGER LWR1, LWR3, LWR4, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), WR3(1:LWR3), WR4(1:LWR4), Y(1:NS), 
     &                 YTOL, U,RHO, CP(1:NS,1:2), P, H(1:NS), 
     &                 OMEGA(1:NS), MASS(1:NS), R, OMEGAVT(1:NVIB+1), 
     &                 OMEGAVE(1:NVIB+1), OMEGAVV(1:NVIB+1), OMEGAET,
     &                 TVEC(1:NVIB+2)
C-----------------------------------------------------------------------
      INTEGER I, IR, IS, IC2, IC, IV, NVIBMODE, IVT, J, K
      DOUBLE PRECISION EPS, EPSP1, H1(1:NS), H2(1:NS), H3(1:NS), 
     &                 H4(1:NS), H5(1:NS), H6(1:NS),
     &                 H1T(1:NS), H2T(1:NS), H3T(1:NS), H4T(1:NS),
     &                 H5T(1:NS), H6T(1:NS), H1P(1:NS), H2P(1:NS),
     &                 H3P(1:NS), H4P(1:NS), H5P(1:NS), H6P(1:NS),
     &                 TT(1:NV), EPSP1TV(1:NV), VT , TH, TE, TV(1:NV) 
C-----------------------------------------------------------------------
      EPS = 1.D-2; EPSP1 = 1.D0 +EPS
      R = WR1(IUR)

C     Translational temperature of electrons
      TH = TVEC(1)
C     Translational temperature of electrons
      TE = TVEC(NVIB+2)
C     Vibrational temperatures
      IC = 0
      DO IS = 1, NS
        NVIBMODE = WI(IVIBI+IS-1)
        DO IV = 1, NVIBMODE
          IC = IC +1
          TV(IC)   = TVEC(WI(IVIBTEMPI+IS-1))
        ENDDO
      ENDDO
C     Vectors of size (1:NV)
      DO I = 1, NV
        TT(I)      = TH
        EPSP1TV(I) = EPSP1 *TV(I)
      ENDDO
C     Pressure
      P = 0.D0
      DO I = 2, NS
        MASS(I) =WR1(IMI+I-1)
        P = P +Y(I) /MASS(I) *TH
      ENDDO 
      I = 1
        MASS(I) =WR1(IMI+I-1)
        P = P +Y(I) /MASS(I) *TE
      P = RHO *R *P
      
C     Specific heat (per unit mass) and enthalpy
      CALL ENTHALPYMASS (WR1, LWR1, WI, LWI, TH, TE, TH, TV, P, H1, 
     &                   H2, H3, H4, H5, H6)
      CALL ENTHALPYMASS (WR1, LWR1, WI, LWI, TH, TH, TH, TT, P, H1T,
     &                   H2T, H3T, H4T, H5T, H6T)
      CALL ENTHALPYMASS (WR1, LWR1, WI, LWI, EPSP1*TH, EPSP1*TE, 
     &                   EPSP1*TH, EPSP1TV, P, H1P, H2P, H3P, H4P, H5P, 
     &                   H6P)

      DO I = 2, NS
        CP(I,1) = (H2P(I) -H2(I) +H4P(I) -H4(I)) /(TH *EPS)
Cthierry electronic levels
C        CP(I,2) = (H5P(I) -H5(I)) /(TVEC(WI(IVIBTEMPI+I-1)) *EPS)
        CP(I,2) = (H5P(I) -H5(I)) /(TVEC(WI(IVIBTEMPI+I-1)) *EPS)
     &          + (H3P(I) -H3(I)) /(TE *EPS)
Cthierry
        H(I)    = H2(I) +H3(I) +H4(I) +H5(I) + H6(I)
      ENDDO 
      I = 1
Cthierry electronic levels + mixture without electrons
C        CP(I,1) = (H2P(I) -H2(I)) /(TE *EPS)
C        CP(I,2) = 0.D0
C        H(I)    = H2(I) +H6(I)
        CP(I,1) = (H2P(I) -H2(I) +H4P(I) -H4(I)) /(TE *EPS)
        CP(I,2) = (H5P(I) -H5(I)) /(TVEC(WI(IVIBTEMPI+I-1)) *EPS)
     &          + (H3P(I) -H3(I)) /(TE *EPS)
        H(I)    = H2(I) +H3(I) +H4(I) +H5(I) + H6(I)
Cthierry

C     Chemical production rates 
      CALL ARRHENIUS (WR1, LWR1, WR3, LWR3, WI, LWI, Y, YTOL, P, TVEC, 
     &                RHO, OMEGA)

C     Energy transfer rates
      IF (NVIB == 0) THEN
        OMEGAVT(1) = 1.D0
        OMEGAVE(1) = 1.D0
        OMEGAVV(1) = 1.D0
      ELSEIF (NVIB > 0) THEN
        CALL OMEGAVTRANSFER (LWR1, WR1, LWR4, WR4, LWI, WI, P, TVEC, Y, 
     &                       RHO, H5T, H5, OMEGAVT, OMEGAVE, OMEGAVV)
      ENDIF
      IF (NELEQ==0) THEN
        OMEGAET = 0.D0
      ELSEIF (NELEQ==1) THEN
        OMEGAET = 1.D0
      ENDIF

C     LTE limit
C      DO I = 1, NS
C        OMEGA(I) = 1.D6*OMEGA(I)
C      ENDDO 

      END SUBROUTINE SOURCE
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE WRITESOL (POSITION, NEQ, Z)  
C-----------------------------------------------------------------------
C     This subroutine writes the solution at the current position.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE '../general/memory.cmn'
      INCLUDE '../general/memoryshock.cmn'
      INCLUDE '../general/memoryradiation.cmn'
C-----------------------------------------------------------------------
      INTEGER NEQ
      DOUBLE PRECISION POSITION, Z(1:NEQ)
C-----------------------------------------------------------------------
      INTEGER I
      DOUBLE PRECISION Y(1:NNS), U, TH, CP(1:NNS,1:2), MM, MACH, P, 
     &                 H(1:NNS), OMEGA(1:NNS), MASS(1:NNS), X(1:NNS),
     &                 RHO, R, OMEGAVT(1:NNVIB+1), OMEGAVE(1:NNVIB+1),
     &                 OMEGAVV(1:NNVIB+1), N(1:NNS), ND, TVEC(1:NT+1),
     &                 TVIB(1:NNVIB), OMEGAET
C-----------------------------------------------------------------------
      DO I = 1, NNS
        Y(I) = Z(I)
      ENDDO
      U  = Z(NNS+1)
      TH = Z(NNS+2)
      IF (NNVIB == 0) THEN
        TVIB(1) = TH
      ELSE
        DO I = 1, NNVIB
          TVIB(I) = Z(NNS+2+I)
        ENDDO
      ENDIF
      TVEC(1)      = TH
      DO I = 1, NNVIB
        TVEC(1+I)  = TVIB(I)
      ENDDO
      TVEC(NNVIB+2) = TVIB(1)
      RHO = MDOT /U

      CALL SOURCE (LWR1, WR1, LWR3, WR3, LWR4, WR4, WI, LWI, Y, YTOL, U,
     &             TVEC, RHO, R, CP, P, H, OMEGA, OMEGAVT, OMEGAVE,
     &             OMEGAVV, OMEGAET, MASS)

      MM = 0.D0
      DO I = 1, NS
        MM = MM + Y(I) /MASS(I)
      ENDDO
      MM = 1.D0 /MM

      DO I = 1, NNS
        X(I) = Y(I) *MM /WR1(IMI+I-1)
      ENDDO

      CALL NUMBERD (WR1, LWR1, P, TH, TVIB(1), X, ND)

      DO I = 1, NNS
C        N(I) = X(I) *ND
        N(I) = Y(I) *RHO
        N(I) = X(I) 
      ENDDO

      K = K +1
      POSITIONK(K) = POSITION
      RHOK(K)      = RHO
      TK(K)        = TH
      TVK(K)       = TVIB(1)
      UK(K)        = U
      IF (K==1) THEN
        TIMEK(K)     = 0.D0
      ELSE
        TIMEK(K)     = TIMEK(K-1)+ 0.5D0 *(POSITIONK(K)-POSITIONK(K-1))
     &                 *(1.D0 /UK(K-1) +1.D0 /UK(K))
      ENDIF
      NTOTK(K)     = ND
      DO I = 1, NNS
        YK(I,K)    = Y(I)
        XK(I,K)    = X(I)
      ENDDO

      WRITE(INOUT1,2010) POSITION, TH, (TVIB(I),I=1,NNVIB), P, U, N
2010  FORMAT (100E14.6)

      END SUBROUTINE WRITESOL
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE READSHOCK (PATH, MIXTURE, REACTION, TRANSFER, P1, T1, 
     &                      US)
C-----------------------------------------------------------------------
C     This subroutine reads the parameters for shocking.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION P1, T1, US
      CHARACTER(10) MIXTURE, REACTION, TRANSFER
      CHARACTER(100) PATH
C-----------------------------------------------------------------------
      INTEGER IN, LPATH, LCHAR
      PARAMETER (IN = 30)
      CHARACTER(4) COM1
      CHARACTER(80) FULLCOM1
C-----------------------------------------------------------------------
      LPATH = LCHAR(PATH)
      OPEN(UNIT=IN,FILE=PATH(1:LPATH)//'/input/shock',STATUS='OLD')
      COM1 = '   '
      DO WHILE (COM1(1:4)/= 'STOP')
        READ(IN,*) FULLCOM1
        COM1 = FULLCOM1(1:4)
        IF (COM1(1:4) == 'Name') THEN
          READ(IN,*) MIXTURE
        ELSEIF (COM1(1:4) == 'Reac') THEN
          READ(IN,*) REACTION
        ELSEIF (COM1(1:4) == 'Ener') THEN
          READ(IN,*) TRANSFER
        ELSEIF (COM1(1:4) == 'Pres') THEN
          READ(IN,*) P1
        ELSEIF (COM1(1:4) == 'Temp') THEN
          READ(IN,*) T1
        ELSEIF (COM1(1:4) == 'Shoc') THEN
          READ(IN,*) US
        ENDIF
      ENDDO
      CLOSE(IN)

      END SUBROUTINE READSHOCK
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE READCOMPO (PATH, WI, LWI, WC, LWC, XTOL, X)
C-----------------------------------------------------------------------
C     This subroutine deals with the species accounted for into the
C     vibrational equation and temperatures for shocking.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      INTEGER LWC, LWI, WI(1:LWI)
      DOUBLE PRECISION X(1:NS), XTOL
      CHARACTER(100) PATH
      CHARACTER WC(1:LWC)
C-----------------------------------------------------------------------
      INTEGER IN, LPATH, LCHAR, I, IC, ISPECIES, ILEN, INDSPECIES(1:NS),
     &        IX
      DOUBLE PRECISION SUM
      PARAMETER (IN = 30)
      CHARACTER(4) COM1
      CHARACTER(10) SPECIES, WORD
      CHARACTER(100) FULLCOM, LINE
C-----------------------------------------------------------------------
      DO I = 1, NS 
        X(I) = 0.D0
      ENDDO
      LPATH = LCHAR(PATH)
      OPEN(UNIT=IN,FILE=PATH(1:LPATH)//'/input/shock',STATUS='OLD')
      COM1 = '   '
      DO WHILE (COM1(1:4)/= 'STOP')
        CALL BLANK (FULLCOM)
        READ(IN,*) FULLCOM
        COM1 = FULLCOM(1:4)
        IF (COM1(1:4) == 'Comp') THEN
          CALL BLANK (SPECIES)
          IC = 0; IX = 0
          CALL BLANK (FULLCOM)
          READ(IN,'(A)') FULLCOM
          CALL RMBLANK (FULLCOM, LINE, ILEN)
          DO I = 1, ILEN
            IC = IC +1
            IF ((LINE(I:I)=='/').OR.(I==ILEN)) THEN
              IF (I==ILEN) THEN
                SPECIES(IC:IC) = LINE(I:I)                
              ENDIF
              CALL TESTSPECIES (SPECIES, WI, LWI, WC, LWC, ISPECIES)
              IF (ISPECIES==0) THEN
                WRITE(*,*) 'Species', SPECIES,
     &                     ' does not exist!!', LINE
                PAUSE
              ELSE
                IX = IX +1
                INDSPECIES(IX) = ISPECIES
              ENDIF
              IC = 0
              CALL BLANK (SPECIES)
            ELSE
              SPECIES(IC:IC) = LINE(I:I)
            ENDIF   
          ENDDO   
          IC = 0; IX = 0
          CALL BLANK (WORD)
          CALL BLANK (FULLCOM)
          READ(IN,'(A)') FULLCOM
          CALL RMBLANK (FULLCOM, LINE, ILEN)
          DO I = 1, ILEN
            IC = IC +1
            IF ((LINE(I:I)=='/').OR.(I==ILEN)) THEN
              IF (I==ILEN) THEN
                WORD(IC:IC) = LINE(I:I)                
              ENDIF
              IX = IX +1
              CALL DTRANSLATE (WORD, X(INDSPECIES(IX)))
              IC = 0
              CALL BLANK (WORD)
            ELSE
              WORD(IC:IC) = LINE(I:I)
            ENDIF      
          ENDDO
        ENDIF
      ENDDO
      SUM = 0
      DO I = 1, NS
        SUM = SUM +X(I)
      ENDDO
      IF (ABS(SUM-1.D0) > 1.D-6) THEN
        WRITE(*,*) SUM, 'Invalid composition'
        PAUSE
      ENDIF
      CLOSE(IN)
      SUM = 0.D0
      DO I = 1, NS
        X(I) = X(I) +XTOL
        SUM = SUM +X(I)
      ENDDO
      DO I = 1, NS
        X(I) = X(I) /SUM
      ENDDO

      END SUBROUTINE READCOMPO
C-----------------------------------------------------------------------


