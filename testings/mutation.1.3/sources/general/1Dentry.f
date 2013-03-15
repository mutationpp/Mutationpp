
C     1DEntry code

C     Thierry E. Magin
C     Bruno Dias

C     von Karman Institute, Rhode-Saint-Genese (Belgium), 01-08-2012


C-----------------------------------------------------------------------
      PROGRAM ONEDENTRY 
C-----------------------------------------------------------------------
C     Starting from free stream conditions, this program computes the 
C     stagnation point heat flux and pressure in an atmospheric entry 
C     flow. Notice that the free stream flow can be in vibrational 
C     nonequilibrium in the case of a ground simulation in a high 
C     enthalpy facility. The equivalent reservoir pressure is also
C     computed.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE '../general/memory.cmn'
      INCLUDE '../general/memoryshock.cmn'
C-----------------------------------------------------------------------
      INTEGER LRW, LIW, N, NB, IWRITE, C
      PARAMETER (LRW= 900000, LIW = 1200)
      INTEGER  NEQ, IS, IMOD, ITERMAX, ITER, INDX(1:3), I, IV, ICUT, MF,
     &         ITOL, ITASK, ISTATE, IOPT, IWORK(LRW), LCHAR, IC1, IC2,
     &         LSPECIES
      DOUBLE PRECISION X1(1:NMAX), X2(1:NMAX), XN(1:NMAX), XINI(1:NMAX),
     &                 P1, T1, US, GAMMA, V2, EPS, YINI(1:NMAX),
     &                 P2, T2, U2, RESMIN, RES, RESINI, R0(1:3), 
     &                 J(1:3,1:3), D, MM1, Y2(1:NMAX), X1TOL(1:NMAX),
     &                 ZP(NMAX+3), RHO, TINI, TOUT, TMAX, TOL, 
     &                 RHO2, RHOV2, RHOE2, TVIB2(1:NMAX), ALPHA,
     &                 RHO1, RHOV1, RHOE1, Y1(1:NMAX), TVIB1(1:NMAX),
     &                 EM1(1:NMAX), EM2(1:NMAX), EM3(1:NMAX), MM2,
     &                 EM4(1:NMAX), EM5(1:NMAX), EM6(1:NMAX), 
     &                 MIXH1, MIXH2, MIXH3, MIXH4, MIXH5, MIXH6,
     &                 MIXH21, MIXH22, MIXH23, MIXH24, MIXH25, MIXH26,
     &                 MIXS1, MIXS2, MIXS3, MIXS4, MIXS5, MIXS6, 
     &                 ND1, RHOI2(1:NMAX), RHOET, PT2, TT2, PT1,
     &                 GAMMAI1, GAMMAE1, DRHODP, SOUNDE1, MACH2, PT2IS,
     &                 PT1IS, MACH1, TT1, MM, DUM, ENTH1, ENTR1, ETA1,
     &                 GAMMAI2, GAMMAE2, TV1, SOUNDE2, ND2, ENTH2, 
     &                 TT2IS, TT1IS
      CHARACTER(10) MIXTURE, REACTION, TRANSFER      
      CHARACTER(100) PATH
      EXTERNAL ENTHALPY
C-----------------------------------------------------------------------
      EPS = 1.D-2
      TOL = 1.D-16
      PATH = '..'
      REACTION = 'empty'
      TRANSFER = 'empty'
      NVIB = 0
      NELEQ = 0

C     Free stream conditions
      MIXTURE = 'nitrogen5'
C      MIXTURE = 'air11'
      GAMMA = 1.4D0
      P1 = 229.D0
      T1 = 41.D0
      TV1 = 1500.D0
      TV1 = T1 
      US = 1850.D0
C     Laboratory reference frame: U1 = 0.D0
C     Shock reference frame: V1 = US

      WRITE(*,*) '*Free stream conditions'
      WRITE(*,*) '-----------------------'
      WRITE(*,*) ' -P1:                  ', P1/101325.,'[atm]'
      WRITE(*,*) ' -P1:                  ', P1,'[Pa]'
      WRITE(*,*) ' -T1:                  ', T1,'[K]'
      WRITE(*,*) ' -TV1:                 ', TV1,'[K]'
      WRITE(*,*) ' -U1:                  ', US, '[m/s]' 
      WRITE(*,*)


C     Link with the MUTATION library
C      flag IMOD = 0: only thermodynamic and chemistry routines
C                = 1: thermodynamic, chemistry and transport routines
C     ---------------------------------------------------------------
      IMOD = 1 

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

C     Initialize data independent on temperature and pressure
      CALL INITIALIZE (PATH, MIXTURE, REACTION, TRANSFER, WR1, LWR1, 
     &                 WR3, LWR3, WR4, LWR4, WI, LWI, WC, LWC, IMOD)

C     Initial composition
C     -------------------
      CALL NUCLEAR (WR1, LWR1, XN)
      DO IS = 1, NS
        XINI(IS) = 1.D0
      ENDDO
      CALL COMPOSITION (WR1, LWR1, WI, LWI, T1, P1, XN, XINI, X1)

C     Conservative variables
C     State 1
      CALL MOLE2MASSFRAC (WR1, LWR1, X1, Y1)
      CALL NUMBERD (WR1, LWR1, P1, T1, TV1, X1, ND1)
      CALL DENSITY (WR1, LWR1, X1, ND1, RHO1)
      CALL MOLARMASS (WR1, LWR1, RHO1, ND1, MM1)
      RHOV1 = RHO1 *US
      DO I = 1, NV
        TVIB1(I) = TV1
      ENDDO
      CALL ENERGYMASS (WR1, LWR1, WI, LWI, T1, TV1, T1, TVIB1, P1, EM1,
     &                 EM2, EM3, EM4, EM5, EM6)
      RHOE1 = 0.D0
      DO I = 1, NS
        RHOE1 = RHOE1 +Y1(I)*EM1(I)
      ENDDO
      RHOE1 = RHOE1 *RHO1 +0.5D0*RHOV1*RHOV1/RHO1

C     State 2
C     Rankine-Hugoniot's jump relations for the normal shock 
C     -------------------------------------------------------
C     Initialization (cold gas)
      CALL RHCOLD (WR1(IUR)/MM1, GAMMA, US, P1, T1, V2, P2, T2)
      U2 = US -V2
      CALL COMPOSITION (WR1, LWR1, WI, LWI, T2, P2, XN, XINI, X2)
      CALL MOLE2MASSFRAC (WR1, LWR1, X2, Y2)
      CALL NUMBERD (WR1, LWR1, P2, T2, T2, X2, ND2)
      CALL DENSITY (WR1, LWR1, X2, ND2, RHO2)
      RHOV2 = RHO2 *V2
      DO I = 1, NV
        TVIB2(I) = T2
      ENDDO
      CALL ENERGYMASS (WR1, LWR1, WI, LWI, T2, T2, T2, TVIB2, P2, EM1,
     &                 EM2, EM3, EM4, EM5, EM6)
      RHOE2 = 0.D0
      DO I = 1, NS
        RHOE2 = RHOE2 +Y2(I)*EM1(I)
      ENDDO
      RHOE2 = RHOE2 *RHO2 +0.5D0*RHOV2*RHOV2/RHO2
      WRITE(*,*) '*Post-shock conditions (cold gas approximation)'
      WRITE(*,*) '-----------------------------'
      WRITE(*,*) ' -P2:                  ', P2/101325.,'[atm]'
      WRITE(*,*) ' -P2:                  ', P2,'[Pa]'
      WRITE(*,*) ' -T2:                  ', T2,'[K]'
      WRITE(*,*) ' -US:                  ', US, '[m/s]' 
      WRITE(*,*) ' -V2:                  ', V2, '[m/s]' 
      WRITE(*,*) ' -U2=US-V2:            ', U2, '[m/s]' 
      WRITE(*,*)

      TINI = T2
      CALL MOLE2MASSFRAC (WR1, LWR1, X2, YINI)
      CALL SYSTEMFAD (WR1, LWR1, WI, LWI,  RHO1, RHOV1, RHOE1, P1,
     &                RHO2, RHOV2, RHOE2, XN, YINI, TINI, T2, Y2, R0)
      CALL JACOBIAN (WR1, LWR1, WI, LWI,  RHO1, RHOV1, RHOE1, P1,
     &               RHO2, RHOV2, RHOE2, XN, YINI, TINI, J)

      RES = R0(1) *R0(1) +R0(2) *R0(2) +R0(3) *R0(3)
      RESINI = RES
      RES = 1.D0

C     Test (Newton's method, residual and maximum number of iterations)
      RESMIN = 1.D-12; ITERMAX = 100; ITER = 0
C     Underrelaxation parameter ALPHA \in ]0, 1]
      ALPHA = 1.D0
      RESMIN = RESMIN *RESMIN

      DO WHILE ( (ITER < ITERMAX) .AND.  (RES > RESMIN) )
        R0(1) = -R0(1)
        R0(2) = -R0(2)
        R0(3) = -R0(3)
        CALL LUDCMP (J, 3, 3, INDX, D)
        CALL LUBKSB (J, 3, 3, INDX, R0)

        RHO2  = RHO2  +ALPHA *R0(1)
        RHOV2 = RHOV2 +ALPHA *R0(2)
        RHOE2 = RHOE2 +ALPHA *R0(3)

        CALL SYSTEMFAD (WR1, LWR1, WI, LWI,  RHO1, RHOV1, RHOE1, P1,
     &                  RHO2, RHOV2, RHOE2, XN, YINI, TINI, T2, Y2, R0)
        CALL JACOBIAN (WR1, LWR1, WI, LWI,  RHO1, RHOV1, RHOE1, P1,
     &                  RHO2, RHOV2, RHOE2, XN, YINI, TINI, J)
        RES  = (R0(1) *R0(1) +R0(2) *R0(2) +R0(3) *R0(3)) /RESINI
        ITER = ITER +1                 
      ENDDO

      V2 = RHOV2 /RHO2
      U2 = US -V2
      CALL MASS2PRESSURE (WR1, LWR1, RHO2, T2, T2, Y2, P2)


C     Stagnation point pressure
C     -------------------------
      CALL COMPOSITION (WR1, LWR1, WI, LWI, T2, P2, XN, XINI, X2)
      CALL NUMBERD (WR1, LWR1, P2, T2, T2, X2, ND2)
      CALL MOLARMASS (WR1, LWR1, RHO2, ND2, MM2)
      DO I = 1, NV
        TVIB2(I) = T2
      ENDDO
      CALL MIXPROPERTY (WR1, LWR1, WI, LWI, T2, T2, T2, TVIB2, DUM, X2,
     &                  MIXH21, MIXH22, MIXH23, MIXH24, MIXH25, MIXH26, 
     &                  ENTHALPY)
      ENTH2 = MIXH21 /MM2 +0.5D0 *V2*V2
      CALL SAINTVENANTPT (WR1, LWR1, WI, LWI, T2, P2, V2, XN, X2, EPS,
     &                    TT2IS, PT2IS, GAMMAE2, GAMMAI2, SOUNDE2)

      CALL TOTAL (WR1, LWR1, WI, LWI, T2, P2, V2, XN, X2, EPS, TT2, 
     &            PT2)

C     Heat flux (Fay-Ridell)
C     ----------------------
      CALL COMPOTOL (X1, TOL, X1TOL)
      CALL COLLISION (WR1, LWR1, WR2, LWR2, T1, T1, ND1, X1)
      CALL ETACG (WR1, LWR1, WR2, LWR2, X1TOL, ETA1)

C     Reservoir pressure
C     ------------------
C     Entropy
      CALL MIXENTROPY (WR1, LWR1, WI, LWI, T1, TV1, T1, TVIB1, P1, X1,
     &                 MIXS1, MIXS2, MIXS3, MIXS4, MIXS5, MIXS6)
      ENTR1 = MIXS1 /MM1

C     Total enthalpy
C     -Mixture enthalpy per unit mole
      CALL MIXPROPERTY (WR1, LWR1, WI, LWI, T1, TV1, T1, TVIB1, DUM, X1,
     &                  MIXH1, MIXH2, MIXH3, MIXH4, MIXH5, MIXH6, 
     &                  ENTHALPY)
C     -Total enthalpy per unit mass
      ENTH1 = MIXH1 /MM1 +0.5D0 *US*US

C     Perfect gas assumption 
      WRITE(*,*)   
      WRITE(*,*)     "ADAPT T1, Tv1 to single T0 for reservoir!"
      WRITE(*,*)   
      CALL SAINTVENANTPT (WR1, LWR1, WI, LWI, T1, P1, US, XN, X1, EPS,
     &                    TT1IS, PT1IS, GAMMAE1, GAMMAI1, SOUNDE1)
      CALL TOTAL (WR1, LWR1, WI, LWI, T1, P1, US, XN, X1, EPS, TT1, 
     &            PT1)


      WRITE(*,*) '*Post-shock conditions (high temperatures)'
      WRITE(*,*) '-------------------------------------------'
      WRITE(*,*) ' -Number of Newton iterations:', ITER
      WRITE(*,*) ' -Residual:            ', DSQRT(RES)
      WRITE(*,*) ' -P2:                  ', P2/101325.,'[atm]'
      WRITE(*,*) ' -P2:                  ', P2,'[Pa]'
      WRITE(*,*) ' -T2:                  ', T2,'[K]'
      WRITE(*,*) ' -US:                  ', US, '[m/s]' 
      WRITE(*,*) ' -V2:                  ', V2, '[m/s]' 
      WRITE(*,*) ' -U2=US-V2:            ', U2, '[m/s]'
      WRITE(*,*) ' -RHO2:                ', RHO2,'[kg/m3]'
      WRITE(*,*) ' -RHO2*V2:             ', RHOV2,'[kg/(m2 s)]'
      WRITE(*,*) ' -RHO2*E2:             ', RHOE2,'[J/m3]'
      WRITE(*,*) ' -PT2 Saint Venant:    ', PT2IS,'[Pa]'
      WRITE(*,*) ' -PT2:                 ', PT2,'[Pa]'
      WRITE(*,*) ' -TT2 Saint Venant:    ', TT2IS,'[Pa]'
      WRITE(*,*) ' -TT2:                 ', TT2,'[K]'
      WRITE(*,*) 'MM2                    ', MM2
      WRITE(*,*) 'CP COLD',  ((MIXH21 /MM2) /T2)
      WRITE(*,*) 'TTOTAL COLD',  ENTH2/((MIXH21 /MM2) /T2)
      WRITE(*,*) 'ENTH2', ENTH2
      WRITE(*,*) ' -GAMMAE2:             ', GAMMAE2,'[-]'
      WRITE(*,*) ' -GAMMAI2:             ', GAMMAI2,'[-]'
      WRITE(*,*) ' -PT1 Saint Venant:    ', PT1IS,'[Pa]'
      WRITE(*,*) ' -PT1:                 ', PT1,'[Pa]'
      WRITE(*,*) ' -TT1 Saint Venant:    ', TT1IS,'[Pa]'
      WRITE(*,*) ' -TT1:                 ', TT1,'[K]'
      WRITE(*,*) 'MM1                    ', MM1
      WRITE(*,*) 'CP COLD',  ((MIXH1 /MM1) /T1)
      WRITE(*,*) 'TTOTAL COLD',  ENTH1/((MIXH1 /MM1) /T1)
      WRITE(*,*) 'ENTH1', ENTH1
      WRITE(*,*) ' -GAMMAI1:             ', GAMMAI1,'[-]'
      WRITE(*,*) ' -GAMMAE1:             ', GAMMAE1,'[-]'
      WRITE(*,*) ' -R2 perfect gas:      ', P2/(RHO2*T2)
      WRITE(*,*) ' -R1 perfect gas:      ', P1/(RHO1*T1)
      WRITE(*,*) ' -ETA1:                ', ETA1,'[kg/(m s)]'

      CALL PRESENTATION (2)

      END PROGRAM ONEDENTRY 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE RHCOLD (R, G, U1, P1, T1, U2, P2, T2)
C-----------------------------------------------------------------------
C     This subroutine computes the variables after the shock based on 
C     Rankine-Hugoniot relations in the case of a cold gas. 
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION R, G, U1, P1, T1, U2, P2, T2
C-----------------------------------------------------------------------
      DOUBLE PRECISION MS, GP1, C1, MS2, A, B, C
C-----------------------------------------------------------------------
      GP1 = G +1.D0

      C1  = SQRT(G * R *T1)
      MS  = U1 / C1
      MS2 = MS *MS

      T2 = T1 *(2.D0 *G *MS2 -G +1.D0)*(G -1.D0 +2.D0 /MS2) /(GP1 *GP1)
      U2 = ((G -1.D0)*MS2 +2.D0) *U1 /(GP1 *MS2) 
      P2 = P1 *(2.D0 * G *MS2 -G +1.D0) /GP1

      END SUBROUTINE RHCOLD
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SYSTEMFAD (WR1, LWR1, WI, LWI,  RHO1, RHOV1, RHOE1, P1,
     &                      RHO2, RHOV2, RHOE2, XN, YINI, TINI, T2, Y2,
     &                      R)
C-----------------------------------------------------------------------
C    This subroutines computes the residual for the Rankine-Hugoniot
C    relations in thermo-chemical equilibrium.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      INTEGER LWR1, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(LWR1), RHO1, RHOV1, RHOE1, P1, RHO2, RHOV2, 
     &                 RHOE2, XN(1:NC), YINI(1:NS), TINI, R(1:3), T2,
     &                 Y2(1:NS)
C-----------------------------------------------------------------------
      INTEGER IS 
      DOUBLE PRECISION V1, V2, P2, RHOET
C-----------------------------------------------------------------------
      V2 = RHOV2 /RHO2
      V1 = RHOV1 /RHO1

      RHOET = RHOE2 - 0.5*RHOV2*V2
      CALL TCEQNEWTON (WR1, LWR1, WI, LWI, RHOET, RHO2, XN, YINI, TINI,
     &                 T2, Y2)
      CALL MASS2PRESSURE (WR1, LWR1, RHO2, T2, T2, Y2, P2)
      
      R(1) = RHOV2- RHOV1 
      R(2) = RHOV2*V2+P2-RHOV1*V1-P1
      R(3) = (RHOE2+P2)*V2-(RHOE1+P1)*V1
      
      END SUBROUTINE SYSTEMFAD
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE JACOBIAN (WR1, LWR1, WI, LWI,  RHO1, RHOV1, RHOE1, 
     &                     P1, RHO2, RHOV2, RHOE2, XN, YINI, TINI, J)
C-----------------------------------------------------------------------
C    This subroutines computes the Jacobian matrix J for the Rankine-
C    Hugoniot relations in thermo-chemical equilibrium.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      INTEGER LWR1, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(LWR1), RHO1, RHOV1, RHOE1, P1, RHO2, RHOV2, 
     &                 RHOE2, XN(1:NC), YINI(1:NS), TINI, J(1:3,1:3) 
C-----------------------------------------------------------------------
      INTEGER I
      DOUBLE PRECISION R(1:3), R0(1:3), EPS, EPSP1, T2, Y2(1:NS) 
C-----------------------------------------------------------------------
      EPS = 1.D-2; EPSP1 = 1.D0 +EPS

      CALL SYSTEMFAD (WR1, LWR1, WI, LWI,  RHO1, RHOV1, RHOE1, P1,
     &                RHO2, RHOV2, RHOE2,  XN, YINI, TINI, T2, Y2, R0)

      CALL SYSTEMFAD (WR1, LWR1, WI, LWI,  RHO1, RHOV1, RHOE1, P1,
     &                RHO2*EPSP1, RHOV2, RHOE2, XN, YINI, TINI, T2, Y2, 
     &                R)
      DO I = 1, 3
        J(I,1) = (R(I)-R0(I))/(RHO2  *EPS)
      ENDDO
      CALL SYSTEMFAD (WR1, LWR1, WI, LWI,  RHO1, RHOV1, RHOE1, P1,
     &                RHO2, RHOV2*EPSP1, RHOE2,  XN, YINI, TINI, T2, Y2,
     &                R)
      DO I = 1, 3
        J(I,2) = (R(I)-R0(I))/(RHOV2  *EPS)
      ENDDO
      CALL SYSTEMFAD (WR1, LWR1, WI, LWI,  RHO1, RHOV1, RHOE1, P1,
     &                RHO2, RHOV2, RHOE2*EPSP1,  XN, YINI, TINI, T2, Y2,
     &                R)
      DO I = 1, 3
        J(I,3) = (R(I)-R0(I))/(RHOE2  *EPS)
      ENDDO
     
      END SUBROUTINE JACOBIAN
C-----------------------------------------------------------------------

