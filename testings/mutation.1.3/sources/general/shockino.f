      PROGRAM SHOCKINO
     
      IMPLICIT NONE   
      INCLUDE 'memorycr.cmn'
      INTEGER LRW, LIW
      PARAMETER (LRW= 900000, LIW = 1200)
      INTEGER I, J, NBIN, ITOL, ITASK, ISTATE, IOPT, MF, NEQ, 
     &        N, NB, IWORK(LRW), LFILE
      DOUBLE PRECISION  ATOL, MM, IDUM, JDUM, 
     &                  US, P1, T1, U2, P2, T2, G, RTOL, Z(1:NMAX),
     &                  NI1(1:NMAX), UE, XN, XN2, RHON,
     &                  RHON2, RHO, TMAX, TINI, TOUT, DT, 
     &                  ELEVELI(1:NMAX), RWORK(LRW)
      CHARACTER(100)  NAMEFILE
      EXTERNAL FUN, DUMMY 

C     CR model
      WRITE(*,*) "Enter CR model (1/2/3)"
      READ(*,*) CRFLAG
      IF (CRFLAG == 1) THEN
        DELTA = 7.D0
        NAMEFILE='HO'
        LFILE = 2
        WRITE(*,*) '=> Vibrational CR model (Trot=T)'
      ELSEIF (CRFLAG == 2) THEN
        DELTA = 5.D0        
        NAMEFILE='100Bin'
        LFILE = 6 
        WRITE(*,*) '=> Bin CR model'
      ELSEIF (CRFLAG == 3) THEN
        DELTA = 7.D0        
        NAMEFILE='Park'
        LFILE = 4 
        WRITE(*,*) '=> Park model'
        CRFLAG = 1
      ELSE
        WRITE(*,*) '=> CR model not implemented'
        STOP
      ENDIF
        WRITE(*,*) 

C     Specific heat ratio
      G = DELTA /(DELTA -2.D0)

C     Molar mass [kg / mole]     
      MN  = 14.0067D-3
      MN2 = 28.0134D-3

C     Formation enthalpy [J / mole] => [J / kg]
      HFN = 470818.D0 /MN 

C     Universal constants
      UKB  = 1.380658D-23   
      UNA  = 6.0221367D23   
      UR   = UKB *UNA
      UE   = 1.602191D-19
      UH   = 6.626075D-34
      UPI  = 3.14159265D0

C     N atom
C     -Degeneracy 
      GN = 4.D0
      GI(1) = GN
C     -Formation enthalpy [J]
      EI(1) = HFN *MN /UNA 

C     Rotational constants of N2
C     -Temperature [K]
      TROTN2 = 2.88D0 
C     -Symmetry number
      SYMN2  = 2.D0
      
C     Read NBIN levels: degeneracy and energy
C     The ground state energy EI(2) is zero
C     Conversion from [eV] to [J]
      OPEN(UNIT=1, 
     &     FILE='../data/chemistry/CR/'//NAMEFILE(1:LFILE)//'Lev.dat',
     &     STATUS='OLD')
        READ(1,*) NBIN
        NNBIN = NBIN
        DO I = 1, NBIN
          READ(1,*) GI(I+1), ELEVELI(I)
          EI(I+1) = (ELEVELI(I) -ELEVELI(1)) *UE 
        ENDDO
      CLOSE(1)

C     Molar mass [kg / mole] and specific heat [J/(kg K)]
      I = 1
        CPI(1) = 5.D0  *UR /(2.D0 *MN) 
        MI(1)  = MN
      DO I = 1, NBIN
        CPI(I+1) = DELTA  *UR /(2.D0 *MN2) 
        MI(I+1)  = MN2
      ENDDO

C     Reaction rates
C     -Dissociation
C      Conversion [cm^3 / s] -> [m^3 / (mole s)]
      OPEN(UNIT=1, 
     &     FILE='../data/chemistry/CR/'//NAMEFILE(1:LFILE)//'Dis.dat',
     &     STATUS='OLD')
        DO I = 1, NBIN
          READ(1,*) IDUM, AD(I), BD(I), CD(I)
            AD(I) = AD(I)*1D-6*UNA
        ENDDO
      CLOSE(1)
C     -Quasi-bound dissociation
      OPEN(UNIT=1, FILE='../data/chemistry/CR/Width_Bin.dat',
     &                  STATUS='OLD')
        DO I = 1, NBIN
          READ(1,*) IDUM, AS(I)
        ENDDO
      CLOSE(1)
C     -Excitation
C      Only the rates for excitation to higher levels are being used
C      Conversion [cm^3 / s] -> [m^3 / (mole s)]
      DO I = 1, NBIN
        DO J = 1, NBIN
          AE(I,J) = 0.D0
          BE(I,J) = 0.D0
          CE(I,J) = 0.D0
        ENDDO
      ENDDO 
        OPEN(UNIT=1, 
     &       FILE='../data/chemistry/CR/'//NAMEFILE(1:LFILE)//'Exc.dat',
     &       STATUS='OLD')
        DO I = 1, NBIN
          DO J = I+1, NBIN
            READ(1,*) IDUM, JDUM, AE(I,J), BE(I,J), CE(I,J)
            AE(I,J) = AE(I,J)*1.D-6*UNA
          ENDDO
        ENDDO
      CLOSE(1)

C     Number of equations (NBIN species mass, mixture momentum and 
C     total energy)
      NEQ = NBIN +3
      IF (NEQ > NMAX) THEN
        WRITE(*,*) 'Not enough memory!'
        STOP
      ENDIF

C     Free stream quantities
      P1 = 13.3D0 
      T1 = 4000.D0
      US = 7000.D0
      XN = 0.144740
      XN2   = 1.D0 -XN
      MM    = XN*MN +XN2*MN2
      RHON  = XN *MN *P1/(UR *T1)
      RHON2 = XN2 *MN2 *P1/(UR *T1)
      RHO   = RHON +RHON2
      MDOT  = RHO *US

C     Boltzmann distribution of the levels at T1
      CALL BOLTLEVELS (NBIN, UKB, T1, GI, EI, NI1)

C     Rankine-Hugoniot jump relations
      CALL  RHJUMP (UR/MM, G, US, P1, T1, U2, P2, T2)
      WRITE(*,*) '*Post-shock conditions (cold gas approximation)'
      WRITE(*,*) '------------------------------------------------'
      WRITE(*,*) ' -P2:                  ', P2/101325.D0,'[atm]'
      WRITE(*,*) ' -P2:                  ', P2,'[Pa]'
      WRITE(*,*) ' -T2:                  ', T2,'[K]'
      WRITE(*,*) ' -US-U2:               ', US-U2, '[m/s]' 

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
      MF   = 22
      Z(1) = RHON /RHO 
      DO I = 1, NBIN
        Z(I+1) = RHON2 *NI1(I) /RHO
      ENDDO
      Z(NBIN+2) = US -U2
      Z(NBIN+3) = T2
C     Loop      
      OPEN(UNIT=1,FILE='../output/CR.dat',STATUS='UNKNOWN')
      TMAX  = 100.D-1
      TINI  = 0.D0 
      TOUT  = TINI
      K = 0
      CALL WRITESOLCR(TOUT, NEQ, Z)
      N = 5; DT   = 10.D0**(-N)
C     Write solution every TOUT = 10^(-3)
C      NB = 10**(N-4)
      NB = 10**(N-2)
      I = 1; TOUT = TOUT +DT
      DOWHILE (TOUT < TMAX)
        CALL DLSODE (FUN, NEQ, Z, TINI, TOUT, ITOL, RTOL, ATOL, ITASK,
     &              ISTATE, IOPT, RWORK, LRW, IWORK, LIW, DUMMY, MF)
        IF (I-INT(I/NB)*NB == 0) THEN          
          CALL WRITESOLCR(TOUT, NEQ, Z)
        ENDIF
        I = I +1; TOUT = TOUT +DT     
        IF (ISTATE <= -2) THEN
          WRITE(*,*) 'ISTATE = ', ISTATE
        STOP
        ENDIF
      ENDDO 
      CLOSE(1)
      WRITE(*,*)
      IF (K>KMAX) THEN
         WRITE(*,*) 'K', K, ':', KMAX
         WRITE(*,*) 'Error: work arrays out of space.'
      ENDIF
 
      END PROGRAM SHOCKINO

C-----------------------------------------------------------------------
      SUBROUTINE RHJUMP (R, G, US, P1, T1, U2, P2, T2)
C-----------------------------------------------------------------------
C     This subroutine computes the variables after the shock based on 
C     Rankine-Hugoniot relations in the case of a cold gas.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION R, G, US, P1, T1, U2, P2, T2
C-----------------------------------------------------------------------
      DOUBLE PRECISION MS, GP1, C1, MS2
C-----------------------------------------------------------------------
      GP1 = G +1.D0

      C1  = SQRT(G * R *T1)
      MS  = US / C1
      MS2 = MS *MS
      T2 = T1 *(2.D0 *G *MS2 -G +1.D0)*(G -1.D0 +2.D0 /MS2) /(GP1 *GP1)

      U2 = C1 *2.D0 /GP1 *(MS - 1.D0 /MS)
      P2 = P1 *(2.D0 * G *MS2 -G +1.D0) /GP1

      END SUBROUTINE RHJUMP
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE  BOLTLEVELS (NBIN, KB, T, GI, EI, NI)
C-----------------------------------------------------------------------
C     This subroutine computes a Boltzmann population of the energy 
C     levels
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NBIN 
      DOUBLE PRECISION KB, T, GI(1:NBIN+1), EI(1:NBIN+1), NI(1:NBIN)
C-----------------------------------------------------------------------
      INTEGER I
      DOUBLE PRECISION Q 
C-----------------------------------------------------------------------
      Q = 0.D0
      DO I = 1, NBIN
        Q = Q + GI(I+1) *DEXP(-EI(I+1)/(KB *T))  
      ENDDO
      DO I = 1, NBIN
        NI(I) = GI(I+1) *DEXP(-EI(I+1)/(KB *T)) /Q
      ENDDO

      END SUBROUTINE BOLTLEVELS
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE WRITESOLCR (POSITION, NEQ, Z)  
C-----------------------------------------------------------------------
C     This subroutine writes the solution at the current position.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'memorycr.cmn'
C-----------------------------------------------------------------------
      INTEGER NEQ
      DOUBLE PRECISION POSITION, Z(1:NEQ)
C-----------------------------------------------------------------------
      INTEGER I
      DOUBLE PRECISION U, T, Y(1:NNBIN+1), RHO, X(1:NNBIN+1), MM
C-----------------------------------------------------------------------
      U = Z(NNBIN+2)
      T = Z(NNBIN+3)
      RHO = MDOT /U
      MM = 0.D0
      DO I = 1, NNBIN+1
        Y(I) = Z(I)
        MM = MM + Y(I) /MI(I)
      ENDDO
      MM = 1.D0 /MM
      DO I = 1, NNBIN+1
        X(I) = Y(I) *MM /MI(I)
      ENDDO


      K = K +1
      POSITIONK(K) = POSITION
      UK(K)        = U
      IF (K==1) THEN
        TIMEK(K)     = 0.D0
      ELSE
        TIMEK(K)     = TIMEK(K-1)+ 0.5D0 *(POSITIONK(K)-POSITIONK(K-1))
     &                 *(1.D0 /UK(K-1) +1.D0 /UK(K))
      ENDIF
      
      WRITE(1,2010) POSITION, U, T, (X(I), I=1,NNBIN+1) 
2010  FORMAT (200E14.6)

      END SUBROUTINE WRITESOLCR
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE FUN (NEQ, X, Z, ZP)  
C-----------------------------------------------------------------------
C     This subroutine computes the source terms of the ODEs (1D Euler 
C     equations). 
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'memorycr.cmn'
C-----------------------------------------------------------------------
      INTEGER NEQ
      DOUBLE PRECISION X, Z(1:NEQ), ZP(1:NEQ)
C-----------------------------------------------------------------------
      INTEGER I 
      DOUBLE PRECISION U, T, Y(1:NNBIN+1), OMEGA(1:NNBIN+1), OMEGASUM, 
     &                 YSUM, CPSUM, HSUM, A, B, C, D, E, F, RHO, 
     &                 H(1:NNBIN+1), DET
C-----------------------------------------------------------------------
      DO I = 1, NNBIN+1
        Y(I) = Z(I)
      ENDDO
      U = Z(NNBIN+2)
      T = Z(NNBIN+3)
      RHO = MDOT /U

      DO I = 1, NNBIN+1
        H(I) = EI(I) *UNA /MI(I) +CPI(I) *T  
      ENDDO
      CALL N2ARRHENIUS (T, Y, RHO, OMEGA)
C     Force the reaction rates to reach equilibrium
C      DO I = 1, NNBIN+1
C        OMEGA(I) = 1.D6 *OMEGA(I)  
C      ENDDO
      
      HSUM = 0.D0; OMEGASUM = 0.D0; CPSUM = 0.D0; YSUM = 0.D0
      DO I = 1, NNBIN+1
        OMEGASUM = OMEGASUM +OMEGA(I) 
        HSUM     = HSUM + H(I) *OMEGA(I) *MI(I)
        CPSUM    = CPSUM +Y(I)*CPI(I)
        YSUM     = YSUM  +Y(I)/MI(I)
      ENDDO

      A = MDOT *U /(UR* T) -RHO *YSUM 
      B = MDOT /T *YSUM
      C = -OMEGASUM
      D = MDOT *U
      E = MDOT *CPSUM
      F = -HSUM
      DET = A *E - B *D

      DO I = 1, NNBIN+1
        ZP(I)    = OMEGA(I) *MI(I) /MDOT  
      ENDDO
      ZP(NNBIN+2) = (C *E -B *F) /DET
      ZP(NNBIN+3) = (A *F -C *D) /DET

      END SUBROUTINE FUN
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE DUMMY()
      END SUBROUTINE DUMMY
C----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE N2ARRHENIUS (T, Y, RHO, OMEGA)
C-----------------------------------------------------------------------
C     This subroutine computes the chemical production rates
C     expressed in [mole m^-3 s^-1].
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'memorycr.cmn'
      DOUBLE PRECISION T, OMEGA(1:NNBIN+1), Y(1:NNBIN+1), RHO
      INTEGER I, J 
      DOUBLE PRECISION ND, KDEQ, KEEQ, LNT, KD, KE, NMOLE(1:NNBIN+1),
     &                 EXP1 
C-----------------------------------------------------------------------
      LNT = DLOG(T)
      DO I = 1, NNBIN+1
        NMOLE(I) = RHO * Y(I) /MI(I)
        OMEGA(I) = 0.D0
      ENDDO

      EXP1 = 0.D0
      DO I = 1, NNBIN +1
        EXP1 = EXP1 +NMOLE(I)
      ENDDO
      ND = EXP1 *UNA

C     Dissociation
      DO I = 1, NNBIN
        KD = AD(I)*DEXP(BD(I)*LNT-CD(I)/T) 
        CALL EQDISSOCIATION (I, ND, T, KDEQ)
        EXP1 = KD *NMOLE(1) *(NMOLE(I+1) -NMOLE(1)**2 /KDEQ)
        OMEGA(1) = OMEGA(1) +2.D0 *EXP1 
        OMEGA(I+1) = -EXP1
      ENDDO

C     Excitation
      DO I = 1, NNBIN-1
        DO J = I+1, NNBIN
          KE = AE(I,J)*DEXP(BE(I,J)*LNT-CE(I,J)/T) 
          CALL EQEXCITATION (I, J, T, KEEQ)
          EXP1 = KE *NMOLE(1) *(NMOLE(I+1) -NMOLE(J+1) /KEEQ)
          OMEGA(I+1) = OMEGA(I+1) -EXP1
          OMEGA(J+1) = OMEGA(J+1) +EXP1
        ENDDO
      ENDDO

      END SUBROUTINE N2ARRHENIUS
C----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE EQDISSOCIATION (I, ND, T, KEQ)
C-----------------------------------------------------------------------
C     This subroutine computes the equilibrium constant for the  
C     dissociation reaction at temperature T:
C     N2(I) + N = 3N 
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'memorycr.cmn'
      INTEGER I 
      DOUBLE PRECISION  ND, T, KEQ 
      DOUBLE PRECISION  QN, QN2, FAC1
C-----------------------------------------------------------------------
      FAC1 = 2.D0 *UPI *UR /(UH*UH *UNA*UNA)
      QN2  = (FAC1 *MI(2) *T )**1.5D0
      IF(CRFLAG == 1) THEN
        QN2 = QN2 *T /(TROTN2 *SYMN2) 
      ENDIF

      QN   = (FAC1 *MI(1) *T )**1.5D0
 
      KEQ = 1.D0 /UNA *(GI(1)*QN)**2 /(GI(I+1) *QN2)
     &     *DEXP(-(2.D0 *EI(1) - EI(I+1))/(UKB *T))     

      END SUBROUTINE EQDISSOCIATION
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE EQEXCITATION (I, J, T, KEQ)
C-----------------------------------------------------------------------
C     This subroutine computes the equilibrium constant for the  
C     excitation reaction at temperature T:
C     N2(I) + N = N2(J) + N 
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'memorycr.cmn'
      INTEGER I, J
      DOUBLE PRECISION  T, KEQ
C-----------------------------------------------------------------------
      KEQ = GI(J+1)/GI(I+1) * EXP(-(EI(J+1) -EI(I+1)) / (UKB*T) )

      END SUBROUTINE EQEXCITATION
C-----------------------------------------------------------------------
