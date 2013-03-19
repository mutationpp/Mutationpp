C-----------------------------------------------------------------------
      SUBROUTINE ARRHENIUS (WR1, LWR1, WR3, LWR3, WI, LWI, Y, YTOL, P, 
     &                      TEF, RHO, OMEGA)
C-----------------------------------------------------------------------
C     This subroutine computes the reaction rates in thermal equilibrium
C     or nonequilibrium situations. The Arrhenius reaction rates are
C     expressed in [mole m^-3 s^-1]. 
C     TEF(1)        : TH = translational temperature of heavy particles.
C     TEF(2:1+NVIB) : Effective vibrational temperatures (if NVIB /= 0).
C     TEF(2+NVIB)   : TE = translational temperature of electrons.
C     -Remark 1. The vibrational modes of a given polyatomic molecule 
C      are assumed to be distributed at the same temperature. Hence,
C      NVIB <= number of molecules.
C     -Remark 2. The vibrational temperature TV of the species IS is 
C      found by means of TV(IS)= TEF(WI(IVIBTEMPI+IS-1)). 
C      If the species IS is not a molecule, WI(IVIBTEMPI+IS-1) = 1.
C     -Remark 3. Why effective temperatures are used?
C         -Reduce the computational cost. In most applications, NVIB is 
C          close to 1. Therefore, the number of evaluations of logarithm
C          and square root functions is reduced in the computation of
C          the rates of dissociation by atomic or molecular impact.
C         -Allow for a more general choice of TV. For instance, TV can 
C          also be the translational temperature of either the free 
C          electrons or the heavy particles.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      INTEGER  LWR1, LWR3, LWI, WI(1:LWI), J
      DOUBLE PRECISION WR1(1:LWR1), WR3(1:LWR3), Y(1:NS), YTOL, P, 
     &                 TEF(1:NVIB+2), RHO, OMEGA(1:NS)
C-----------------------------------------------------------------------
      DOUBLE PRECISION XI(1:NREA), KF(1:NREA), KB(1:NREA), KEQ, KFEQ,
     &                 GTH(1:NS), GTE(1:NS), YT(1:NS), THV(1:NV+1), 
     &                 TEV(1:NV+1), PRODP, PRODR, MOLE(1:NS), 
     &                 LNMOLE(1:NS), TH, TV(1:NV), TE, SUM, LNR,
     &                 RTH, LNTH, LNRTH, RTE, LNTE,  LNRTE, 
     &                 STHTV(1:NS), LNSTHTV(1:NS), TVIB(1:NS), 
     &                 SE(1:NVIB+2), LNSE(1:NVIB+2),
     &                 UTE, TRATIO, Ia, KeV, TREFO2, TREFN2, TREFNO
      INTEGER IS, IR, IV, ISDIS, NVIBMODE, IC, JR
C-----------------------------------------------------------------------
C     Translational temperature of heavy particles
      TH = TEF(1)
!      IF(TH.GT.30000.D0) THEN
!           TH = 30000.D0
!       ENDIF
C     Translational temperature of electrons
      TE = TEF(NVIB+2)
      IF(TE.GT.15000.D0) THEN
          TE = 15000.D0
      ENDIF
C     Vibrational temperatures
      IC = 0
      DO IS = 1, NS
        NVIBMODE = WI(IVIBI+IS-1)
        DO IV = 1, NVIBMODE
          IC = IC +1
          TV(IC)   = TEF(WI(IVIBTEMPI+IS-1))
!          IF(TV(IC).GT.15000.D0) THEN
!              TV(IC)= 15000.D0
!          ENDIF
        ENDDO
      ENDDO

C     Translational temperature vectors of size (1:NV)
      DO IV = 1, NV
        THV(IV) = TH
        TEV(IV) = TE
      ENDDO
      
C     Logarithms
      LNR    = DLOG(WR1(IUR))
      RTH    = WR1(IUR) *TH
      LNTH   = DLOG(TH)
      LNRTH  = LNR +LNTH
      RTE    = WR1(IUR) *TE
      LNTE   = DLOG(TE)
      LNRTE  = LNR +LNTE
      SE(1)   = TH
      LNSE(1) = LNTH
      DO IV = 1, NVIB
        SE(IV+1)   = DSQRT(TH *TEF(IV+1))
        LNSE(IV+1) = DLOG(SE(IV+1))
      ENDDO
      SE(NVIB+2)   = DSQRT(TH *TE)
      LNSE(NVIB+2) = DLOG(SE(NVIB+2))
      DO IS = 1, NS
        STHTV(IS)   = SE(WI(IVIBTEMPI+IS-1))
        LNSTHTV(IS) = LNSE(WI(IVIBTEMPI+IS-1))
      ENDDO

C     Tolerance on mass fractions
      DO IS = 1, NS
        IF (Y(IS) < YTOL) THEN
          YT(IS) = YTOL
        ELSE
          YT(IS) = Y(IS)
        ENDIF
      ENDDO

C     Gibbs's free energy (per unit mole) at p = 1.D0 Pa
      CALL GIBBS (WR1, LWR1, WI, LWI, TH, TH, TH, THV, 1.D0, GTH)
      CALL GIBBS (WR1, LWR1, WI, LWI, TE, TE, TE, TEV, 1.D0, GTE)
      
C     Molar density [mole /m^3]
      DO IS = 1, NS
        MOLE(IS)   = RHO *YT(IS) /WR1(IMI+IS-1)
        LNMOLE(IS) = LOG(MOLE(IS))
      ENDDO

C     Reaction rates
      DO IR = 1, NREA
        SELECT CASE(WI(ITNEQ+IR-1))
C     -1. Dissociation by atomic or molecular impact / recombination 
C         => KF(SQRT(TH*TV)), KB(TH)
          CASE(1)
            ISDIS = WI(IDISSOCIATION+IR-1)
            KF(IR) = WR3(IRATE+(IR-1)*4) 
     &              *DEXP( WR3(IRATE+(IR-1)*4+1) *LNSTHTV(ISDIS)
     &              -WR3(IRATE+(IR-1)*4+2) /STHTV(ISDIS))
            KFEQ   = WR3(IRATE+(IR-1)*4) 
     &              *DEXP( WR3(IRATE+(IR-1)*4+1) *LNTH
     &              -WR3(IRATE+(IR-1)*4+2) /TH)
            CALL KEQUILIBRIUM (NS, NREA, GTH, WR3(ISTOIR), WR3(ISTOIP), 
     &                         WR1(IMI), RTH, LNRTH, KEQ, IR)
            KB(IR) = KFEQ /KEQ
C            WRITE(*,*)'c', KF(IR)
C            STOP
C     -2. Dissociation by electron impact / recombination
C     -3. Ionization by electron impact / recombination
C         => KF(TE), KB(TE)
          CASE(2,3)
            KF(IR) = WR3(IRATE+(IR-1)*4) 
     &              *DEXP( WR3(IRATE+(IR-1)*4+1) *LNTE
     &              -WR3(IRATE+(IR-1)*4+2) /TE)
            KFEQ   = KF(IR)
            CALL KEQUILIBRIUM (NS, NREA, GTE, WR3(ISTOIR), WR3(ISTOIP), 
     &                         WR1(IMI), RTE, LNRTE, KEQ, IR)
            KB(IR) = KFEQ /KEQ
C     -4. Associative ionization
C         => KF(TH), KB(TE)
          CASE(4)
            KF(IR) = WR3(IRATE+(IR-1)*4) 
     &              *DEXP( WR3(IRATE+(IR-1)*4+1) *LNTH
     &              -WR3(IRATE+(IR-1)*4+2) /TH)
            KFEQ   = WR3(IRATE+(IR-1)*4) 
     &              *DEXP( WR3(IRATE+(IR-1)*4+1) *LNTE
     &              -WR3(IRATE+(IR-1)*4+2) /TE)
            CALL KEQUILIBRIUM (NS, NREA, GTE, WR3(ISTOIR), WR3(ISTOIP), 
     &                         WR1(IMI), RTE, LNRTE, KEQ, IR)
            KB(IR) = KFEQ /KEQ

C     -7. Dissociative recombination 
C         => KF(TE), KB(TH)
          CASE(7)
            KF(IR) = WR3(IRATE+(IR-1)*4) 
     &              *DEXP( WR3(IRATE+(IR-1)*4+1) *LNTE
     &              -WR3(IRATE+(IR-1)*4+2) /TE)
            KFEQ   = WR3(IRATE+(IR-1)*4) 
     &              *DEXP( WR3(IRATE+(IR-1)*4+1) *LNTH
     &              -WR3(IRATE+(IR-1)*4+2) /TH)
            CALL KEQUILIBRIUM (NS, NREA, GTH, WR3(ISTOIR), WR3(ISTOIP), 
     &                         WR1(IMI), RTH, LNRTH, KEQ, IR)
            KB(IR) = KFEQ /KEQ

C     -5. Radical reactions (including Zeldovich reactions)
C     -6. Charge exchange
C         => KF(TH), KB(TH)
          CASE(5,6)
            KF(IR) = WR3(IRATE+(IR-1)*4) 
     &              *DEXP( WR3(IRATE+(IR-1)*4+1) *LNTH
     &              -WR3(IRATE+(IR-1)*4+2) /TH)
            KFEQ   = KF(IR)
            CALL KEQUILIBRIUM (NS, NREA, GTH, WR3(ISTOIR),  WR3(ISTOIP), 
     &                         WR1(IMI), RTH, LNRTH, KEQ, IR)
            KB(IR) = KFEQ /KEQ

C         Electron Impact excitation reactions for Nitrogen and Oxygen
          CASE(8)  ! Allowed Transitions (l_i /= l_j)
           ! IF (TE.GT.15000.D0) THEN
           !       TE = 15000.D0
           ! ENDIF
!      TE = 16000.D0
!
!      OPEN(UNIT=2,FILE='../output/RatesArnauda.dat',STATUS='UNKNOWN')
!      DO WHILE (TE.LT.14000.D0) 
            UTE= DSQRT((8.D0*WR1(IUKB)*TE)/(WR1(IUPI)
     &                             *WR1(IMI)/WR1(IUNA)))

            TRATIO = (WR3(IRATE+(IR-1)*4+2) - WR3(IRATE+(IR-1)*4+1))*
     &              1.43876866D0/TE

            Ia = 0.63255D0 * (TRATIO**(-1.6454D0)) * DEXP(-TRATIO) 

            KF(IR) = (UTE*4.D0*WR1(IUPI)*((0.529D-10)**2.D0)
     &               *0.05D0*(157820.40447D0/TE)**2.D0 *Ia)*WR1(IUNA)

            KFEQ   = KF(IR)
!      WRITE(2,*)TE,KF(IR)*1.D6/WR1(IUNA)
!      TE = TE + 100.D0
!      ENDDO
!      STOP
            CALL KEQUILIBRIUM (NS, NREA, GTE, WR3(ISTOIR), WR3(ISTOIP), 
     &                         WR1(IMI), RTE, LNRTE, KEQ, IR)
            KB(IR) = KFEQ /KEQ
!            WRITE(*,*)'BACKWORD',KB(IR)
!            KB(IR) = 9.D0 * KFEQ* DEXP(TRATIO)
!            WRITE(*,*)'BACKWORD',KB(IR)
!            STOP

          CASE(9)  ! Optically Forbidden Transitions (l_i = l_j)
           ! IF (TE.GT.15000.D0) THEN
           !       TE = 15000.D0
           ! ENDIF

!      TE = 6000.D0
!      OPEN(UNIT=2,FILE='../output/RatesArnaudF.dat',STATUS='UNKNOWN')
!      DO WHILE (TE.LT.14000.D0) 
            UTE= DSQRT((8.D0*WR1(IUKB)*TE)/(WR1(IUPI)
     &                             *WR1(IMI)/WR1(IUNA)))

            TRATIO = (WR3(IRATE+(IR-1)*4+2) - WR3(IRATE+(IR-1)*4+1))*
     &              1.43876866D0/TE

            Ia = 0.23933D0 * (TRATIO**(-1.4933D0)) * DEXP(-TRATIO) 

            KF(IR) = (UTE*4.D0*WR1(IUPI)*((0.529D-10)**2.D0)
     &               *0.05D0*(TRATIO)**2.D0 * Ia)*WR1(IUNA)
!      WRITE(2,*)TE,KF(IR)*1.D6/WR1(IUNA)
!      TE = TE + 100.D0
!      ENDDO
!      STOP

            KFEQ   = KF(IR)
            CALL KEQUILIBRIUM (NS, NREA, GTE, WR3(ISTOIR), WR3(ISTOIP), 
     &                         WR1(IMI), RTE, LNRTE, KEQ, IR)

            KB(IR) = KFEQ /KEQ

C         Electron Impact Ionisation reactions for Nitrogen and Oxygen


          CASE(10)  !Electron Impact Ionisation
!      TE = 6000.D0
!      OPEN(UNIT=2,FILE='../output/RatesArnaudF.dat',STATUS='UNKNOWN')
!      DO WHILE (TE.LT.14000.D0) 
            !IF (TE.GT.15000.D0) THEN
            !      TE = 15000.D0
            !ENDIF
            UTE= DSQRT((8.D0*WR1(IUKB)*TE)/(WR1(IUPI)
     &                             *WR1(IMI)/WR1(IUNA)))


          TRATIO = (WR3(IRATE+(IR-1)*4+2) - WR3(IRATE+(IR-1)*4+1))*
     &              1.43876866D0/TE

            Ia = 0.63255D0 * (TRATIO**(-1.6454D0)) * DEXP(-TRATIO) 

            KF(IR) = (UTE*4.D0*WR1(IUPI)*((0.529D-10)**2.D0)
     &              *(157820.40447D0/TE)**2.D0 * Ia)*WR1(IUNA)

            KFEQ   = KF(IR)
!      WRITE(2,*)TE,KF(IR)*1.D6/WR1(IUNA)
!      TE = TE + 100.D0
!      ENDDO
!      STOP

            CALL KEQUILIBRIUM (NS, NREA, GTE, WR3(ISTOIR), WR3(ISTOIP), 
     &                         WR1(IMI), RTE, LNRTE, KEQ, IR)
            KB(IR) = KFEQ /KEQ

          CASE(11)  ! Spontaneus Emission !! sec^-1 mole
           
            KF(IR) = WR3(IRATE+(IR-1)*3+1) !*WR1(IUNA)
            KB(IR) = 1.D-16
            !WRITE(*,*) 'Mechanism not implemented YET ! '  
            !STOP

C          ARNAUD's Model
          CASE(20)
C            WARNING ! This CASE has not been TESTED !! Usually is used only for 
C                      for negative IONS.. Watch out !! Arnaud put + instead of minus
C                      in the formula .. second term !!                   
             KF(IR) = WR3(IRATE+(IR-1)*4)
     &              *DEXP( WR3(IRATE+(IR-1)*4+1) *LNTH
     &              -WR3(IRATE+(IR-1)*4+2) /TH)
     &              *DEXP(  WR3(IRATE+(IR-1)*4+3)*(TE-TH)/(TE*TH))

            KFEQ   = KF(IR)
            CALL KEQUILIBRIUM (NS, NREA, GTH, WR3(ISTOIR),  WR3(ISTOIP),
     &                         WR1(IMI), RTH, LNRTH, KEQ, IR)

             KB(IR) = KFEQ /KEQ

          CASE(21)
            TREFN2 =   6000.D0
!            IF(TH.GT.30000.D0) THEN
!               TH = 30000.D0
!            ENDIF
            IF (TH.LT.TREFN2) THEN
               KF(IR) = WR3(IRATE+(IR-1)*4)
     &              *DEXP(- WR3(IRATE+(IR-1)*4+1)/TH)
     &              *(1.D0 - DEXP(- WR3(IRATE+(IR-1)*4+2)/TH))
     &              *WR1(IUNA)*1.D6

            ELSE
               KF(IR) = WR3(IRATE+(IR-1)*4)
     &              *DEXP(- WR3(IRATE+(IR-1)*4+1)/TH)
     &              *(1.D0 - DEXP(- WR3(IRATE+(IR-1)*4+2)/TE))
     &              *(DEXP(-(WR3(IRATE+(IR-1)*4+1)
     &              - WR3(IRATE+(IR-1)*4+3)*TH) 
     &              *(1.D0/TE - 1.D0/TH)))
     &              *WR1(IUNA)*1.D6
            ENDIF
                KFEQ   = KF(IR)
                CALL KEQUILIBRIUM(NS,NREA,GTH,WR3(ISTOIR), WR3(ISTOIP),
     &                         WR1(IMI), RTH, LNRTH, KEQ, IR)
                KB(IR) = KFEQ /KEQ

          CASE(22)
             TREFNO =   4000.D0
             IF (TH.LT.TREFNO) THEN
                KF(IR) = WR3(IRATE+(IR-1)*4)
     &              *DEXP(- WR3(IRATE+(IR-1)*4+1)/TH)
     &              *(1.D0 - DEXP(- WR3(IRATE+(IR-1)*4+2)/TH))
     &              *WR1(IUNA)*1.D6
            ELSE
                KF(IR) = WR3(IRATE+(IR-1)*4)
     &              *DEXP(- WR3(IRATE+(IR-1)*4+1)/TH)
     &              *(1.D0 - DEXP(- WR3(IRATE+(IR-1)*4+2)/TE))
     &              /(1.D0 - DEXP(- WR3(IRATE+(IR-1)*4+2)/TH))
     &              *(DEXP(-(WR3(IRATE+(IR-1)*4+1)
     &              - WR3(IRATE+(IR-1)*4+3)*TH) 
     &              *(1.D0/TE - 1.D0/TH))) 
     &              *WR1(IUNA)*1.D6
            ENDIF
                KFEQ   = KF(IR)
                CALL KEQUILIBRIUM(NS,NREA,GTH,WR3(ISTOIR), WR3(ISTOIP),
     &                         WR1(IMI), RTH, LNRTH, KEQ, IR)
                KB(IR) = KFEQ /KEQ
          CASE(23)
             TREFO2 =   5000.D0
             IF (TH.LT.TREFO2) THEN
                KF(IR) = WR3(IRATE+(IR-1)*4)
     &              *DEXP(- WR3(IRATE+(IR-1)*4+1)/TH)
     &              *(1.D0 - DEXP(- WR3(IRATE+(IR-1)*4+2)/TH))
     &              *WR1(IUNA)*1.D6
            ELSE
                KF(IR) = WR3(IRATE+(IR-1)*4)
     &              *DEXP(- WR3(IRATE+(IR-1)*4+1)/TH)
     &              *(1.D0 - DEXP(- WR3(IRATE+(IR-1)*4+2)/TE))
     &              *(DEXP(-(WR3(IRATE+(IR-1)*4+1)
     &              - WR3(IRATE+(IR-1)*4+3)*TH) 
     &              *(1.D0/TE - 1.D0/TH))) 
     &              *WR1(IUNA)*1.D6
            ENDIF
                KFEQ   = KF(IR)
                CALL KEQUILIBRIUM(NS,NREA,GTH,WR3(ISTOIR), WR3(ISTOIP),
     &                         WR1(IMI), RTH, LNRTH, KEQ, IR)
                KB(IR) = KFEQ /KEQ
C     -5. Radical reactions (including Zeldovich reactions)
C     -6. Charge exchange
C         => KF(TH), KB(TH)
          CASE(24)
            KF(IR) = WR3(IRATE+(IR-1)*4) 
     &              *DEXP( -WR3(IRATE+(IR-1)*4+1) *LNTH
     &              -WR3(IRATE+(IR-1)*4+2) /TH)*WR1(IUNA)*1.D6
            KFEQ   = KF(IR)
            CALL KEQUILIBRIUM (NS, NREA, GTH, WR3(ISTOIR),  WR3(ISTOIP),
     &                         WR1(IMI), RTH, LNRTH, KEQ, IR)
            KB(IR) = KFEQ /KEQ

          CASE(25)
            KF(IR) = WR3(IRATE+(IR-1)*4) 
     &              *DEXP( -WR3(IRATE+(IR-1)*4+1) *LNTE
     &              -WR3(IRATE+(IR-1)*4+2) /TE)*WR1(IUNA)*1.D6
            KFEQ   = KF(IR)
            CALL KEQUILIBRIUM (NS, NREA, GTE, WR3(ISTOIR),  WR3(ISTOIP),
     &                         WR1(IMI), RTE, LNRTE, KEQ, IR)
            KB(IR) = KFEQ /KEQ
C     -4. Associative ionization
C         => KF(TH), KB(TE)

          CASE(26)

            KF(IR)   = WR3(IRATE+(IR-1)*4) 
     &              *DEXP(-WR3(IRATE+(IR-1)*4+1) *LNTE
     &              -WR3(IRATE+(IR-1)*4+2) /TE)
     &              *WR1(IUNA)*1.D6

            KFEQ = WR3(IRATE+(IR-1)*4) 
     &              *DEXP(-WR3(IRATE+(IR-1)*4+1) *LNTH
     &              -WR3(IRATE+(IR-1)*4+2) /TH)
     &              *WR1(IUNA)*1.D6

            CALL KEQUILIBRIUM (NS, NREA, GTH, WR3(ISTOIR), WR3(ISTOIP), 
     &                         WR1(IMI), RTH, LNRTH, KEQ, IR)
            KB(IR) = KFEQ /KEQ
          CASE DEFAULT
            WRITE(*,*) 'Mechanism not implemented'  
            STOP
        ENDSELECT
      ENDDO

C     Reaction source terms [mole /(m^3 s)]
      DO IR = 1, NREA
        IF (WI(ICATA+IR-1)==0) THEN
          SUM   = 1.D0
        ELSE
C     Contribution of catalytic species (reduced cost implementation)
          SUM = 0.D0
          DO IS = 1, NS  
            SUM = SUM + WR3(IALPHA+(IR-1)*NS+IS-1) *MOLE(IS)     
          ENDDO
        ENDIF
        PRODR = 0.D0
        PRODP = 0.D0
        DO IS = 1, NS
          PRODR = PRODR +WR3(ISTOIR+(IR-1)*NS+IS-1) *LNMOLE(IS)
          PRODP = PRODP +WR3(ISTOIP+(IR-1)*NS+IS-1) *LNMOLE(IS) 
        ENDDO
        PRODR  = DEXP(PRODR)
        PRODP  = DEXP(PRODP)
        XI(IR) = SUM *(KF(IR) *PRODR -KB(IR) *PRODP)
      ENDDO
      DO IS = 1, NS
        OMEGA(IS) = 0.D0
        DO IR = 1, NREA
          OMEGA(IS) = OMEGA(IS) +(WR3(ISTOIP+(IR-1)*NS+IS-1)
     &                            -WR3(ISTOIR+(IR-1)*NS+IS-1)) *XI(IR)
        ENDDO
      ENDDO
      
      END SUBROUTINE ARRHENIUS
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE  KEQUILIBRIUM (NS, NREA, G, STOIR, STOIP, MI, RT, LNRT,
     &                          KEQ, IR)
C-----------------------------------------------------------------------
C     This subroutine computes the equilibrium constants. Gibbs 's free
C     energy is given per unit mole.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NS, NREA, IR
      DOUBLE PRECISION G(1:NS), MI(1:NS), RT, LNRT, KEQ, 
     &                 STOIR(1:NS*NREA), STOIP(1:NS*NREA)
C-----------------------------------------------------------------------
      INTEGER IS
      DOUBLE PRECISION SUM, RTLNRT
C-----------------------------------------------------------------------
      RTLNRT = RT *LNRT

      SUM = 0
      DO IS = 1, NS
        SUM = SUM +(STOIP((IR-1)*NS+IS) -STOIR((IR-1)*NS+IS)) 
     &        *(G(IS) +RTLNRT)
      ENDDO
      KEQ = DEXP(-SUM /RT)

      END SUBROUTINE  KEQUILIBRIUM
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE READRATES (PATH, REACTION, WI, LWI, WC, LWC, WR1, LWR1,
     &                      WR3, LWR3)
C-----------------------------------------------------------------------
C     This subroutine reads the reactions and corresponding rates in 
C     data files. The routine detects whether the same reaction is
C     specified several times. The catalytic reactions can be grouped
C     together by using a generic M symbol in order to reduce the CPU
C     cost. Finally, stoichiometry of the reactions is checked.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      INTEGER LWI, WI(1:LWI), LWR1, LWR3, LWC
      DOUBLE PRECISION WR1(1:LWR1), WR3(1:LWR3)
      CHARACTER(10) REACTION
      CHARACTER(100) PATH
      CHARACTER WC(1:LWC)
C-----------------------------------------------------------------------
      INTEGER LPATH, LREACTION, ILEN, I, IC, IM, ITYPE, IARRHENIUS, IR, 
     &        IS, J, LCHAR, ISPECIES, NTHERMAL, IS1, IS2, IS3, IS4, IS5,
     &        ISDIS, JR, JJ, N1, N2, IA1, IA2, NALPHA(1:NREA), 
     &        ALPHA(1:NREA,1:NS)
      DOUBLE PRECISION COEFFICIENT, ARRHENIUS(1:4), 
     &                 REACTANT(1:NS), PRODUCT(1:NS),
     &                 ELEMENTR(1:NC), ELEMENTP(1:NC), SUM,
     &                 VEC1(0:NS), VEC2(0:NS), VEC3(0:NS), VEC4(0:NS),
     &                 NVEC1, NVEC2, DIF1, DIF2, DIF
      LOGICAL PROP
      CHARACTER(4) COM
      CHARACTER(10) WORD, SPECIES
      CHARACTER(100) FULLCOM, LINE
C-----------------------------------------------------------------------
      LPATH     = LCHAR (PATH)
      LREACTION = LCHAR (REACTION)
      DO I = 1, LWR3
        WR3(I) = 0.D0
      ENDDO
      IR = 1

      OPEN(UNIT=INOUT3,FILE=PATH(1:LPATH)//'/data/chemistry/gasreact/'
     &     //REACTION(1:LREACTION),STATUS='OLD')
      CALL BLANK (FULLCOM)
      READ(INOUT3,'(A)') FULLCOM
      CALL BLANK (LINE)
      CALL RMBLANK (FULLCOM, LINE, ILEN)
      COM = LINE(1:4)
      DO WHILE (COM(1:4) /= 'STOP')
        DO I = 1, NS
          REACTANT(I) = 0.D0
          PRODUCT(I)  = 0.D0
        ENDDO
C     Start with a word
        CALL BLANK (WORD)
        IC = 0
C     Distinction between reactant, product, and rate
        ITYPE = 0
C     First species /= third body
        IM = 0
        DO I = 1, ILEN
          IC = IC +1
C     Coefficient + species (reactant or product)
          IF (((LINE(I:I)=='+').AND.(ITYPE/=2)).OR.(LINE(I:I)=='=')
     &        .OR.(LINE(I:I)==':')) THEN
            CALL CUTWORD (WORD, COEFFICIENT, SPECIES)
            CALL TESTSPECIES (SPECIES, WI, LWI, WC, LWC, ISPECIES)
            IF (ISPECIES==-1) THEN
              IM = IM +1
            ELSEIF (ISPECIES==0) THEN
              WRITE(*,*) 'Species ', SPECIES,
     &                   ' does not exist in reaction ', LINE
            ELSE
              IF (ITYPE==0) THEN
                REACTANT(ISPECIES) = REACTANT(ISPECIES) +COEFFICIENT
              ELSE
                PRODUCT(ISPECIES)  = PRODUCT(ISPECIES)  +COEFFICIENT
              ENDIF
            ENDIF
            CALL BLANK (WORD) 
            IC = 0
C     Next word = products
            IF (LINE(I:I)=='=') THEN 
              ITYPE = 1
            ENDIF
C     Next word = rates
            IF (LINE(I:I)==':') THEN 
              ITYPE = 2
              IARRHENIUS = 0
            ENDIF
C     Arrhenius rates
          ELSEIF ((LINE(I:I)=='/').OR.(I==ILEN)) THEN
            IF (I==ILEN) THEN
              WORD(IC:IC) = LINE(I:I)
            ENDIF
            IARRHENIUS = IARRHENIUS +1
            IF (IARRHENIUS>5) THEN
              WRITE(*,*) 'More than 4 rates + 1 thermal nonequilibrium', 
     &                   ' mechanism in reaction ', LINE
            ELSEIF (IARRHENIUS==5) THEN
              CALL ITRANSLATE (WORD, NTHERMAL)
            ELSE
              CALL DTRANSLATE (WORD, ARRHENIUS(IARRHENIUS))
            ENDIF
            CALL BLANK (WORD)
            IC = 0
          ENDIF
          IF (IC/=0) THEN
            WORD(IC:IC) = LINE(I:I)
          ENDIF
        ENDDO
        CALL BLANK (FULLCOM)
        READ(INOUT3,'(A)') FULLCOM
        CALL BLANK (LINE)
        CALL RMBLANK (FULLCOM, LINE, ILEN)
        COM = FULLCOM(1:4)
        IF (IM>0) THEN
C     M = specificied species
          IF ((LINE(1:1)=='M').AND.(LINE(2:2)=='=')) THEN
            DO I = 3, ILEN
              IC = IC +1
              IF ((LINE(I:I)==',').OR.(I==ILEN)) THEN
                IF (I==ILEN) THEN
                  WORD(IC:IC) = LINE(I:I)
                ENDIF     
                CALL TESTSPECIES (WORD, WI, LWI, WC, LWC, ISPECIES)
                DO IS = 1, NS
                  WR3(ISTOIR+(IR-1)*NS+IS-1) = REACTANT(IS)
                  WR3(ISTOIP+(IR-1)*NS+IS-1) = PRODUCT(IS)
                ENDDO
                WR3(IALPHA+(IR-1)*NS+ISPECIES-1) = 1.D0
                DO J = 1, 4
                  WR3(IRATE+(IR-1)*4+J-1) = ARRHENIUS(J)
                ENDDO
                WI(ICATA+IR-1) = 1
                WI(ITNEQ+IR-1) = NTHERMAL
                CALL BLANK (WORD)   
                IC = 0
              ENDIF
              IF (IC/=0) THEN
                WORD(IC:IC) = LINE(I:I)
              ENDIF
            ENDDO
            IR = IR +1
            CALL BLANK (FULLCOM)
            CALL BLANK (LINE)
            READ(INOUT3,'(A)') FULLCOM
            CALL RMBLANK (FULLCOM, LINE, ILEN)
            COM = FULLCOM(1:4)
C     M = all species
          ELSE
            DO IS = 1, NS
              WR3(ISTOIR+(IR-1)*NS+IS-1) = REACTANT(IS)
              WR3(ISTOIP+(IR-1)*NS+IS-1) = PRODUCT(IS)
              WR3(IALPHA+(IR-1)*NS+IS-1) = 1.D0
            ENDDO
            DO J = 1, 4
              WR3(IRATE+(IR-1)*4+J-1) = ARRHENIUS(J)
            ENDDO         
            WI(ICATA+IR-1) = 1
            WI(ITNEQ+IR-1) = NTHERMAL
            IR = IR +1
          ENDIF
          IM = 0
C     Noncatalytic reaction (or catalytic reaction without any special 
C     treatment to reduce CPU cost)
        ELSE
          DO IS = 1, NS
            WR3(ISTOIR+(IR-1)*NS+IS-1) = REACTANT(IS)
            WR3(ISTOIP+(IR-1)*NS+IS-1) = PRODUCT(IS)
          ENDDO
          DO J = 1, 4
            WR3(IRATE+(IR-1)*4+J-1) = ARRHENIUS(J)
          ENDDO
          WI(ICATA+IR-1) = 0
          WI(ITNEQ+IR-1) = NTHERMAL
          IR = IR +1
        ENDIF
        IF (IARRHENIUS<5) THEN
          WRITE(*,*) 'Only', IARRHENIUS, 
     &               ' coefficients in reaction ', IR
          PAUSE
        ENDIF
      ENDDO
      CLOSE(INOUT3)

C     Check stoichiometry of reactions
      DO IR = 1, NREA
        DO IC = 1, NC
          ELEMENTR(IC) = 0.D0
          ELEMENTP(IC) = 0.D0
        ENDDO
        DO IC = 1, NC
          IS = IC
          ELEMENTR(IC) = ELEMENTR(IC) +WR3(ISTOIR+(IR-1)*NS+IS-1)
          ELEMENTP(IC) = ELEMENTP(IC) +WR3(ISTOIP+(IR-1)*NS+IS-1)
        ENDDO
        DO IS = NC+1, NS
          DO IC = 1, NC
            ELEMENTR(IC) = ELEMENTR(IC) -WR3(ISTOIR+(IR-1)*NS+IS-1)
     &                     *WI(INUIK+(IS-1)*NS+IC-1)
     &                     /WI(INUIK+(IS-1)*NS+IS-1)
            ELEMENTP(IC) = ELEMENTP(IC) -WR3(ISTOIP+(IR-1)*NS+IS-1)
     &                     *WI(INUIK+(IS-1)*NS+IC-1)
     &                     /WI(INUIK+(IS-1)*NS+IS-1)
          ENDDO
        ENDDO
        DO IC = 1, NC
          IF (ELEMENTR(IC)/=ELEMENTP(IC)) THEN
            WRITE(*,*) 'Error in stoichiometry of element ', IC, 
     &                 ' of reaction ', IR, ' over',
     &                 NREA, ' reactions.'
            WRITE(*,*) 'Reactants-Products'
            IS = IC
            WRITE(*,*) IS, WR3(ISTOIR+(IR-1)*NS+IS-1), 
     &                 WR3(ISTOIP+(IR-1)*NS+IS-1)
            PAUSE
          ENDIF
        ENDDO
      ENDDO
      WRITE(*,*) '    Chemical reactions satisfy stoichiometry.'
      WRITE(*,*) 

C     Check whether the reactions are not written twice.
C      PROP = .FALSE.
      PROP = .TRUE.
      IF (PROP) THEN
        DO IR = 1, NREA
          NALPHA(IR) = 0
          DO IS = 1, NS
              ALPHA(IR,IS) = 0
          ENDDO
          DO IS = 1, NS
            IF (WR3(IALPHA+(IR-1)*NS+IS-1)>1.D-6) THEN
              NALPHA(IR) = NALPHA(IR) +1
              ALPHA(IR,NALPHA(IR)) = IS
            ENDIF 
          ENDDO
          IF (NALPHA(IR) == 0) THEN
            NALPHA(IR) = 1
          ENDIF
        ENDDO
        DO IR = 1, NREA-1
          DO JR = IR+1, NREA
            NVEC1 = 0.D0; NVEC2 = 0.D0
            DO IS = 1, NS
              VEC1(IS) = WR3(ISTOIR+(IR-1)*NS+IS-1)
     &                  -WR3(ISTOIP+(IR-1)*NS+IS-1)         
              VEC2(IS) = WR3(ISTOIR+(JR-1)*NS+IS-1)
     &                  -WR3(ISTOIP+(JR-1)*NS+IS-1)
              NVEC1    = NVEC1 + VEC1(IS)*VEC1(IS)
              NVEC2    = NVEC2 + VEC2(IS)*VEC2(IS)
            ENDDO
            NVEC1 = DSQRT(NVEC1)
            NVEC2 = DSQRT(NVEC2)
            DIF1 = 0.D0; DIF2 = 0.D0
            DO IS = 1, NS
              DIF  = VEC1(IS) /NVEC1 -VEC2(IS) /NVEC2
              DIF1 = DIF1 +DIF*DIF
              DIF  = VEC1(IS) /NVEC1 +VEC2(IS) /NVEC2
              DIF2 = DIF2 +DIF*DIF
            ENDDO 
            DIF = DIF1*DIF2
C    Additional test for catalytic reactions
            IF ((DIF<1.D-10).AND.(IR/=JR)) THEN
              DO I = 1, NALPHA(IR)
              DO J = 1, NALPHA(JR)
                VEC1(0) = 0.D0; VEC2(0) = 0.D0
                VEC3(0) = 0.D0; VEC4(0) = 0.D0
                DO IS = 1, NS
                  VEC1(IS) = WR3(ISTOIR+(IR-1)*NS+IS-1)
                  VEC2(IS) = WR3(ISTOIP+(IR-1)*NS+IS-1)
                  VEC3(IS) = WR3(ISTOIR+(JR-1)*NS+IS-1)
                  VEC4(IS) = WR3(ISTOIP+(JR-1)*NS+IS-1)
                ENDDO
                VEC1(ALPHA(IR,I)) = VEC1(ALPHA(IR,I)) +1
                VEC2(ALPHA(IR,I)) = VEC2(ALPHA(IR,I)) +1
                VEC3(ALPHA(JR,J)) = VEC3(ALPHA(JR,J)) +1
                VEC4(ALPHA(JR,J)) = VEC4(ALPHA(JR,J)) +1
                NVEC1 = 0.D0; NVEC2 = 0.D0
                DO IS = 1, NS
                  NVEC1    = NVEC1 +VEC1(IS)*VEC1(IS) +VEC2(IS)*VEC2(IS)
                  NVEC2    = NVEC2 +VEC3(IS)*VEC3(IS) +VEC4(IS)*VEC4(IS)
                ENDDO
                NVEC1 = DSQRT(NVEC1)
                NVEC2 = DSQRT(NVEC2)
                DIF1 = 0.D0; DIF2 = 0.D0
                DO IS = 1, NS
                  DIF  =  (VEC1(IS) +VEC2(IS))/NVEC1
     &                   -(VEC3(IS) +VEC4(IS))/NVEC2
                  DIF1 = DIF1 +DIF*DIF
                  DIF  =  (VEC1(IS) +VEC2(IS))/NVEC1
     &                   +(VEC3(IS) +VEC4(IS))/NVEC2
                  DIF2 = DIF2 +DIF*DIF
                ENDDO
                DIF = DIF1*DIF2
                IF (DIF<1.D-10) THEN
                  WRITE(*,*) 'Reactions',IR,' and', JR,' identical'
                  PAUSE
                  PROP = .FALSE.
                ENDIF
              ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        IF (PROP) THEN
          WRITE(*,*) '    Chemical reactions are all different.'
        ELSE
          WRITE(*,*) '    Some chemical reactions are identical.'
          PAUSE
        ENDIF
      ELSE
        WRITE(*,*) '    Chemical reactions difference not checked.'
      ENDIF
      WRITE(*,*) 

C     Selection of the molecule dissociated by atomic or molecular 
C     impact in a reaction of dissociation. This molecule is
C     denoted here by IDIS and stored in WI(IDISSOCIATION+IR-1).
      DO IR = 1, NREA
        IF (WI(ITNEQ+IR-1)==1) THEN
          ISDIS = 0; IC = 0; IS1 = 0; IS2 = 0
          DO IS = 1, NS
            IF (INT(WR3(ISTOIR+(IR-1)*NS+IS-1))==1) THEN
              IF (IC == 0) THEN 
                IS1 = IS
              ELSE
                IS2 = IS
              ENDIF
              IC = IC +1
            ELSEIF (INT(WR3(ISTOIR+(IR-1)*NS+IS-1))==2) THEN
C     Reaction: 2 AB = A +B +AB
              ISDIS = IS             
            ENDIF 
          ENDDO
          IF (ISDIS == 0) THEN
            IF (IC == 1) THEN
C     Catalytic reaction with special treatment to reduce the CPU cost
              ISDIS = IS1
            ELSE
              IC = 0
              DO IS = 1, NS
                IF (INT(WR3(ISTOIP+(IR-1)*NS+IS-1))==3) THEN
C     Reaction: AA + A = 3A
                  IF (IS1==IS) THEN
                    ISDIS = IS2
                  ELSE
                    ISDIS = IS1
                  ENDIF
                ELSEIF (INT(WR3(ISTOIP+(IR-1)*NS+IS-1))==2) THEN
                  IF (IC == 0) THEN
                    IS3 = IS
                  ELSEIF (IC == 1) THEN
                    IS4 = IS
                  ENDIF
                  IC = IC +1
                ELSEIF (INT(WR3(ISTOIP+(IR-1)*NS+IS-1))==1) THEN
                  IF (IC == 0) THEN
                    IS3 = IS
                  ELSEIF (IC == 1) THEN
                    IS4 = IS
                  ELSEIF (IC == 2) THEN
                    IS5 = IS
                  ENDIF
                  IC = IC +1
                ENDIF
              ENDDO
              IF (ISDIS == 0) THEN
                IF (IC == 2) THEN
C     Reactions:
C                AA +X = 2A  +X (X/=A)
C                AB +A =  B +2A (A/=B) 
C                AB +B =  A +2B (A/=B) 
C     ISDIS is the reactant not retrieved as a product
                  IF ((IS1/=IS3).AND.(IS1/=IS4)) THEN
                    ISDIS = IS1
                  ELSE
                    ISDIS = IS2  
                  ENDIF       
                ELSEIF (IC == 3) THEN
C     Reaction: AB +X = A +B +X (A/=B, X/=A, X/=B)
C     ISDIS is the reactant not retrieved as a product
                  IF ((IS1/=IS3).AND.(IS1/=IS4).AND.(IS1/=IS5)) THEN
                    ISDIS = IS1
                  ELSE
                    ISDIS = IS2  
                  ENDIF               
                ENDIF
              ENDIF
            ENDIF
          ENDIF
          IF (ISDIS/=0) THEN
            WI(IDISSOCIATION+IR-1)= ISDIS
          ELSE
            WRITE(*,*) 'Error in the dissociation reaction', IR
            PAUSE
          ENDIF         
        ENDIF
      ENDDO

C     Unit conversion for the rates
      DO IR = 1, NREA
C     cm^3, mol, and s => m^3, mol, and s 
        SUM = 0
        DO IS = 1, NS 
          SUM = SUM +WR3(ISTOIR+(IR-1)*NS+IS-1)
        ENDDO
        IF (WI(ICATA+IR-1)==1) THEN
          SUM = SUM +1
        ENDIF
C        J = 1     
C        WR3(IRATE+(IR-1)*4+J-1) = WR3(IRATE+(IR-1)*4+J-1)*WR1(IUNA)
C     &                            *(1.D+6)**(SUM -1.D0) 

        J = 1     
        WR3(IRATE+(IR-1)*4+J-1) = WR3(IRATE+(IR-1)*4+J-1) 
     &                            *(1.D-6)**(SUM -1.D0)
C     kcal/mole      => K
        J = 3         
C        WR3(IRATE+(IR-1)*4+J-1) = WR3(IRATE+(IR-1)*4+J-1) 
C     &                            *4186.8 /WR1(IUR)
      ENDDO  

      END SUBROUTINE READRATES
C----------------------------------------------------------------------- 
C-----------------------------------------------------------------------
      SUBROUTINE CUTWORD (WORD, COEFFICIENT, SPECIES)
C-----------------------------------------------------------------------
C     This subroutine cuts a word composed of one coefficient and one 
C     species. The default value of the coefficient is 1.D0.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(*) WORD, SPECIES
      DOUBLE PRECISION COEFFICIENT
C-----------------------------------------------------------------------
      INTEGER I, J, LWORD, ICO
      INTEGER LCHAR
C-----------------------------------------------------------------------
      DO I = 1, LEN(SPECIES)
        SPECIES(I:I) = " "
      ENDDO
      LWORD = LCHAR(WORD)
      COEFFICIENT = 0.D0
      ICO = 0

      J = 0
      DO I = 1, LWORD
        IF ((J==0).AND.(ICHAR(WORD(I:I))>=48)
     &      .AND.(ICHAR(WORD(I:I))<=57)) THEN
          COEFFICIENT = 10.D0 *COEFFICIENT +ICHAR(WORD(I:I)) -48.D0
          ICO = 1
        ELSE
          IF ((LWORD-I) > LEN(SPECIES)) THEN
            WRITE(*,*) 'Problem in CUTWORD'
          ENDIF
          J = J+1
          SPECIES(J:J) = WORD(I:I)
        ENDIF
      ENDDO
      IF (ICO == 0) THEN
        COEFFICIENT = 1.D0
      ENDIF

      END SUBROUTINE CUTWORD
C----------------------------------------------------------------------- 
C-----------------------------------------------------------------------
      SUBROUTINE TESTSPECIES (SPECIES, WI, LWI, WC, LWC, ISPECIES)
C-----------------------------------------------------------------------
C     This subroutine finds the index (=ISPECIES) of the species 
C     (=SPECIES) in the mixture. 
C     ISPECIES = -1 if SPECIES = M (third body)
C              =  0 if SPECIES does not belong to the mixture.
C              =  I, I \in {1,NS} otherwise.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER LWI, LWC, WI(1:LWI), ISPECIES
      CHARACTER(10) SPECIES
C     The following declaration allows for WC to be out of range      
      CHARACTER WC(1:LWC)
C-----------------------------------------------------------------------
      INTEGER I, J, LTOT1, LTOT2, IC1, LCHAR
      LOGICAL PROP
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      LTOT1 = LCHAR(SPECIES)
      ISPECIES = 0
      IC1 = 0
      DO I = 1, NS
        LTOT2 = WI(ILNAMEI+I-1)
        IF (LTOT1==LTOT2) THEN
          J = 1
          PROP = .FALSE.
          IF ((INAMEI+IC1+J-1)<=LWC) THEN
            IF (SPECIES(J:J) == WC(INAMEI+IC1+J-1)) THEN
              PROP = .TRUE.
            ENDIF
          ENDIF
          DOWHILE (PROP.AND.(J<=LTOT2))  
            J = J+1       
            PROP = .FALSE.
            IF ((INAMEI+IC1+J-1)<=LWC) THEN
              IF (SPECIES(J:J) == WC(INAMEI+IC1+J-1)) THEN
                PROP = .TRUE.
              ENDIF
            ENDIF
          ENDDO 
          IF (J > LTOT2) THEN 
            ISPECIES = I
          ENDIF
        ENDIF
        IC1 = IC1 + LTOT2
      ENDDO

      IF ((ISPECIES == 0).AND.(LTOT1==1).AND.(SPECIES(1:1)=='M')) THEN
        ISPECIES = -1
      ENDIF

      END SUBROUTINE TESTSPECIES
C----------------------------------------------------------------------- 
C-----------------------------------------------------------------------
      SUBROUTINE DTRANSLATE (WORD, NUMBER)
C-----------------------------------------------------------------------
C     This subroutine translates a number written in CHARACTER to 
C     DOUBLE PRECISION.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(*) WORD
      DOUBLE PRECISION NUMBER
C-----------------------------------------------------------------------
      DOUBLE PRECISION NUMBER1, NUMBER2, DECI1, DECI2
      INTEGER I, LCHAR, PM1, PM2, POWER
C-----------------------------------------------------------------------
      NUMBER1 = 0.D0  ; NUMBER2 = 0.D0
      PM1     = 1     ; PM2     = 1
      DECI1    = 10.D0; DECI2   = 1.D0 
      POWER   = 0
      DO I = 1, LCHAR(WORD)
        IF (WORD(I:I)=='.') THEN
          DECI1  = 1.D0
          DECI2  = 0.1D0
        ELSEIF ((WORD(I:I)=='E').OR.(WORD(I:I)=='e')) THEN
          POWER = 1
          DECI1    = 10.D0; DECI2   = 1.D0 
        ELSEIF (WORD(I:I)=='+') THEN
          CONTINUE
        ELSEIF (WORD(I:I)=='-') THEN
          IF (POWER == 0) THEN
            PM1 = -1
          ELSE
            PM2 = -1
          ENDIF
        ELSEIF ((ICHAR(WORD(I:I))>=48)
     &          .AND.(ICHAR(WORD(I:I))<=57)) THEN
          IF (POWER == 0) THEN
            NUMBER1 = DECI1 *NUMBER1 + DECI2 *(ICHAR(WORD(I:I)) -48.D0)
          ELSE
            NUMBER2 = DECI1 *NUMBER2 + DECI2 *(ICHAR(WORD(I:I)) -48.D0)
          ENDIF
          IF (ABS(DECI2-1.D0)<1.D-3) THEN
            DECI2 = 1.D0
          ELSE
            DECI2 = DECI2*0.1D0
          ENDIF
        ELSEIF (WORD(I:I)==' ') THEN
        ELSE
          WRITE(*,*) 'Character not recognized in DTRANSLATE'
          PAUSE
        ENDIF
      ENDDO

      IF (POWER==0) THEN
        NUMBER = PM1 *NUMBER1
      ELSE
        NUMBER = PM1 *NUMBER1 *EXP(PM2 *NUMBER2 *LOG(10.D0))
      ENDIF

      END SUBROUTINE DTRANSLATE
C----------------------------------------------------------------------- 
C-----------------------------------------------------------------------
      SUBROUTINE ITRANSLATE (WORD, NUMBER)
C-----------------------------------------------------------------------
C     This subroutine translates a number written in CHARACTER to 
C     INTEGER.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(*) WORD
      INTEGER NUMBER, LCHAR, I
C-----------------------------------------------------------------------
      NUMBER = 0
      DO I = 1, LCHAR(WORD)
            NUMBER = 10 *NUMBER + (ICHAR(WORD(I:I)) -48.D0)
      ENDDO

      END SUBROUTINE ITRANSLATE
C----------------------------------------------------------------------- 