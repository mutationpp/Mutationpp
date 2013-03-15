C-----------------------------------------------------------------------
      SUBROUTINE ARRHENIUS (WR1, LWR1, WR3, LWR3, WI, LWI, Y, YTOL, P, 
     &                     TEFS, RHO, KF, ESCAPE, OMEGA, OMEGACV,
     &                     OMEGAI)
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
      INCLUDE '../general/memorycr.cmn'
C-----------------------------------------------------------------------
      INTEGER  LWR1, LWR3, LWI, WI(1:LWI), J, I, FLAGION, FLAGI, FLAGJ,
     &         LIAUX(1:20),LJAUX(1:20)
      DOUBLE PRECISION WR1(1:LWR1), WR3(1:LWR3), Y(1:NS), YTOL, P, 
     &              TEF(1:NVIB+2), RHO, OMEGA(1:NS), TA, TEFS(1:NVIB+2)
C-----------------------------------------------------------------------
      DOUBLE PRECISION XI(1:NREA), KF(1:NREA), KB(1:NREA), KEQ, KFEQ,
     &                 GTH(1:NS), GTE(1:NS), YT(1:NS), THV(1:NV+1), 
     &                 TEV(1:NV+1), PRODP, PRODR, MOLE(1:NS), 
     &                 LNMOLE(1:NS), TH, TV(1:NV), TE, SUM, LNR,
     &                 RTH, LNTH, LNRTH, RTE, LNTE,  LNRTE, 
     &                 STHTV(1:NS), LNSTHTV(1:NS), TVIB(1:NS), 
     &                 SE(1:NVIB+2), LNSE(1:NVIB+2),
     &                 UTE, TRATIO, Ia, KeV, TREFO2, TREFN2, TREFNO,
     &                 RTA, LNTA, LNRTA, AJI, ESCAPE
      INTEGER IS, IR, IV, ISDIS, NVIBMODE, IC, JR, S, IS1, IS2, IVS,
     &                 GAUX(1:20),GAUX2(1:20), NA, NB, IVT
      DOUBLE PRECISION  THCLIP, LNTHCLIP, RTHCLIP, LNRTHCLIP,
     &                  GTECLIP(1:NS), TECLIP, LNTECLIP, RTECLIP, 
     &                  LNRTECLIP, GTHCLIP(1:NS),xCNBTg ,xCNXTg,xCNATg,
     &                  PSI4, PSI3, PSI5, TVI, EAUX(1:20),EAUX2(1:20)
      DOUBLE PRECISION H1(1:NS), H2(1:NS), H3(1:NS), H4(1:NS), H5(1:NS),
     &                 H6(1:NS), GVAI, GVAI1, GVAI2, OMEGACV(1:NVIB),
     &                 OMEGAI, EI, EJ, EJJ, KIJ, NION, OION, fy, psi,
     &                 EIN, IH, EIKELVIN, ENU, ENL, 
     &                 TVN2, ALPHA, PSITM, U, TF, PAR, PSIK(1:6),
     &                 PSIHAM, PSIHAS, PSIHANS, PHI, PSIMACH, CV
      CHARACTER   NUMBERC
      EXTERNAL KNAB, TREANOR, LOSEV, HAMMERLING, HASSAN, HANSEN,
     &         MACHERET

C-----------------------------------------------------------------------
C     Fix to avoid problems with the precision for low temperatures
C     This may kill the convergency of Newton iterator

      PAR = 3.D0
      ALPHA = 0.7D0
      DO I=1,NVIB+2
         TEF(I) = TEFS(I)
         IF(TEF(I).LT.250.D0) THEN
            TEF(I) = 250.D0
         ENDIF
!          IF(TEF(I).GT.30000.D0) THEN
!             TEF(I) = 30000.D0
!          ENDIF
      ENDDO
C-----------------------------------------------------------------------

C     Translational temperature of heavy particles
      TH = TEF(1)
      THCLIP = TEF(1)
C     Translational temperature of electrons
      TE = TEF(NVIB+2)
      TECLIP = TEF(NVIB+2)
C      IF (THCLIP.GT.15000.D0) THCLIP = 15000.D0
C      IF (THCLIP.GT.30000.D0) THCLIP = 30000.D0
C      IF (TECLIP.GT.30000.D0) TECLIP = 30000.D0
C     Vibrational temperatures
      IC = 0
      DO IS = 1, NS
        NVIBMODE = WI(IVIBI+IS-1)
        DO IV = 1, NVIBMODE
          IC = IC +1
          TV(IC)   = TEF(WI(IVIBTEMPI+IS-1))
        ENDDO
      ENDDO

C     Translational temperature vectors of size (1:NV)
      DO IV = 1, NV
        THV(IV) = TH
        TEV(IV) = TE
      ENDDO
      TA = DSQRT(TH*TE)
C     Logarithms
      LNR    = DLOG(WR1(IUR))
      RTH    = WR1(IUR) *TH
      LNTH   = DLOG(TH)
      LNRTH  = LNR +LNTH
      RTHCLIP= WR1(IUR) *THCLIP
      LNTHCLIP= DLOG(THCLIP)
      LNRTHCLIP  = LNR +LNTHCLIP
      RTE    = WR1(IUR) *TE
      LNTE   = DLOG(TE)
      LNRTE  = LNR +LNTE
      RTECLIP= WR1(IUR) *TECLIP
      LNTECLIP= DLOG(TECLIP)
      LNRTECLIP= LNR +LNTECLIP
      RTA    = WR1(IUR) *TA
      LNTA   = DLOG(TA)
      LNRTA  = LNR +LNTA
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
C     Gibbs's free energy (per unit mole) at p = 1.D0 Pa
      CALL GIBBS (WR1, LWR1, WI, LWI, THCLIP, THCLIP, TH, THV, 
     &            1.D0, GTHCLIP)
      CALL GIBBS (WR1, LWR1, WI, LWI, TECLIP, TECLIP, TECLIP, TEV,
     &            1.D0, GTECLIP)
      
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
C         MARCO TAKE A LOOK PARK BOOK (pag 272)
C         => KF(TA), KB(TE)
          CASE(4)
C            KF(IR) = WR3(IRATE+(IR-1)*4) 
C     &              *DEXP( WR3(IRATE+(IR-1)*4+1) *LNTA
C     &              -WR3(IRATE+(IR-1)*4+2) /TA)
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
          CASE(5)
            KF(IR) = WR3(IRATE+(IR-1)*4) 
     &              *DEXP( WR3(IRATE+(IR-1)*4+1) *LNTH
     &              -WR3(IRATE+(IR-1)*4+2) /TH)
            KFEQ   = KF(IR)
            CALL KEQUILIBRIUM (NS, NREA, GTH, WR3(ISTOIR),  WR3(ISTOIP),
     &                         WR1(IMI), RTH, LNRTH, KEQ, IR)
            KB(IR) = KFEQ /KEQ

C     -5. Radical reactions (including Zeldovich reactions)
C      CLIPPED !!!!
          CASE(6)
            KF(IR) = WR3(IRATE+(IR-1)*4) 
     &              *DEXP( WR3(IRATE+(IR-1)*4+1) *LNTHCLIP
     &              -WR3(IRATE+(IR-1)*4+2) /THCLIP)
            KFEQ   = KF(IR)
            CALL KEQUILIBRIUM (NS, NREA, GTH, WR3(ISTOIR),  WR3(ISTOIP),
     &                         WR1(IMI), RTHCLIP, LNRTHCLIP, KEQ, IR)
            KB(IR) = KFEQ /KEQ

C         Electron Impact excitation reactions for Nitrogen and Oxygen
          CASE(8)  ! Allowed Transitions (l_i /= l_j)

            UTE= DSQRT((8.D0*WR1(IUKB)*TE)/(WR1(IUPI)
     &                             *WR1(IMI)/WR1(IUNA)))

            TRATIO = (WR3(IRATE+(IR-1)*4+2) - WR3(IRATE+(IR-1)*4+1))*
     &              1.43876866D0/TE

            Ia = 0.63255D0 * (TRATIO**(-1.6454D0)) * DEXP(-TRATIO) 

            KF(IR) = (UTE*4.D0*WR1(IUPI)*((0.529D-10)**2.D0)
     &               *0.05D0*(157820.40447D0/TE)**2.D0 *Ia)*WR1(IUNA)

            KFEQ   = KF(IR)

            CALL KEQUILIBRIUM (NS, NREA, GTE, WR3(ISTOIR), WR3(ISTOIP), 
     &                         WR1(IMI), RTE, LNRTE, KEQ, IR)
            KB(IR) = KFEQ /KEQ

          CASE(9)  ! Optically Forbidden Transitions (l_i = l_j)
            UTE= DSQRT((8.D0*WR1(IUKB)*TE)/(WR1(IUPI)
     &                             *WR1(IMI)/WR1(IUNA)))

            TRATIO = (WR3(IRATE+(IR-1)*4+2) - WR3(IRATE+(IR-1)*4+1))*
     &              1.43876866D0/TE

            Ia = 0.23933D0 * (TRATIO**(-1.4933D0)) * DEXP(-TRATIO) 

            KF(IR) = (UTE*4.D0*WR1(IUPI)*((0.529D-10)**2.D0)
     &               *0.05D0*(TRATIO)**2.D0 * Ia)*WR1(IUNA)

            KFEQ   = KF(IR)
            CALL KEQUILIBRIUM (NS, NREA, GTE, WR3(ISTOIR), WR3(ISTOIP), 
     &                         WR1(IMI), RTE, LNRTE, KEQ, IR)

            KB(IR) = KFEQ /KEQ

C         Electron Impact Ionisation reactions for Nitrogen and Oxygen


          CASE(10)  !Electron Impact Ionisation

            UTE= DSQRT((8.D0*WR1(IUKB)*TE)/(WR1(IUPI)
     &                             *WR1(IMI)/WR1(IUNA)))

          TRATIO = (WR3(IRATE+(IR-1)*4+2) - WR3(IRATE+(IR-1)*4+1))*
     &              1.43876866D0/TE

            Ia = 0.63255D0 * (TRATIO**(-1.6454D0)) * DEXP(-TRATIO) 

            KF(IR) = (UTE*4.D0*WR1(IUPI)*((0.529D-10)**2.D0)
     &              *(157820.40447D0/TE)**2.D0 * Ia)*WR1(IUNA)


            KFEQ   = KF(IR)

            CALL KEQUILIBRIUM (NS, NREA, GTE, WR3(ISTOIR), WR3(ISTOIP), 
     &                         WR1(IMI), RTE, LNRTE, KEQ, IR)
            KB(IR) = KFEQ /KEQ

          CASE(11)  ! Spontaneus Emission !! sec^-1 
!           IF((WR3(IRATE+(IR-1)*4+2)-WR3(IRATE+(IR-1)*4+3))
!     &            .GT.50000.0) THEN
!            KF(IR) = ESCAPE*1.D-40
!              ELSE           
            KF(IR) = (ESCAPE)*WR3(IRATE+(IR-1)*4+1) 
!           ENDIF

            KB(IR) = 1.D-40

          CASE(12)  ! Spontaneus Molecules Emission !! sec^-1 
            I = IDINT(WR3(IRATE+(IR-1)*4+2))
            J = IDINT(WR3(IRATE+(IR-1)*4+1))
            CALL ProbaEmiSpontMolec(I,J,TE,AJI)
            KF(IR) = (ESCAPE)*AJI
            KB(IR) = 1.D-40


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
            TREFN2 = 6000.D0
            IF (TH.LE.TREFN2) THEN
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
C     -3. Ionization by electron impact / recombination
          CASE(25)
            KF(IR) = WR3(IRATE+(IR-1)*4) 
     &              *DEXP(-WR3(IRATE+(IR-1)*4+1) *LNTECLIP
     &              -WR3(IRATE+(IR-1)*4+2) /TECLIP)*WR1(IUNA)*1.D6
            KFEQ   = KF(IR)
            CALL KEQUILIBRIUM (NS, NREA, GTECLIP, WR3(ISTOIR),  
     &               WR3(ISTOIP), WR1(IMI), RTECLIP, LNRTECLIP, KEQ, IR)
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
C     - Heavy particle impact (three body reaction)
          CASE(27)
            KF(IR) = WR3(IRATE+(IR-1)*4) 
     &              *DEXP(-WR3(IRATE+(IR-1)*4+1) *LNTHCLIP
     &              -WR3(IRATE+(IR-1)*4+2) /THCLIP)*WR1(IUNA)
     &              *WR1(IUNA)*1.D6*1.D6
            KFEQ   = KF(IR)
            CALL KEQUILIBRIUM (NS, NREA, GTHCLIP, WR3(ISTOIR),  
     &                       WR3(ISTOIP), WR1(IMI), RTHCLIP, LNRTHCLIP,
     &                       KEQ, IR)
            KB(IR) = KFEQ /KEQ
C     - Heavy particle impact (two body reaction)
          CASE(28)
            KF(IR) = WR3(IRATE+(IR-1)*4) 
     &              *DEXP(-WR3(IRATE+(IR-1)*4+1) *LNTHCLIP
     &              -WR3(IRATE+(IR-1)*4+2) /THCLIP)*WR1(IUNA)*1.D6
            KFEQ   = KF(IR)
            CALL KEQUILIBRIUM (NS, NREA, GTHCLIP, WR3(ISTOIR),  
     &            WR3(ISTOIP), WR1(IMI), RTHCLIP, LNRTHCLIP, KEQ, IR)
            KB(IR) = KFEQ /KEQ
          CASE(29)
C     - Heavy particle impact (three body reaction)
            KF(IR) = WR3(IRATE+(IR-1)*4) 
     &              *DEXP(-WR3(IRATE+(IR-1)*4+1) *LNTECLIP
     &              -WR3(IRATE+(IR-1)*4+2) /TECLIP)*WR1(IUNA)*1.D6
     &              *WR1(IUNA)*1.D6
            KFEQ   = KF(IR)
            CALL KEQUILIBRIUM (NS, NREA, GTECLIP, WR3(ISTOIR),  
     &           WR3(ISTOIP), WR1(IMI), RTECLIP, LNRTECLIP, KEQ, IR)
            KB(IR) = KFEQ /KEQ
          CASE(30)
C     - TITAN ONLY
            KF(IR) = WR3(IRATE+(IR-1)*4)
            CALL CNBOLTZMANN (TH, TH, TH, xCNXTg, xCNBTg, xCNATg)
            KB(IR) = KF(IR)/xCNATg*xCNXTg*DEXP(-4.*1.44 *2358.0 /TH)

          CASE(31)
C     - TITAN ONLY
           KF(IR) = WR3(IRATE+(IR-1)*4)
            CALL CNBOLTZMANN (TH, TH, TH, xCNXTg, xCNBTg, xCNATg)
           KB(IR) = KF(IR)/xCNBTg*xCNXTg*DEXP(-11.*1.44 *2358.0 /TH)

CC   - Knab Model for Dissociation !
C          CASE(32)
C            CALL DEVIATIONFACTOR4(WR1, LWR1, WR3, LWR3, WI, LWI, 
C     &                            WR3(IRATE+(IR-1)*4+2), TH, TV, PSI4,
C     &                            WI(ITNEQ+IR-1), IR, TVI)
C 
C                     KF(IR) = PSI4*WR3(IRATE+(IR-1)*4) 
C     &                        *DEXP(WR3(IRATE+(IR-1)*4+1)*LNTH
C     &                        -WR3(IRATE+(IR-1)*4+2)/TH) 
C     
C                     KFEQ   = WR3(IRATE+(IR-1)*4) 
C     &                        *DEXP(WR3(IRATE+(IR-1)*4+1)*LNTH
C     &                        -WR3(IRATE+(IR-1)*4+2)/TH)  
C
C            CALL KEQUILIBRIUM (NS, NREA, GTH, WR3(ISTOIR), WR3(ISTOIP),
C     &                         WR1(IMI), RTH, LNRTH, KEQ, IR)
C            KB(IR) = KFEQ/KEQ            
C
C-----------------------------------------------------------------------
C     Electronic excitation by heavy particle impact  (two body reaction)
C     Nitrogen excited by impact with Nitrogen atoms
C-----------------------------------------------------------------------
         CASE(33)
            CALL CIEALOURDS(WR1(IUKB), WR1(IUPI), WR1(IUNA), TH,
     &                     WR3(IRATE+(IR-1)*4+1), WR3(IRATE+(IR-1)*4+2),
     &                     WR1(IMI+1),WR1(IMI+1), KF(IR))
            KF(IR) = KF(IR)*WR1(IUNA)
            KFEQ = KF(IR)
            CALL KEQUILIBRIUM (NS, NREA, GTH, WR3(ISTOIR), WR3(ISTOIP), 
     &                        WR1(IMI), RTH, LNRTH, KEQ, IR)
            KB(IR) = KFEQ/KEQ
C-----------------------------------------------------------------------
C     Electronic excitation by heavy particle impact  (two body reaction)
C     Nitrogen excited by impact with oxygen atoms
C-----------------------------------------------------------------------
          CASE(34)
            CALL CIEALOURDS(WR1(IUKB), WR1(IUPI), WR1(IUNA), TH,
     &                     WR3(IRATE+(IR-1)*4+1), WR3(IRATE+(IR-1)*4+2),
     &                     WR1(IMI+1), WR1(IMI+2), KF(IR))
            KF(IR) = KF(IR)*WR1(IUNA)
            KFEQ = KF(IR)
            CALL KEQUILIBRIUM (NS, NREA, GTH, WR3(ISTOIR), WR3(ISTOIP), 
     &                        WR1(IMI), RTH, LNRTH, KEQ, IR)
            KB(IR) = KFEQ/KEQ
C-----------------------------------------------------------------------
C     Electronic excitation by heavy particle impact  (two body reaction)
C     Oxygen excited by impact with other oxygen atoms
C-----------------------------------------------------------------------
          CASE(35) 
            CALL CIEALOURDS(WR1(IUKB), WR1(IUPI), WR1(IUNA), TH,
     &                     WR3(IRATE+(IR-1)*4+1), WR3(IRATE+(IR-1)*4+2),
     &                     WR1(IMI+2), WR1(IMI+2), KF(IR))
            KF(IR) = KF(IR)*WR1(IUNA)
            KFEQ = KF(IR)
            CALL KEQUILIBRIUM (NS, NREA, GTH, WR3(ISTOIR), WR3(ISTOIP), 
     &                        WR1(IMI), RTH, LNRTH, KEQ, IR)
            KB(IR) = KFEQ/KEQ
C-----------------------------------------------------------------------
C     Electronic excitation by heavy particle impact  (two body reaction)
C     Oxygen excited by impact with nitrogen atoms
C-----------------------------------------------------------------------
          CASE(36)
            CALL CIEALOURDS(WR1(IUKB), WR1(IUPI), WR1(IUNA), TH,
     &                     WR3(IRATE+(IR-1)*4+1), WR3(IRATE+(IR-1)*4+2),
     &                     WR1(IMI+2), WR1(IMI+1), KF(IR))
            KF(IR) = KF(IR)*WR1(IUNA)
            KFEQ = KF(IR)
            CALL KEQUILIBRIUM (NS, NREA, GTH, WR3(ISTOIR), WR3(ISTOIP), 
     &                        WR1(IMI), RTH, LNRTH, KEQ, IR)
            KB(IR) = KFEQ/KEQ


C-----------------------------------------------------------------------
C     Electron impact excitation 
C     Formula given by Park (See the book)
C-----------------------------------------------------------------------
          CASE(37)
            KF(IR) = WR3(IRATE+(IR-1)*4)*DEXP(WR3(IRATE+(IR-1)*4+1)*
     &               DLOG(TE/10000.D0) - WR3(IRATE+(IR-1)*4+2)*
     &               1.43876866D0/TE)*WR1(IUNA)
            KFEQ   = KF(IR)
            CALL KEQUILIBRIUM (NS, NREA, GTE, WR3(ISTOIR), WR3(ISTOIP), 
     &                         WR1(IMI), RTE, LNRTE, KEQ, IR)
            KB(IR) = KFEQ /KEQ

c          CASE(38)
c
c           y = (WR3(IRATE+(IR-1)*4+2)-WR3(IRATE+(IR-1)*4+1))*
c     &              1.43876866D0/TE
c
c            KF(IR) = (UTE*4.D0*WR1(IUPI)*((0.529D-10)**2.D0)
c     &              *(157820.40447D0/TE)**2.D0 * Ia)*WR1(IUNA)
c            
c           PSIy =  DEXP(-y)/(1.+y)*(1./(20.+y)+DLOG(1.25*(1+1/y)))

         
c           KF(IR) = 1.46D-10*TE^0.5*(157820.40447D0/(WR3(IRATE+(IR-1)*4+2)-WR3(IRATE+(IR-1)*4+1))*())**2.D0

c            KFEQ   = KF(IR)

c            CALL KEQUILIBRIUM (NS, NREA, GTE, WR3(ISTOIR), WR3(ISTOIP), 
c     &                         WR1(IMI), RTE, LNRTE, KEQ, IR)
c            KB(IR) = KFEQ /KEQ

C-----------------------------------------------------------------------
C     Electron impact excitation 
C     Formula given by Gryzinski
C-----------------------------------------------------------------------

          CASE(40)

            KF(IR) = WR3(IRATE+(IR-1)*4)*DEXP(WR3(IRATE+(IR-1)*4+1)*
     &               DLOG(TE) - WR3(IRATE+(IR-1)*4+2)/TE)*WR1(IUNA) 
            KFEQ   = KF(IR)

            CALL KEQUILIBRIUM (NS, NREA, GTE, WR3(ISTOIR), WR3(ISTOIP), 
     &                         WR1(IMI), RTE, LNRTE, KEQ, IR)
            KB(IR) = KFEQ /KEQ

C-----------------------------------------------------------------------
C     Electron impact Ionization 
C     Formula given by Johnston Drawin 1968
C-----------------------------------------------------------------------
          CASE(41)

           fy = (WR3(IRATE+(IR-1)*4+2) - WR3(IRATE+(IR-1)*4+1))*
     &              1.43876866D0/TE

            
           psi = DEXP(-fy)/(1.D0+fy)*(1.D0/(20.D0 + fy)+
     &           DLOG(1.25*(1.+1./fy)))


            KF(IR) = (1.46D-10*sqrt(TE)*(157820.40447D0/
     &       ((WR3(IRATE+(IR-1)*4+2) - WR3(IRATE+(IR-1)*4+1))
     &        *1.43876866D0))**2.D0*fy*psi)/1.D6*WR1(IUNA)
 
            KFEQ   = KF(IR)

            CALL KEQUILIBRIUM (NS, NREA, GTE, WR3(ISTOIR), WR3(ISTOIP), 
     &                         WR1(IMI), RTE, LNRTE, KEQ, IR)
            KB(IR) = KFEQ /KEQ

            CASE(50)

              FLAGION = INT(WR3(IRATE+(IR-1)*4+2)) 

              IF (FLAGION.EQ.1) THEN
                CALL AVERAGEION(WR1,LWR1,TE,CRE1,CRG1,3,KF(IR))
              ELSE IF (FLAGION.EQ.2) THEN
                CALL AVERAGEION(WR1,LWR1,TE,CRE2,CRG2,7,KF(IR))
              ELSE IF (FLAGION.EQ.3) THEN
                CALL AVERAGEION(WR1,LWR1,TE,CRE3,CRG3,8,KF(IR))
              ELSE IF (FLAGION.EQ.4) THEN
                CALL AVERAGEION(WR1,LWR1,TE,CRE4,CRG4,6,KF(IR))
              ELSE IF (FLAGION.EQ.5) THEN
                CALL AVERAGEION(WR1,LWR1,TE,CRE5,CRG5,19,KF(IR))
              ENDIF

            KFEQ   = KF(IR)

            CALL KEQUILIBRIUM (NS, NREA, GTE, WR3(ISTOIR), WR3(ISTOIP), 
     &                         WR1(IMI), RTE, LNRTE, KEQ, IR)
            KB(IR) = KFEQ /KEQ
 
             CASE(51)

             FLAGI = INT(WR3(IRATE+(IR-1)*4+1)) 
             FLAGJ = INT(WR3(IRATE+(IR-1)*4+2)) 

            IF (FLAGI.EQ.1) THEN
              DO I=1,20
                    EAUX(I)=CRE1(I); GAUX(I)=CRG1(I)
                    LIAUX(I)=SQN1(I)
              ENDDO
              NA = 3
            ELSE IF (FLAGI.EQ.2) THEN
              DO I=1,20
                    EAUX(I)=CRE2(I); GAUX(I)=CRG2(I)
                    LIAUX(I)=SQN2(I)
              ENDDO
              NA = 7
            ELSE IF (FLAGI.EQ.3) THEN
              DO I=1,20
                    EAUX(I)=CRE3(I); GAUX(I)=CRG3(I)
                    LIAUX(I)=SQN3(I)
              ENDDO
              NA = 8
            ELSE IF (FLAGI.EQ.4) THEN
              DO I=1,20
                    EAUX(I)=CRE4(I); GAUX(I)=CRG4(I)
                    LIAUX(I)=SQN4(I)
              ENDDO
              NA = 6
            ELSE IF (FLAGI.EQ.5) THEN
              DO I=1,20
                    EAUX(I)=CRE5(I); GAUX(I)=CRG5(I)
                    LIAUX(I)=SQN5(I)
              ENDDO
              NA = 19
            ENDIF

            IF (FLAGJ.EQ.1) THEN
              DO I=1,20
                    EAUX2(I)=CRE1(I); GAUX2(I)=CRG1(I)
                    LJAUX(I)=SQN1(I)
              ENDDO
              NB = 3
            ELSE IF (FLAGJ.EQ.2) THEN
              DO I=1,20
                    EAUX2(I)=CRE2(I); GAUX2(I)=CRG2(I)
                    LJAUX(I)=SQN2(I)
              ENDDO
              NB = 7
            ELSE IF (FLAGJ.EQ.3) THEN
              DO I=1,20
                    EAUX2(I)=CRE3(I); GAUX2(I)=CRG3(I)
                    LJAUX(I)=SQN3(I)
              ENDDO
              NB = 8
            ELSE IF (FLAGJ.EQ.4) THEN
              DO I=1,20
                    EAUX2(I)=CRE4(I); GAUX2(I)=CRG4(I)
                    LJAUX(I)=SQN4(I)
              ENDDO
              NB = 6
            ELSE IF (FLAGJ.EQ.5) THEN
              DO I=1,20
                    EAUX2(I)=CRE5(I); GAUX2(I)=CRG5(I)
                    LJAUX(I)=SQN5(I)
              ENDDO
              NB = 19
            ENDIF

           CALL AVERAGEEXC(WR1,LWR1,TE,GAUX,GAUX2,LIAUX,LJAUX, 
     &                   EAUX,EAUX2,NA, NB,KF(IR))

            KFEQ   = KF(IR)

            CALL KEQUILIBRIUM (NS, NREA, GTE, WR3(ISTOIR), WR3(ISTOIP), 
     &                         WR1(IMI), RTE, LNRTE, KEQ, IR)
            KB(IR) = KFEQ /KEQ

            CASE(52)
              
              FLAGION = INT(WR3(IRATE+(IR-1)*4+2)) 
              EIKELVIN= WR3(IRATE+(IR-1)*4+1)*1.43876866D0

              IF (FLAGION.EQ.1) THEN
                CALL AVERAGEEXC2(WR1,LWR1,TE,CRG1,EIKELVIN,
     &                         SQN1,CRE1,3,KF(IR))
             ELSE IF (FLAGION.EQ.2) THEN
               CALL AVERAGEEXC2(WR1,LWR1,TE,CRG2,EIKELVIN,
     &                         SQN2,CRE2,7,KF(IR))
             ELSE IF (FLAGION.EQ.3) THEN
               CALL AVERAGEEXC2(WR1,LWR1,TE,CRG3,EIKELVIN,
     &                         SQN3,CRE3,8,KF(IR))
             ELSE IF (FLAGION.EQ.4) THEN
               CALL AVERAGEEXC2(WR1,LWR1,TE,CRG4,EIKELVIN,
     &                         SQN4,CRE4,6,KF(IR))
             ELSE IF (FLAGION.EQ.5) THEN
               CALL AVERAGEEXC2(WR1,LWR1,TE,CRG5,EIKELVIN,
     &                         SQN5,CRE5,19,KF(IR))
             ENDIF

           KFEQ   = KF(IR)

            CALL KEQUILIBRIUM (NS, NREA, GTE, WR3(ISTOIR), WR3(ISTOIP), 
     &                         WR1(IMI), RTE, LNRTE, KEQ, IR)
            KB(IR) = KFEQ /KEQ
C-----------------------------------------------------------------------
C     60. Associative ionization reactions 
C     =>  Losev's formulation for KF(TH,TV), KB(TH,TV)
C-----------------------------------------------------------------------
          CASE(60)
            KF(IR) = WR3(IRATE+(IR-1)*4)*
     &               (WR3(IRATE+(IR-1)*4+2)*TH + 1.D0)*
     &                DEXP(WR3(IRATE+(IR-1)*4+1)*LNTH -
     &                WR3(IRATE+(IR-1)*4+3)/TH)
           
            CALL KEQUILIBRIUM (NS, NREA, GTH, WR3(ISTOIR),  WR3(ISTOIP),
     &                        WR1(IMI), RTH, LNRTH, KEQ, IR)
            KB(IR) = KF(IR)/KEQ 
C-----------------------------------------------------------------------
C     -61. Dissociation by atomic or molecular impact / recombination 
C          Knab et al. formulation for KF(TH,TV), KB(TH,TV)
C-----------------------------------------------------------------------
          CASE(61)        
            ISDIS = WI(IDISSOCIATION+IR-1)

            CALL KNAB(WR1,LWR1,WR3,LWR3,WI,LWI,IR,ISDIS,WI(ITNEQ+IR-1),
     &               ALPHA,PAR, WR3(IRATE+(IR-1)*4+2),TEF,PSIK)
            
            KFEQ = WR3(IRATE+(IR-1)*4) 
     &             *DEXP(WR3(IRATE+(IR-1)*4+1)*LNTH
     &             -WR3(IRATE+(IR-1)*4+2)/TH)

            KF(IR) = KFEQ*PSIK(1)                               
         
            CALL KEQUILIBRIUM (NS, NREA, GTH, WR3(ISTOIR), WR3(ISTOIP), 
     &                         WR1(IMI), RTH, LNRTH, KEQ, IR)
            KB(IR) = KFEQ/KEQ
C-----------------------------------------------------------------------
C     -62. Radical reactions (including Zeldovich reactions)
C          Knab et al. formulation for KF(TH,TV), KB(TH,TV)
C-----------------------------------------------------------------------
          CASE(62)            
            CALL KNAB(WR1,LWR1,WR3,LWR3,WI,LWI,IR,ISDIS,
     &              WI(ITNEQ+IR-1), ALPHA,PAR, 
     &              WR3(IRATE+(IR-1)*4+2),TEF,PSIK)
                     
            KFEQ = WR3(IRATE+(IR-1)*4) 
     &           *DEXP(WR3(IRATE+(IR-1)*4+1)*LNTH
     &           -WR3(IRATE+(IR-1)*4+2)/TH)

            KF(IR) = PSIK(2)*KFEQ
            
            CALL KEQUILIBRIUM (NS, NREA, GTH, WR3(ISTOIR),  WR3(ISTOIP),
     &                        WR1(IMI), RTH, LNRTH, KEQ, IR)

            KB(IR) = PSIK(3)*KFEQ/KEQ
C-----------------------------------------------------------------------
C     -63. Associative ionization reactions
C          Knab et al. formulation for KF(TH,TV), KB(TE,TV)
C-----------------------------------------------------------------------
         CASE(63)            
           CALL KNAB(WR1,LWR1,WR3,LWR3,WI,LWI,IR,ISDIS,
     &               WI(ITNEQ+IR-1), ALPHA,PAR, 
     &               WR3(IRATE+(IR-1)*4+2),TEF,PSIK)

           KF(IR) = WR3(IRATE+(IR-1)*4) 
     &             *DEXP(WR3(IRATE+(IR-1)*4+1)*LNTH
     &             -WR3(IRATE+(IR-1)*4+2)/TH)
           
           KFEQ = WR3(IRATE+(IR-1)*4) 
     &            *DEXP(WR3(IRATE+(IR-1)*4+1)*LNTE
     &            -WR3(IRATE+(IR-1)*4+2)/TE)

           CALL KEQUILIBRIUM (NS, NREA, GTE, WR3(ISTOIR), WR3(ISTOIP), 
     &                       WR1(IMI), RTE, LNRTE, KEQ, IR)
          
           KB(IR) = PSIK(4)*KFEQ/KEQ
C-----------------------------------------------------------------------
C     -64. Electron impact dissociation 
C          Knab et al. formulation for KF(TE,TV), KB(TE,TV)
C-----------------------------------------------------------------------
         CASE(64)
            CALL KNAB(WR1,LWR1,WR3,LWR3,WI,LWI,IR,ISDIS,
     &               WI(ITNEQ+IR-1),ALPHA,PAR, 
     &               WR3(IRATE+(IR-1)*4+2),TEF,PSIK)     
             
            KFEQ   = WR3(IRATE+(IR-1)*4) 
     &              *DEXP(WR3(IRATE+(IR-1)*4+1)*LNTE
     &              -WR3(IRATE+(IR-1)*4+2)/TE)

            KF(IR) = KFEQ*PSIK(5)

            CALL KEQUILIBRIUM (NS, NREA, GTE, WR3(ISTOIR), WR3(ISTOIP), 
     &                         WR1(IMI), RTE, LNRTE, KEQ, IR)
            KB(IR) = KFEQ/KEQ
C-----------------------------------------------------------------------
C     -65. Charge exchange reaction 
C          Knab et al. formulation for KF(TH,TV), KB(TH,TV)
C-----------------------------------------------------------------------
         CASE(65)
            CALL KNAB(WR1,LWR1,WR3,LWR3,WI,LWI,IR,ISDIS,
     &                WI(ITNEQ+IR-1),ALPHA,PAR, 
     &                WR3(IRATE+(IR-1)*4+2),TEF,PSIK)

            KF(IR) = WR3(IRATE+(IR-1)*4) 
     &              *DEXP(WR3(IRATE+(IR-1)*4+1)*LNTH
     &              -WR3(IRATE+(IR-1)*4+2)/TH)

            KFEQ = KF(IR)

            CALL KEQUILIBRIUM (NS, NREA, GTH, WR3(ISTOIR), WR3(ISTOIP), 
     &                        WR1(IMI), RTH, LNRTE, KEQ, IR)
         
            KB(IR) = PSIK(6)*KFEQ/KEQ 
C-----------------------------------------------------------------------
C     -66. Dissociation / Recombination reaction
C          Treanor and Marrone formulation for  KF(TH,TV), KB(TH)
C-----------------------------------------------------------------------
         CASE(66)
            ISDIS = WI(IDISSOCIATION+IR-1)
            TVI = TEF(WI(IVIBTEMPI+ISDIS-1)) 
             
            CALL TREANOR(WR1,LWR1,WI,LWI,ISDIS,ALPHA,6.D0,
     &                   WR3(IRATE+(IR-1)*4+2),TH,TVI,PSITM)  
          
            KFEQ   = WR3(IRATE+(IR-1)*4) 
     &              *DEXP(WR3(IRATE+(IR-1)*4+1)*LNTH
     &              -WR3(IRATE+(IR-1)*4+2)/TH)
       
            KF(IR) = PSITM*KFEQ   

            CALL KEQUILIBRIUM (NS, NREA, GTH, WR3(ISTOIR), WR3(ISTOIP), 
     &                        WR1(IMI), RTH, LNRTH, KEQ, IR)

            KB(IR) = KFEQ/KEQ
C-----------------------------------------------------------------------
C-67,68,69.Dissociation / Recombination reaction
C          Losev formulation for  KF(TH,TV), KB(TH)
C          This model is label in literature as "Beta" model and represents
C          a particular case of that of Kutznetsov 
C-----------------------------------------------------------------------
         CASE(67,68,69)
            ISDIS = WI(IDISSOCIATION+IR-1)

            CALL LOSEV (WR1,LWR1,WR3,LWR3,WI,LWI,TEF,ISDIS,
     &                 WI(ITNEQ+IR-1),IR,PHI,KFEQ)
               
            KFEQ = WR3(IRATE+(IR-1)*4)*
     &             (1.D0 - EXP(-WR3(IRATE+(IR-1)*4+1)/TH))*
     &             DEXP(-WR3(IRATE+(IR-1)*4+2)/TH)

            KF(IR) = PHI*KFEQ   
          
            CALL KEQUILIBRIUM (NS, NREA, GTH, WR3(ISTOIR), WR3(ISTOIP), 
     &                      WR1(IMI), RTH, LNRTH, KEQ, IR)
            KB(IR) = KFEQ/KEQ   
C-----------------------------------------------------------------------
C   -70,71.Exchange reactions involving the formation of nitric oxyde
C          The model for this kind of chemical reactions is illustrated 
C          in Losev paper and comes from the study made by Levitsky 
C          based on quasi-classical trajectory method computation
C          N2(v) + O = NO + N
C          O2(v) + N = NO + O 
C-----------------------------------------------------------------------
         CASE(70,71)   
            ISDIS = 0.D0

            CALL LOSEV (WR1,LWR1,WR3,LWR3,WI,LWI,TEF,ISDIS,
     &                 WI(ITNEQ+IR-1),IR,PHI,KFEQ)
          
            KFEQ = WR3(IRATE+(IR-1)*4)*KFEQ
            KF(IR) = PHI*KFEQ
 
            CALL KEQUILIBRIUM (NS, NREA, GTH, WR3(ISTOIR), WR3(ISTOIP), 
     &                      WR1(IMI), RTH, LNRTH, KEQ, IR)
            KB(IR) = KFEQ/KEQ   
C-----------------------------------------------------------------------
C     -72 Hammerling formulation for dissociation-recombination reactions  
C-----------------------------------------------------------------------
         CASE(72)
            ISDIS = WI(IDISSOCIATION+IR-1)
            TVI = TEF(WI(IVIBTEMPI+ISDIS-1))

            CALL HAMMERLING (WR1,LWR1,WI,LWI,ISDIS,
     &                      WR3(IRATE+(IR-1)*4+2),TH,TVI,PSIHAM)

            KFEQ   = WR3(IRATE+(IR-1)*4)
     &              *DEXP(WR3(IRATE+(IR-1)*4+1)*LNTH
     &              -WR3(IRATE+(IR-1)*4+2)/TH)     

            KF(IR) = PSIHAM*KFEQ

            CALL KEQUILIBRIUM (NS, NREA, GTH, WR3(ISTOIR), WR3(ISTOIP),
     &                        WR1(IMI), RTH, LNRTH, KEQ, IR)

            KB(IR) = KFEQ/KEQ
C-----------------------------------------------------------------------
C     -73 Hassan formulation for dissociation-recombination reactions  
C-----------------------------------------------------------------------
         CASE(73)
            ISDIS = WI(IDISSOCIATION+IR-1)
            TVI = TEF(WI(IVIBTEMPI+ISDIS-1))

            CALL HASSAN (WR1,LWR1,WI,LWI,ISDIS,
     &                   WR3(IRATE+(IR-1)*4+2),TH,TVI,PSIHAS)

            KFEQ   = WR3(IRATE+(IR-1)*4)
     &              *DEXP(WR3(IRATE+(IR-1)*4+1)*LNTH
     &              -WR3(IRATE+(IR-1)*4+2)/TH)

            KF(IR) = PSIHAS*KFEQ

            CALL KEQUILIBRIUM (NS, NREA, GTH, WR3(ISTOIR), WR3(ISTOIP),
     &                        WR1(IMI), RTH, LNRTH, KEQ, IR)
            KB(IR) = KFEQ/KEQ
C-----------------------------------------------------------------------
C     -74 Hansen formulation for dissociation-recombination reactions  
C-----------------------------------------------------------------------
         CASE(74)
            ISDIS = WI(IDISSOCIATION+IR-1)
            TVI = TEF(WI(IVIBTEMPI+ISDIS-1))  

            CALL HANSEN (ISDIS,WR3(IRATE+(IR-1)*4+2),TH,TVI,PSIHANS)

            KFEQ   = WR3(IRATE+(IR-1)*4)
     &              *DEXP(WR3(IRATE+(IR-1)*4+1)*LNTH
     &              -WR3(IRATE+(IR-1)*4+2)/TH) 

            KF(IR) = PSIHANS*KFEQ
             
            CALL KEQUILIBRIUM (NS, NREA, GTH, WR3(ISTOIR), WR3(ISTOIP),
     &                        WR1(IMI), RTH, LNRTH, KEQ, IR)
            KB(IR) = KFEQ/KEQ 
C-----------------------------------------------------------------------
C     -75,76,77 Macheret formulation for dissociation-recombination 
C         reactions  
C-----------------------------------------------------------------------
         CASE(75,76,77)
            ISDIS = WI(IDISSOCIATION+IR-1)
            TVI = TEF(WI(IVIBTEMPI+ISDIS-1)) 
         
            CALL MACHERET (WR1,LWR1,WI,LWI,ISDIS,TH,TVI,
     &                    WR3(IRATE+(IR-1)*4+2),WR3(IRATE+(IR-1)*4+1),
     &                    WI(ITNEQ+IR-1),PSIMACH)  

            KFEQ   = WR3(IRATE+(IR-1)*4)
     &              *DEXP(WR3(IRATE+(IR-1)*4+1)*DLOG(TH)
     &              -WR3(IRATE+(IR-1)*4+2)/TH) 

            KF(IR) = PSIMACH*KFEQ

            CALL KEQUILIBRIUM (NS, NREA, GTH, WR3(ISTOIR), WR3(ISTOIP),
     &                        WR1(IMI), RTH, LNRTH, KEQ, IR)
            KB(IR) = KFEQ/KEQ

          CASE DEFAULT
            WRITE(*,*) 'Mechanism not implemented'  
            STOP
        ENDSELECT
      ENDDO

C     Reaction source terms [mole /(m^3 s)]
      OMEGAI = 0.D0
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
        SELECT CASE(WI(ITNEQ+IR-1))
          CASE(2)
!            OMEGAI = OMEGAI - WR3(IRATE+(IR-1)*4+2)*WR1(IUR)*XI(IR)

          CASE(3)
            OMEGAI = OMEGAI - WR3(IRATE+(IR-1)*4+2)*WR1(IUR)*XI(IR)
C   TO BE ADDED CASE 8 and 9
          CASE(8)
            OMEGAI = OMEGAI - XI(IR)*(WR3(IRATE+(IR-1)*4+2) -
     &                         WR3(IRATE+(IR-1)*4+1))*1.43876866D0*
     &                         WR1(IUR)
          CASE(9)
            OMEGAI = OMEGAI - XI(IR)*(WR3(IRATE+(IR-1)*4+2) -
     &                         WR3(IRATE+(IR-1)*4+1))*1.43876866D0*
     &                         WR1(IUR)
          CASE(10)
            OMEGAI = OMEGAI - XI(IR)*(WR3(IRATE+(IR-1)*4+2) -
     &                         WR3(IRATE+(IR-1)*4+1))*1.43876866D0*
     &                         WR1(IUR)


         CASE(50)
            ENU = 117345.
            FLAGI =  INT(WR3(IRATE+(IR-1)*4+2))

            IF (FLAGI.EQ.1) THEN
                ENL = 86197.
            ELSEIF (FLAGI.EQ.2) THEN
                ENL = 96755.
            ELSEIF (FLAGI.EQ.3) THEN
                ENL = 104723.
            ELSEIF (FLAGI.EQ.4) THEN
                ENL = 107224.
            ELSEIF (FLAGI.EQ.5) THEN
                ENL = 110079.
            ENDIF

            OMEGAI = OMEGAI - XI(IR)*(ENU - ENL)*1.43876866D0*
     &                         WR1(IUR)

          CASE(51)
            FLAGI =  INT(WR3(IRATE+(IR-1)*4+1))
            FLAGJ =  INT(WR3(IRATE+(IR-1)*4+2))

            IF (FLAGI.EQ.1) THEN
                ENL = 86197.
            ELSEIF (FLAGI.EQ.2) THEN
                ENL = 96755.
            ELSEIF (FLAGI.EQ.3) THEN
                ENL = 104723.
            ELSEIF (FLAGI.EQ.4) THEN
                ENL = 107224.
            ELSEIF (FLAGI.EQ.5) THEN
                ENL = 110079.
            ENDIF

            IF (FLAGJ.EQ.1) THEN
                ENU = 86197.
            ELSEIF (FLAGJ.EQ.2) THEN
                ENU = 96755.
            ELSEIF (FLAGJ.EQ.3) THEN
                ENU = 104723.
            ELSEIF (FLAGJ.EQ.4) THEN
                ENU = 107224.
            ELSEIF (FLAGJ.EQ.5) THEN
                ENU = 110079.
            ENDIF

            OMEGAI = OMEGAI - XI(IR)*(ENU - ENL)*1.43876866D0*
     &                         WR1(IUR)
          CASE(52)
            ENL =  WR3(IRATE+(IR-1)*4+1)
            FLAGJ =  INT(WR3(IRATE+(IR-1)*4+2))

            IF (FLAGJ.EQ.1) THEN
                ENU = 86197.
            ELSEIF (FLAGJ.EQ.2) THEN
                ENU = 96755.
            ELSEIF (FLAGJ.EQ.3) THEN
                ENU = 104723.
            ELSEIF (FLAGJ.EQ.4) THEN
                ENU = 107224.
            ELSEIF (FLAGJ.EQ.5) THEN
                ENU = 110079.
            ENDIF

            OMEGAI = OMEGAI - XI(IR)*(ENU - ENL)*1.43876866D0*
     &                         WR1(IUR)

!          CASE(37)
!            OMEGAI = OMEGAI - XI(IR)*WR3(IRATE+(IR-1)*4+2)
!     &                         *1.43876866D0*WR1(IUR)
        END SELECT
      ENDDO
      DO IS = 1, NS
        OMEGA(IS) = 0.D0
        DO IR = 1, NREA
          OMEGA(IS) = OMEGA(IS) +(WR3(ISTOIP+(IR-1)*NS+IS-1)
     &                            -WR3(ISTOIR+(IR-1)*NS+IS-1)) *XI(IR)
       ENDDO
      ENDDO
      IC = 0
      DO I = 1, NVIB
        OMEGACV(I)=0.D0
        DO J = 1, WI(INVIBSPEI+I-1)
          IC = IC +1
          IVT = WI(IVIBSPEI+IC-1)

        CALL CVTRANSFER (LWR1, WR1, LWR3, WR3, 
     &                   LWI, WI, IVT, WR1(IMI+IVT-1), KF, 
     &                   KB, XI, TEF, ALPHA, PAR, MOLE, LNMOLE, CV) 
        OMEGACV(I)=OMEGACV(I)+CV
       ENDDO
      ENDDO

C-----------------------------------------------------------------------
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


      SUM = 0.D0
      DO IS = 1, NS
        SUM = SUM +(STOIP((IR-1)*NS+IS) -STOIR((IR-1)*NS+IS)) 
     &        *(G(IS) +RTLNRT)
      ENDDO
      KEQ = DEXP(-SUM /RT)+1.D-50

      END SUBROUTINE  KEQUILIBRIUM
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      SUBROUTINE READRATES (MIXTURE, PATH, REACTION, WI, LWI, WC, LWC, 
     &                      WR1, LWR1, WR3, LWR3)
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
      INTEGER LWI, WI(1:LWI), LWR1, LWR3, LWC, NCR, IC1, IC2, LSPECIES
      INTEGER II, ISPEC
      DOUBLE PRECISION WR1(1:LWR1), WR3(1:LWR3)
      CHARACTER(*) REACTION, MIXTURE
      CHARACTER(*) PATH
      CHARACTER WC(1:LWC)
C-----------------------------------------------------------------------
      INTEGER LPATH, LREACTION, ILEN, I, IC, IM, ITYPE, IARRHENIUS, IR, 
     &        IS, J, LCHAR, ISPECIES, NTHERMAL, IS1, IS2, IS3, IS4, IS5,
     &        ISDIS, JR, JJ, N1, N2, IA1, IA2, NALPHA(1:NREA), 
     &        ALPHA(1:NREA,1:NS), LMIXTURE, ICR, IAI, DCR
      DOUBLE PRECISION COEFFICIENT, ARRHENIUS(1:5),DUM3, 
     &                 REACTANT(1:NS), PRODUCT(1:NS),
     &                 ELEMENTR(1:NC), ELEMENTP(1:NC), SUM, DUM2,
     &                 VEC1(0:NS), VEC2(0:NS), VEC3(0:NS), VEC4(0:NS),
     &                 NVEC1, NVEC2, DIF1, DIF2, DIF, REACTTEMP,DUM1
      LOGICAL PROP
      CHARACTER(4) COM
      CHARACTER(10) WORD, SPECIES, DUM, IONSPEC
      CHARACTER(100) FULLCOM, LINE
      CHARACTER(80) COM1, COM2, FULLCOM1, FULLCOM2
C-----------------------------------------------------------------------
      LPATH     = LCHAR (PATH)
      LREACTION = LCHAR (REACTION)
      LMIXTURE  = LCHAR (MIXTURE)

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
              CALL DTRANSLATE (WORD, ARRHENIUS(IARRHENIUS))
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

C     Selection of the molecule dissociated by atomic or molecular 
C     impact in a reaction of dissociation. This molecule is
C     denoted here by IDIS and stored in WI(IDISSOCIATION+IR-1).
      DO IR = 1, NREA
        IF (WI(ITNEQ+IR-1).EQ.1.OR.WI(ITNEQ+IR-1).EQ.60.OR.
     &     WI(ITNEQ+IR-1).EQ.61.OR.WI(ITNEQ+IR-1).EQ.64.OR.
     &     WI(ITNEQ+IR-1).EQ.66.OR.WI(ITNEQ+IR-1).EQ.67.OR.
     &     WI(ITNEQ+IR-1).EQ.68.OR.WI(ITNEQ+IR-1).EQ.69.OR.
     &     WI(ITNEQ+IR-1).EQ.72.OR.WI(ITNEQ+IR-1).EQ.73.OR.
     &     WI(ITNEQ+IR-1).EQ.74.OR.WI(ITNEQ+IR-1).EQ.75.OR.
     &     WI(ITNEQ+IR-1).EQ.76.OR.WI(ITNEQ+IR-1).EQ.77)THEN
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

      
C    Reaction Enthalpies for all the reactions (Marco 12/11/2008)
      DO IR = 1, NREA
        WR3(IREACTENTH+IR -1) = 0.D0
        REACTTEMP = 0.D0
        DO IS = 1, NS
          REACTTEMP = REACTTEMP + WR3(ISTOIP+(IR-1)*NS+IS-1)
     &               * WR1(IHFORI+IS-1) - WR3(ISTOIR+(IR-1)*NS+IS-1)
     &               * WR1(IHFORI+IS-1)
        ENDDO
        WR3(IREACTENTH+IR -1) = REACTTEMP
      ENDDO
CC     Marco Check this out
C     Mixture file
      ICR = 0; NCR = 0;
      OPEN(UNIT=INOUT1,FILE=PATH(1:LPATH)//'/data/mixture/'
     &     //MIXTURE(1:LMIXTURE)//'.mix',STATUS='OLD')
      COM1 = '   '
      IC1 = 0; II = 0; IC2 = 0; IAI=0;
      DCR = 0
C     1. Nuclei
      DO WHILE (COM1(1:4)/= 'STOP')
        READ(INOUT1,*) FULLCOM1
        COM1 = FULLCOM1(1:4)
        IF (COM1(1:4) == 'Spec') THEN
          DO I = 1, NC
            READ(INOUT1,1003) SPECIES 
            LSPECIES  = LCHAR (SPECIES)
            IF (LSPECIES>5) THEN
              J = LSPECIES-4
              IF (SPECIES(J:LSPECIES)=='-star') THEN 
C     Treatment of the ground state of electronic specific nuclei
                OPEN(UNIT=INOUT2,FILE=PATH(1:LPATH)//'/data/thermo/'//
     &               SPECIES(1:J-1),STATUS='OLD')
                COM2 = '   '
                DO WHILE (COM2(1:4)/= 'STOP')
                  READ(INOUT2,*) FULLCOM2
                  COM2 = FULLCOM2(1:4)
                  IF (COM2(1:4)=='Elec') THEN
                    READ(INOUT2,*) NCR
                    READ(INOUT2,*) 
                    DO JJ=1,NCR
                      READ(INOUT2,*)DUM, WR3(IION+I+JJ-1-NE-1+IC1)
C                     cm-1 to --> K
                      WR3(IION+I+JJ-1-NE+IC1-1)= 
     &                          WR3(IION+I+JJ-1-NE-1+IC1)*1.43876866D0
C                      WI(INEUION+JJ-1) = ISPEC
                    ENDDO
                  ENDIF
                ENDDO
                REWIND(INOUT2)
                COM2 = '   '
                DO WHILE (COM2(1:4)/= 'STOP')
                  READ(INOUT2,*) FULLCOM2
                  COM2 = FULLCOM2(1:4)
                  IF (COM2(1:4)=='Ioni') THEN
                    READ(INOUT2,*)DUM1
                    DO JJ=1,NCR
                     WR3(IION+I+JJ-1-NE+IC1-1) = DUM1 - 
     &                                     WR3(IION+I+JJ-1-NE-1+IC1)
                     WR3(IION+I+JJ-1-NE+IC1-1) = 
     &                            WR3(IION+I+JJ-1-NE-1+IC1)*WR1(IUR)
                    ENDDO
                  ENDIF
                ENDDO
                CLOSE(INOUT2)
                IC1 = IC1 + NCR -1 
                DCR = DCR + NCR -1 
              ELSE 
                OPEN(UNIT=INOUT2,FILE=PATH(1:LPATH)//'/data/thermo/'//
     &               SPECIES(1:LSPECIES),STATUS='OLD')
                COM2 = '   '
                DO WHILE (COM2(1:4)/= 'STOP')
                  READ(INOUT2,*) FULLCOM2
                  COM2 = FULLCOM2(1:4)
                  IF (COM2(1:4)=='Ioni') THEN
                    READ(INOUT2,*)WR3(IION+I-1-NE+IC1)
                  ENDIF
                ENDDO
                CLOSE(INOUT2)
              ENDIF
            ELSE
              IF (I.GT.NE) THEN
                OPEN(UNIT=INOUT2,FILE=PATH(1:LPATH)//'/data/thermo/'//
     &                 SPECIES(1:LSPECIES),STATUS='OLD')
                COM2 = '   '
                DO WHILE (COM2(1:4)/= 'STOP')
                  READ(INOUT2,*) FULLCOM2
                  COM2 = FULLCOM2(1:4)
                  IF (COM2(1:4)=='Ioni') THEN
                    READ(INOUT2,*)WR3(IION+I-1-NE+IC1)
                    WR3(IION+I-1-NE+IC1) = WR3(IION+I-1-NE+IC1)*WR1(IUR)
                  ENDIF
                ENDDO
                CLOSE(INOUT2)
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      REWIND(UNIT=INOUT1)
 
      
C     2. Other species
      CALL BLANK(COM1)
      DO WHILE (COM1(1:4)/= 'STOP')
        CALL BLANK(FULLCOM1)
        READ(INOUT1,*) FULLCOM1
        COM1 = FULLCOM1(1:4)
        IF (COM1(1:4) == 'Spec') THEN
          DO I = 1, NC
            READ(INOUT1,1003) SPECIES 
          ENDDO

          DO I = 1+NC, NS-NSTAR
            IF (WI(IQI+I+DCR-1).EQ.0) THEN
              CALL BLANK(SPECIES)
              READ(INOUT1,1003) SPECIES 
              LSPECIES  = LCHAR (SPECIES)
              IF (LSPECIES>5) THEN
                J = LSPECIES-4
                IF (SPECIES(J:LSPECIES)=='-star') THEN 
C     Excited specific species
                  OPEN(UNIT=INOUT2,FILE=PATH(1:LPATH)//'/data/thermo/'//
     &                 SPECIES(1:J-1),STATUS='OLD')
                  COM2 = '   '
                  DO WHILE (COM2(1:4)/= 'STOP')
                    READ(INOUT2,*) FULLCOM2
                    COM2 = FULLCOM2(1:4)
                    IF (COM2(1:4)=='Elec') THEN
                      READ(INOUT2,*) NCR
                      READ(INOUT2,*) 
                      DO JJ=1,NCR
                        READ(INOUT2,*)DUM, DUM2
C                       cm-1 to --> K
                        WR3(IION+I+JJ-1-NE+IC1-1)= 
     &                            DUM2*1.43876866D0
                      ENDDO
                    ENDIF
                  ENDDO
                  REWIND(INOUT2)
                  COM2 = '   '
                  DO WHILE (COM2(1:4)/= 'STOP')
                    READ(INOUT2,*) FULLCOM2
                    COM2 = FULLCOM2(1:4)
                    IF (COM2(1:4)=='Ioni') THEN
                      READ(INOUT2,*)DUM1
                      DO JJ=1,NCR
                       WR3(IION+I+JJ-1-NE+IC1-1) = DUM1 - 
     &                                       WR3(IION+I+JJ-1-NE-1+IC1)
                       WR3(IION+I+JJ-1-NE+IC1-1) = 
     &                              WR3(IION+I+JJ-1-NE-1+IC1)*WR1(IUR)
                      ENDDO
                    ENDIF
                  ENDDO
                  IF (WI(IATOMI+DCR+I-1).GT.1) THEN
                    REWIND(INOUT2)
                    COM2 = '   '
                    DO WHILE (COM2(1:4)/= 'STOP')
                      READ(INOUT2,*) FULLCOM2
                      COM2 = FULLCOM2(1:4)
                      IF (COM2(1:4)=='Elec') THEN
                        READ(INOUT2,*)NCR
                        READ(INOUT2,*)
                        DO JJ=1,NCR
                          READ(INOUT2,*) DUM, DUM, DUM, DUM3
                          WR3(IDIS+I+JJ-2+IC2-NC-IAI) = DUM3 
     &                                               *1.43876866D0
                          WR3(IDIS+I+JJ-2+IC2-NC-IAI) = 
     &                             WR3(IDIS+I+JJ-2+IC2-NC-IAI)*WR1(IUR)
                        ENDDO
                      ENDIF
                    ENDDO
                    IC2 = IC2 + NCR -1
                  ENDIF
                  CLOSE(INOUT2)
                  IC1 = IC1 + NCR -1 
                  DCR = DCR + NCR -1 
                ELSE 
                  OPEN(UNIT=INOUT2,FILE=PATH(1:LPATH)//'/data/thermo/'//
     &                 SPECIES(1:LSPECIES),STATUS='OLD')
                  COM2 = '   '
                  DO WHILE (COM2(1:4)/= 'STOP')
                    READ(INOUT2,*) FULLCOM2
                    COM2 = FULLCOM2(1:4)
                    IF (COM2(1:4)=='Ioni') THEN
                      READ(INOUT2,*)WR3(IION+I-1-NE+IC1)
                    ENDIF
                  ENDDO
                  IF (WI(IATOMI+DCR+I-1).GT.1) THEN
                    REWIND(INOUT2)
                    COM2 = '   '
                    DO WHILE (COM2(1:4)/= 'STOP')
                      READ(INOUT2,*) FULLCOM2
                      COM2 = FULLCOM2(1:4)
                      IF (COM2(1:4)=='Diss') THEN
                        READ(INOUT2,*) WR3(IDIS+I-1+IC2-NC-IAI)
                            WR3(IDIS+I-1+IC2-NC-IAI) = 
     &                             WR3(IDIS+I-1+IC2-NC-IAI)*WR1(IUR)
                      ENDIF
                    ENDDO
                    ENDIF
                    CLOSE(INOUT2)
                 ENDIF
              ELSE
                OPEN(UNIT=INOUT2,FILE=PATH(1:LPATH)//'/data/thermo/'//
     &                 SPECIES(1:LSPECIES),STATUS='OLD')
                COM2 = '   '
                DO WHILE (COM2(1:4)/= 'STOP')
                  READ(INOUT2,*) FULLCOM2
                  COM2 = FULLCOM2(1:4)
                  IF (COM2(1:4)=='Ioni') THEN
                    READ(INOUT2,*)WR3(IION+I-1-NE+IC1)
                    WR3(IION+I-1-NE+IC1) = WR3(IION+I-1-NE+IC1)*WR1(IUR)
                  ENDIF
                ENDDO

                IF (WI(IATOMI+DCR+I-1).GT.1) THEN
                  REWIND(INOUT2)
                  COM2 = '   '
                  DO WHILE (COM2(1:4)/= 'STOP')
                    READ(INOUT2,*) FULLCOM2
                    COM2 = FULLCOM2(1:4)
                    IF (COM2(1:4)=='Diss') THEN
                      READ(INOUT2,*) WR3(IDIS+I-1+IC2-NC-IAI)
                          WR3(IDIS+I-1+IC2-NC-IAI) = 
     &                             WR3(IDIS+I-1+IC2-NC-IAI)*WR1(IUR)
                    ENDIF
                  ENDDO
                ELSE
                  IAI = IAI +1
                ENDIF
                CLOSE(INOUT2)
              ENDIF
            ELSE
C           Ionized Species 
            CALL BLANK(SPECIES)
            READ(INOUT1,1003) SPECIES 
            LSPECIES  = LCHAR (SPECIES)
            IF (WI(IATOMI+DCR+I-1).GT.1) THEN
                IF (LSPECIES>5) THEN
                  J = LSPECIES-4
                  IF (SPECIES(J:LSPECIES)=='-star') THEN 
C     Excited specific species
                  OPEN(UNIT=INOUT2,FILE=PATH(1:LPATH)//'/data/thermo/'//
     &                   SPECIES(1:J-1),STATUS='OLD')
                    COM2 = '   '
                    DO WHILE (COM2(1:4)/= 'STOP')
                      READ(INOUT2,*) FULLCOM2
                      COM2 = FULLCOM2(1:4)
                      IF (COM2(1:4)=='Elec') THEN
                        READ(INOUT2,*) NCR
                        READ(INOUT2,*) 
                        DO JJ=1,NCR
                          READ(INOUT2,*) DUM, DUM, DUM, DUM3
                          WR3(IDIS+I+JJ-2+IC2-NC-IAI) = DUM3 
     &                                               *1.43876866D0
                          WR3(IDIS+I+JJ-2+IC2-NC-IAI) = 
     &                             WR3(IDIS+I+JJ-2+IC2-NC-IAI)*WR1(IUR)
                        ENDDO
                      ENDIF
                    ENDDO
                    IC2 = IC2 + NCR -1
                    DCR = DCR + NCR -1 
                    CLOSE(INOUT2)
                  ELSE
                    OPEN(UNIT=INOUT2,FILE=PATH(1:LPATH)//
     &                '/data/thermo/'//SPECIES(1:LSPECIES),STATUS='OLD')
                    COM2 = '   '
                    DO WHILE (COM2(1:4)/= 'STOP')
                      READ(INOUT2,*) FULLCOM2
                      COM2 = FULLCOM2(1:4)
                      IF (COM2(1:4)=='Diss') THEN
                        READ(INOUT2,*) WR3(IDIS+I-1+IC2-NC-IAI)
                        WR3(IDIS+I-1+IC2-NC-IAI) = 
     &                      WR3(IDIS+I-1+IC2-NC-IAI)*WR1(IUR)
                      ENDIF
                    ENDDO
                    CLOSE(INOUT2)
                  ENDIF
                ELSE
                  OPEN(UNIT=INOUT2,FILE=PATH(1:LPATH)//
     &                '/data/thermo/'//SPECIES(1:LSPECIES),STATUS='OLD')
                  COM2 = '   '
                  DO WHILE (COM2(1:4)/= 'STOP')
                    READ(INOUT2,*) FULLCOM2
                    COM2 = FULLCOM2(1:4)
                    IF (COM2(1:4)=='Diss') THEN
                      READ(INOUT2,*) WR3(IDIS+I-1+IC2-NC-IAI)
                      WR3(IDIS+I-1+IC2-NC-IAI) = 
     &                       WR3(IDIS+I-1+IC2-NC-IAI)*WR1(IUR)
                    ENDIF
                  ENDDO
                  CLOSE(INOUT2)
                ENDIF
             ELSE
                IAI = IAI +1
             ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      CLOSE(INOUT1)

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
1003  FORMAT (A)
   
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
          WRITE(*,*) 'Character not recognized in DTRANSLATE >> ',
     &              WORD(I:I)
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
C----------------------------------------------------------------------- 
      END SUBROUTINE ITRANSLATE
C----------------------------------------------------------------------- 
C-----------------------------------------------------------------------
      SUBROUTINE ARRHENIUSRAD (WR1, LWR1, WR3, LWR3, WI, LWI, Y, YTOL, 
     &                      RHO, TE, HINT, ESCAPE,OMEGARAD)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      INTEGER  LWR1, LWR3, LWI, WI(1:LWI), J, I
      DOUBLE PRECISION WR1(1:LWR1), WR3(1:LWR3), Y(1:NS), YTOL, P, 
     &                 RHO, OMEGARAD
C-----------------------------------------------------------------------
      DOUBLE PRECISION XI(1:NREA), KF(1:NREA), YT(1:NS), PRODR, 
     &         MOLE(1:NS), AJI, HINT(1:NS), TE, DELTAE(1:NREA),
     &         DELTAER, DELTAEP, ESCAPE
      INTEGER IS, IR, IV, ISDIS, NVIBMODE, IC, JR
C-----------------------------------------------------------------------
C     Tolerance on mass fractions
      DO IS = 1, NS
        IF (Y(IS) < YTOL) THEN
          YT(IS) = YTOL
        ELSE
          YT(IS) = Y(IS)
        ENDIF
      ENDDO

C     Molar density [mole /m^3]
      DO IS = 1, NS
        MOLE(IS)   = RHO *YT(IS) /WR1(IMI+IS-1)
      ENDDO

C     Reaction rates
      DO IR = 1, NREA
         DELTAE(IR)  = 0.D0
         IF(WI(ITNEQ+IR-1)== 11) THEN
            KF(IR) =(ESCAPE)*WR3(IRATE+(IR-1)*4+1) 
         ELSE 
            IF(WI(ITNEQ+IR-1)== 12) THEN
              I = IDINT(WR3(IRATE+(IR-1)*4+2))
              J = IDINT(WR3(IRATE+(IR-1)*4+1))
              CALL ProbaEmiSpontMolec(I,J,TE,AJI)
              KF(IR) = (ESCAPE)*AJI
              DELTAE(IR)= 0.D0
C              PRODR = 0.D0
C              DO IS = 1, NS
C                 PRODR = PRODR +WR3(ISTOIR+(IR-1)*NS+IS-1)*MOLE(IS)
C              ENDDO
              DELTAER = 0.D0
              DELTAEP = 0.D0
              DO IS = 1, NS
                DELTAER = DELTAER +WR3(ISTOIR+(IR-1)*NS+IS-1)*HINT(IS)
                DELTAEP = DELTAEP +WR3(ISTOIP+(IR-1)*NS+IS-1)*HINT(IS)
              ENDDO
              DELTAE(IR)  = DELTAER-DELTAEP

            ELSE
              KF(IR)= 0.D0
            ENDIF
         ENDIF
      ENDDO

C     Radiation source terms [J/(m^3 s)]
      OMEGARAD = 0.D0
      DO IR = 1, NREA
        PRODR = 0.D0
        DO IS = 1, NS
          PRODR = PRODR +WR3(ISTOIR+(IR-1)*NS+IS-1)*MOLE(IS)
        ENDDO
        XI(IR) = (KF(IR) *PRODR)
        IF(WI(ITNEQ+IR-1)== 11) THEN
!              IF((WR3(IRATE+(IR-1)*4+2)-WR3(IRATE+(IR-1)*4+3))
!     &            .GT.50000.0) THEN
!                  OMEGARAD = OMEGARAD + 1.D-30
!              ELSE
            OMEGARAD = OMEGARAD +XI(IR)*(WR3(IRATE+(IR-1)*4+2)
     &         -WR3(IRATE+(IR-1)*4+3))*1.43876866D0*WR1(IUR)
!              ENDIF
        ELSE
          OMEGARAD = OMEGARAD +XI(IR)*DELTAE(IR)
        ENDIF
      ENDDO
C-----------------------------------------------------------------------
      END SUBROUTINE ARRHENIUSRAD
C-----------------------------------------------------------------------
        subroutine ProbaEmiSpontMolec(i,j,Tee,Aji)
        implicit none
        integer i, j, ii
        double precision A(0:10), Aji, Te, kref, Tref, som, Tee
        Te = Tee
        IF(Te.le.2000.D0) Te = 2000.D0

c-------Cas de NO beta:
        if(i.eq.21.and.j.eq.23) then
          A(0) = -0.00236888213
          A(1) = 3563.24299
          A(2) = -1240543.35
          A(3) = 7.25178582d+010
          A(4) = -1.65285487d+015
          A(5) = 9.71925285d+018
          A(6) = -2.67424905d+022
          A(7) = 3.55933533d+025
          A(8) = -1.30459694d+028
          A(9) = -1.78911054d+031
          A(10) = 1.52352628d+034
          kref = 5.175d+005
          Tref = 2000.
        endif
c-------Cas de NO gamma:
          if(i.eq.21.and.j.eq.22) then
          A(0) = 0.0188102144
          A(1) = 3334.41698
          A(2) = 47170368.5
          A(3) = -1.1573668d+012
          A(4) = 8.08095109d+015
          A(5) = -2.49649264d+019
          A(6) = 1.96071635d+022
          A(7) = 9.64253718d+025
          A(8) = -3.1398653d+029
          A(9) = 3.79607944d+032
          A(10) = -1.72337581d+035
          kref = 5.104d+006
          Tref = 2000.
        endif

c-------Cas de N2(1er  positif):
        if(i.eq.38.and.j.eq.39) then
          A(0) = -0.00133793817
          A(1) = 1992.23723
          A(2) = 852462.587
          A(3) = 5.85224217d+010
          A(4) = -1.52128718d+015
          A(5) = 1.11985551d+019
          A(6) = -4.33090791d+022
          A(7) = 1.00068707d+026
          A(8) = -1.39757749d+029
          A(9) = 1.09449172d+032
          A(10) = -3.70418948d+034
          kref = 8.107d+004
          Tref = 2000.
        endif

c	Cas de N2(2eme positif):
        if(i.eq.39.and.j.eq.41) then
          A(0) = 2.13357072d-006
          A(1) = -87.0701493
          A(2) = 138956.36
          A(3) = 698266078.
          A(4) = 3.68595835d+011
          A(5) = -1.73223046d+016
          A(6) = 5.66563877d+019
          A(7) = -7.66215715d+022
          A(8) = 2.98157628d+025
          A(9) = 3.04011041d+028
          A(10) = -2.58886484d+031
          kref = 2.725d+007
          Tref = 2000.
        endif

c-------Cas de N2+(1er negatif):
        if(i.eq.43.and.j.eq.45) then
          A(0) = -0.00191545157
          A(1) = -1765.05083
          A(2) = -284122.849
          A(3) = 4.68399019d+011
          A(4) = -6.08541977d+015
          A(5) = 3.88464525d+019
          A(6) = -1.48428612d+023
          A(7) = 3.55807797d+026
          A(8) = -5.25999722d+029
          A(9) = 4.39270858d+032
          A(10) = -1.58739412d+035
          kref = 1.707d+007
          Tref = 2000.
        endif

c-------Calcul de la proba d'emission spontanee:
        som = dlog(kref)+A(0)*dlog(Te/Tref)
        do ii=1,10
         som = som-A(ii)/dfloat(ii)/Tref**ii*((Tref/Te)**ii-1.)
       enddo
       Aji = dexp(som)
      end subroutine ProbaEmiSpontMolec
C-----------------------------------------------------------------------
      SUBROUTINE CNBOLTZMANN (T, TV, TE, CNX, CNB, CNA)
C-----------------------------------------------------------------------
C     This subroutines computes the Boltzmann populations of the 
C     electronic excited states of CN radical in thermal nonequilibrium.
C     Data from K.P. Huberg and G. Herzberg, Molecular spectra and 
C     molecular structure constants of diatomic molecules.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION T, TV, TE, CNX, CNB, CNA
C-----------------------------------------------------------------------
      DOUBLE PRECISION QTOT, QVIBROTX, QVIBROTA, QVIBROTB, QVIBROTD
C-----------------------------------------------------------------------
          QVIBROTX  =  1. /(1. -DEXP(-1.44 *2068.59 /TV)) 
     &                 *T /(1.44 *1.8997)
          QVIBROTA  =  1. /(1. -DEXP(-1.44 *1812.56 /TV)) 
     &                 *T /(1.44 *1.7151)
          QVIBROTB  =  1. /(1. -DEXP(-1.44 *2163.9  /TV)) 
     &                 *T /(1.44 *1.973)
          QVIBROTD  =  1. /(1. -DEXP(-1.44 *1004.7  /TV)) 
     &                 *T /(1.44 *1.162)
          QTOT      =  2. *QVIBROTX
     &                +4. *QVIBROTA *DEXP(-1.44 *9245.28 /TE)
     &                +2. *QVIBROTB *DEXP(-1.44 *25752.  /TE)
     &                +4. *QVIBROTD *DEXP(-1.44 *54486.3 /TE)
      CNX       =  2. *QVIBROTX /QTOT
      CNB       =  2. *QVIBROTB *DEXP(-1.44 *25752.  /TE) /QTOT
      CNA       =  4. *QVIBROTA *DEXP(-1.44 *9245.28 /TE) /QTOT

      END SUBROUTINE CNBOLTZMANN
C-----------------------------------------------------------------------
CC-----------------------------------------------------------------------
C      SUBROUTINE ARRHENIUSVT(WR1, LWR1, WR3, LWR3, WI, LWI, Y, YTOL, P,
C     &                       T, TV, RHO, OMEGA, OMEGACV)
CC-----------------------------------------------------------------------
CC             T O   BE    R E M O V E D   !!!!!!!!!!!!
CC     This subroutine computes the reaction rates in thermal
CC     nonequilibrium based on Arrhenius's formula. The reaction rates
CC     are expressed in [mole m^-3 s^-1]. 
CC-----------------------------------------------------------------------
C      IMPLICIT NONE
C      INCLUDE '../general/memory.cmn'
C      INTEGER LWR1, LWR3, LWI, WI(1:LWI), S, J, IVS, IS1, IS2, K, M,
C     &        IK, IT1, IT2, IS,IR,I,IC  
C      DOUBLE PRECISION WR1(1:LWR1), WR3(1:LWR3), Y(1:NS),P , T,
C     &                 TV(1:NVIB), TE, RHO, OMEGA(1:NS), YTOL,
C     &                 XI(1:NREA), KF(1:NREA), KB(1:NREA), KEQ, KFEQ,
C     &                 GT(1:NS), GTE(1:NS), YT(1:NS), OMEGACV(1:NVIB), 
C     &                 PRODP, PRODR, MOLE(1:NS), LNMOLE(1:NS),
C     &                 RT, RTE, LNTE, SUM, GVAI, GVAI1, GVAI2,
C     &                 LNT, PSI4, PSI3, PSI5, THV(1:NV), TEV(1:NV),
C     &                 REACF(1:NS*NREA), REACB(1:NS*NREA), LNRT, LNRTE,
C     &                 THETA1, THETA2, TVI, TMED, LNTMED
CC-----------------------------------------------------------------------    
CC     Vectors containing the stoichiometric reaction numbers: forward and reward 
CC-----------------------------------------------------------------------
C      IC = 0
C      DO I = 1, NREA
C       DO IS = 1,NS
C         IC = IC + 1
C         REACF(IC)=WR3(ISTOIP+(I-1)*NS+IS-1)
C         REACB(IC)=WR3(ISTOIR+(I-1)*NS+IS-1)
C       ENDDO
C      ENDDO
C      TE = TV(1)
CC-----------------------------------------------------------------------      
CC     Thermal nonequilibrium
CC-----------------------------------------------------------------------
C      RT     = WR1(IUR)*T
C      LNT    = DLOG(T)
C      LNRT   = DLOG(RT)
C      RTE    = WR1(IUR)*TE
C      LNTE   = DLOG(TE)
C      LNRTE  = DLOG(RTE)
CC-----------------------------------------------------------------------      
CC     Tolerance on mass fractions
CC-----------------------------------------------------------------------
C      DO IS = 1, NS
C        IF (Y(IS) < YTOL) THEN
C          YT(IS) = YTOL
C        ELSE
C          YT(IS) = Y(IS)
C        ENDIF
C      ENDDO
C      THV(:) = T 
C      TEV(:) = TE
CC-----------------------------------------------------------------------      
CC     Gibbs's free energy (per unit mole) at p = 1.D0 Pa
CC-----------------------------------------------------------------------
C      CALL GIBBS (WR1, LWR1, WI, LWI, T, T, T, THV, 1.D0, GT)
C      CALL GIBBS (WR1, LWR1, WI, LWI, TE, TE, TE, TEV, 1.D0, GTE)      
CC-----------------------------------------------------------------------      
CC     Molar density [mole /m^3]
CC-----------------------------------------------------------------------
C      DO IS = 1, NS
C        MOLE(IS)   = RHO*YT(IS) /WR1(IMI+IS-1)
C        LNMOLE(IS) = LOG(MOLE(IS))
C      ENDDO
CC-----------------------------------------------------------------------      
CC     Reaction rates
CC-----------------------------------------------------------------------
C      DO IR = 1, NREA
C        SELECT CASE(WI(ITNEQ+IR-1))
CC         KF(SQRT(T*TV)), KB(T)
C          CASE(1)
C          WRITE(*,*)'YOU ARE AN IDIOT! ASK MARCO '
C          STOP
CC         KF(TE), KB(TE)
C          CASE(2)
C            KF(IR) = WR3(IRATE+(IR-1)*4) 
C     &              *DEXP(WR3(IRATE+(IR-1)*4+1) *LNTE
C     &              -WR3(IRATE+(IR-1)*4+2) /TE)
C            KFEQ   = KF(IR)        
C            CALL KEQUILIBRIUM (NS, NREA, GTE, REACB, REACF, 
C     &                         WR1(IMI:IMI+NS-1), RTE, LNRTE, KEQ, IR)
C        KB(IR) = KFEQ /KEQ
CC         KF(T), KB(T)
CC-----------------------------------------------------------------------
CC       ELECTRON-EXCHANGE reactions
CC-----------------------------------------------------------------------  
C          CASE(3)
C            KF(IR) = WR3(IRATE+(IR-1)*4) 
C     &               *DEXP( WR3(IRATE+(IR-1)*4+1)*LNT
C     &               -WR3(IRATE+(IR-1)*4+2) /T)
C
C            KFEQ   = KF(IR)
CC
C            CALL KEQUILIBRIUM (NS, NREA, GT, REACB, REACF,
C     &                         WR1(IMI:IMI+NS-1), RT, LNRT, KEQ, IR)
C           
C            KB(IR) = KFEQ /KEQ
C-----------------------------------------------------------------------
C       ASSOCIATIVE-IONIZATION reactions
C-----------------------------------------------------------------------  
C          CASE(4)
C            KF(IR) = WR3(IRATE+(IR-1)*4) 
C     &               *DEXP(WR3(IRATE+(IR-1)*4+1)*LNT
C     &               -WR3(IRATE+(IR-1)*4+2) /T)
C
C         CALL KEQUILIBRIUM (NS,NREA, GTE, REACB, REACF,
C     &                      WR1(IMI:IMI+NS-1), RTE, LNRTE, KEQ, IR)  
C       
C             SELECT CASE(PHYSICS)
C                 CASE(2)
C                    PSI5 = 1.D0
C                 CASE(3) 
C                    CALL DEVIATIONFACTOR5(WR1,LWR1, WR3, LWR3, WI, LWI,
C     &                                    WR3(IRATE+(IR-1)*4+2), T, TV,
C     &                                    PSI5, WI(ITNEQ+IR-1), IR)
C             END SELECT          
C            
C            KFEQ   = WR3(IRATE+(IR-1)*4) 
C     &              *DEXP( WR3(IRATE+(IR-1)*4+1)*LNTE
C     &              -WR3(IRATE+(IR-1)*4+2)/TE)
C                    
C            KB(IR) = PSI5*KFEQ /KEQ            
CC-----------------------------------------------------------------------
CC       DISSOCIATION-RECOMBINATION reactions
CC-----------------------------------------------------------------------
C          CASE(5)
C            SELECT CASE(PHYSICS)
C                 CASE(2)
C                    TMED = 0.D0
C                    LNTMED = 0.D0
C                    DO I = 1 ,NS
C                         IF (WR3(ISTOIR+(IR-1)*NS+I-1).EQ.1.D0) THEN
C                             S = I    
C                         ENDIF
C                    ENDDO
C                    IC = 0.D0
C                    DO I = 1,NS
C                      IF (WI(IATOMI+I-1).GT.1) THEN
C                          IC = IC + 1
C                            IF (I.EQ.S) THEN
C                                TVI = TV(IC)
C                            ENDIF
C                      ENDIF
C                    ENDDO
C                    TMED = (T*TVI)**0.5D0
C                    LNTMED = DLOG(TMED)
C
C                    KF(IR) =  WR3(IRATE+(IR-1)*4) 
C     &                       *DEXP(WR3(IRATE+(IR-1)*4+1)*LNTMED
C     &                       -WR3(IRATE+(IR-1)*4+2)/TMED)
C                   
C                    KFEQ = KF(IR)
C 
C                 CASE(3) 
C            CALL DEVIATIONFACTOR4(WR1, LWR1, WR3, LWR3, WI, LWI, 
C     &                            WR3(IRATE+(IR-1)*4+2), T, TV, PSI4,
C     &                            WI(ITNEQ+IR-1), IR, TVI)
C 
C                     KF(IR) = PSI4*WR3(IRATE+(IR-1)*4) 
C     &                        *DEXP(WR3(IRATE+(IR-1)*4+1)*LNT
C     &                        -WR3(IRATE+(IR-1)*4+2)/T) 
C     
C                     KFEQ   = WR3(IRATE+(IR-1)*4) 
C     &                        *DEXP(WR3(IRATE+(IR-1)*4+1)*LNT
C     &                        -WR3(IRATE+(IR-1)*4+2)/T)  
C            END SELECT  
C                 
C            CALL KEQUILIBRIUM (NS, NREA, GT, REACB, REACF,
C     &                         WR1(IMI:IMI+NS-1), RT, LNRT, KEQ, IR)
C            KB(IR) = KFEQ/KEQ            
CC-----------------------------------------------------------------------
CC       EXCHANGE reactions
CC-----------------------------------------------------------------------             
C          CASE(6)
C            SELECT CASE(PHYSICS)
C                CASE(2)
C                   PSI5 = 1.D0
C                   TMED = 0.D0
C                   LNTMED = 0.D0
C                   DO I = 1 ,NS
C                        IF (WR3(ISTOIR+(IR-1)*NS+I-1).EQ.1.D0.AND.
C     &                      WI(IATOMI+I-1).GT.1) THEN
C                            S = I     
C                        ENDIF
C                   ENDDO
C                   IC = 0.D0
C                   DO I = 1,NS
C                     IF (WI(IATOMI+I-1).GT.1) THEN
C                         IC = IC + 1
C                           IF (I.EQ.S) THEN
C                               TVI = TV(IC)
C                           ENDIF
C                     ENDIF
C                   ENDDO
C                   TMED = (T*TVI)**0.5D0
C                   LNTMED = DLOG(TMED)
C
C                   KF(IR) =  WR3(IRATE+(IR-1)*4) 
C     &                       *DEXP(WR3(IRATE+(IR-1)*4+1)*LNTMED
C     &                       -WR3(IRATE+(IR-1)*4+2)/TMED)
C                  
C                   KFEQ = KF(IR)
C
C                CASE(3)
C          CALL DEVIATIONFACTOR3(WR1, LWR1, WR3, LWR3, WI, LWI, 
C     &                          WR3(IRATE+(IR-1)*4+2), T, TV, PSI3, 
C     &                          WI(ITNEQ+IR-1), IR, TVI)
C             
C          CALL DEVIATIONFACTOR5(WR1, LWR1, WR3, LWR3, WI, LWI, 
C     &                          WR3(IRATE+(IR-1)*4+2),T , TV, PSI5, 
C     &                          WI(ITNEQ+IR-1), IR)
C
C                   KF(IR) = PSI3*WR3(IRATE+(IR-1)*4) 
C     &                      *DEXP( WR3(IRATE+(IR-1)*4+1)*LNT
C     &                      -WR3(IRATE+(IR-1)*4+2)/T)
C 
C                   KFEQ   = WR3(IRATE+(IR-1)*4

