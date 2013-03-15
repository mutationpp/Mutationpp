C-----------------------------------------------------------------------
      SUBROUTINE FILLARRAY (NMIX, MIX, NS, NSPC, SPC, ARRAY,
     &                      TH, P, RHO, ND, MM, TT, PT, PTB, PTSV, PTIS,
     &                      MIXH1, MIXH2, MIXH3, MIXH4, MIXH5, MIXH6,
     &                      MIXS1, MIXS2, MIXS3, MIXS4, MIXS5, MIXS6,
     &                      MIXCP, MIXFCP, CPE, CPR, CPV, CPINT,
     &                      GAMMAF, GAMMAE, GAMMAI, DRHODP, 
     &                      SOUNDF, SOUNDE, MACH, MIXMCP, 
     &                      ETA, LAMBDATH, LAMBDATE,
     &                      LAMBDAINTE, LAMBDAINTR, LAMBDAINTV, 
     &                      LAMBDAINTH, LAMBDAREA, LAMBDASM, LAMBDAR,
     &                      LAMBDAFI, LAMBDAFR, LAMBDAK, LAMBDATDH, 
     &                      LAMBDATDE, LAMBDATD, 
     &                      LAMBDATOT, ETAI, LAMBDAI, SIGMA, SIGMAI, 
     &                      MFP,
     &                      X, Y, CHIH, CHIE, JDIF, JDIFR, JDIFK, DF, 
     &                      EAMB, EAMBR, EAMBK, H1, H2, H3, H4, H5, H6,
     &                      S1, S2, S3, S4, S5, MI)
C-----------------------------------------------------------------------
C     This subroutine fills in the output array for mutation.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NMIX, MIX(1:NMIX), NS, NSPC, SPC(1:NSPC)
      DOUBLE PRECISION ARRAY(1:NMIX+NS*NSPC), 
     &                 TH, P, RHO, ND, MM, TT, PT, PTB, PTSV, PTIS,
     &                 MIXH1, MIXH2, MIXH3, MIXH4, MIXH5, MIXH6,
     &                 MIXS1, MIXS2, MIXS3, MIXS4, MIXS5, MIXS6,
     &                 MIXCP, MIXFCP, CPE(1:NS), CPR(1:NS), CPV(1:NS),
     &                 CPINT(1:NS), GAMMAF, GAMMAE, GAMMAI, DRHODP, 
     &                 SOUNDF, SOUNDE, MACH, MIXMCP,
     &                 ETA, LAMBDATH, LAMBDATE, LAMBDAFI, LAMBDAFR,
     &                 LAMBDAINTE, LAMBDAINTR, LAMBDAINTV, LAMBDAINTH, 
     &                 LAMBDAREA, LAMBDASM, LAMBDAR, LAMBDAK, SIGMA, 
     &                 SIGMAI, MFP, LAMBDATDH, LAMBDATDE, LAMBDATD, 
     &                 LAMBDATOT, X(1:NS), Y(1:NS), CHIH(1:NS), 
     &                 CHIE(1:NS), JDIF(1:NS), JDIFR(1:NS), JDIFK(1:NS),
     &                 DF(1:NS), EAMB, EAMBR, EAMBK, ETAI(1:NS), 
     &                 LAMBDAI(1:NS), H1(1:NS), H2(1:NS), H3(1:NS), 
     &                 H4(1:NS), H5(1:NS), H6(1:NS), S1(1:NS), S2(1:NS),
     &                 S3(1:NS), S4(1:NS), S5(1:NS), S6(1:NS), MI(1:NS)
C-----------------------------------------------------------------------
      INTEGER IND, I, S 
C-----------------------------------------------------------------------
      I = 1
C     Mixture properties
      DO IND = 1, NMIX
        SELECT CASE (MIX(I))
          CASE(1) 
C           Equilibrium temperature
            ARRAY(I) = TH 
            I = I +1
          CASE(2) 
C           Pressure
            ARRAY(I) = P
            I = I +1
          CASE(3) 
C           Number Density
            ARRAY(I) = ND 
            I = I +1
          CASE(4) 
C           Density
            ARRAY(I) = RHO 
            I = I +1
          CASE(5)
C           Molar masss 
            ARRAY(I) = MM
            I = I +1
          CASE(6)
C           Total enthalpy (per unit mole)
            ARRAY(I) = MIXH1
            I = I +1
          CASE(7)
C           Translational enthalpy (per unit mole)
            ARRAY(I) = MIXH2
            I = I +1
          CASE(8)
C           Electronic enthalpy (per unit mole)
            ARRAY(I) = MIXH3
            I = I +1
          CASE(9)
C           Rotational enthalpy (per unit mole)
            ARRAY(I) = MIXH4
            I = I +1
          CASE(10)
C           Vibrational enthalpy (per unit mole)
            ARRAY(I) = MIXH5
            I = I +1
          CASE(11)
C           Formation enthalpy (per unit mole)
            ARRAY(I) = MIXH6
            I = I +1
          CASE(12)
C           Total enthalpy (per unit mass)
            ARRAY(I) = MIXH1 /MM
            I = I +1
          CASE(13)
C           Translational enthalpy (per unit mass)
            ARRAY(I) = MIXH2 /MM
            I = I +1
          CASE(14)
C           Electronic enthalpy (per unit mass)
            ARRAY(I) = MIXH3 /MM
            I = I +1
          CASE(15)
C           Rotational enthalpy (per unit mass)
            ARRAY(I) = MIXH4 /MM
            I = I +1
          CASE(16)
C           Vibrational enthalpy (per unit mass)
            ARRAY(I) = MIXH5 /MM
            I = I +1
          CASE(17)
C           Formation enthalpy (per unit mass)
            ARRAY(I) = MIXH6 /MM
            I = I +1
          CASE(18)
C           Total entropy (per unit mole)
            ARRAY(I) = MIXS1
            I = I +1
          CASE(19)
C           Translational entropy (per unit mole)
            ARRAY(I) = MIXS2
            I = I +1
          CASE(20)
C           Electronic entropy (per unit mole)
            ARRAY(I) = MIXS3
            I = I +1
          CASE(21)
C           Rotational entropy (per unit mole)
            ARRAY(I) = MIXS4
            I = I +1
          CASE(22)
C           Vibrational entropy (per unit mole)
            ARRAY(I) = MIXS5
            I = I +1
          CASE(23)
C           Entropy of mixing (per unit mole)
            ARRAY(I) = MIXS6
            I = I +1
          CASE(24)
C           Total entropy (per unit mass)
            ARRAY(I) = MIXS1 /MM
            I = I +1
          CASE(25)
C           Translational entropy (per unit mass)
            ARRAY(I) = MIXS2 /MM
            I = I +1
          CASE(26)
C           Electronic entropy (per unit mass)
            ARRAY(I) = MIXS3 /MM
            I = I +1
          CASE(27)
C           Rotational entropy (per unit mass)
            ARRAY(I) = MIXS4 /MM
            I = I +1
          CASE(28)
C           Vibrational entropy (per unit mass)
            ARRAY(I) = MIXS5 /MM
            I = I +1
          CASE(29)
C           Entropy of mixing (per unit mass)
            ARRAY(I) = MIXS6 /MM
            I = I +1
          CASE(30)
C           Specific heat (per unit mole)
            ARRAY(I) = MIXCP
            I = I +1
          CASE(31)
C           Frozen specific heat (per unit mole)
            ARRAY(I) = MIXFCP
            I = I +1
          CASE(32)
C           Specific heat (per unit mass)
            ARRAY(I) = MIXMCP
            I = I +1
          CASE(33)
C           Frozen specific heat (per unit mass)
            ARRAY(I) = MIXFCP /MM
            I = I +1
          CASE(34)
C           Frozen specific heat ratio
            ARRAY(I) = GAMMAF
            I = I +1
          CASE(35)
C           Equilibrium specific heat ratio (per unit mass)
            ARRAY(I) = GAMMAE
            I = I +1
          CASE(36)
C           Isentropic exponent (per unit mass)
            ARRAY(I) = GAMMAI
            I = I +1
          CASE(37)
C           Frozen sound speed
            ARRAY(I) = SOUNDF
            I = I +1
          CASE(38)
C           Equilibrium sound speed
            ARRAY(I) = SOUNDE
            I = I +1
          CASE(39)
C           Equilibrium Mach number
            ARRAY(I) = MACH
            I = I +1
          CASE(40)
C           Total temperature
            ARRAY(I) = TT
            I = I +1
          CASE(41)
C           Total pressure
            ARRAY(I) = PT
            I = I +1
          CASE(42)
C           Total pressure (Bernoulli)
            ARRAY(I) = PTB
            I = I +1
          CASE(43)
C           Total pressure (Saint-Venant, specific heat ratio)
            ARRAY(I) = PTSV
            I = I +1
          CASE(44)
C           Total pressure (Saint-Venant, isentropic exponent)
            ARRAY(I) = PTIS
            I = I +1
          CASE(100)
C           Viscosity
            ARRAY(I) = ETA
            I = I +1
          CASE(101)
C           Translational thermal conductivity (heavy particles)
            ARRAY(I) = LAMBDATH
            I = I +1
          CASE(102)
C           Translational thermal conductivity (free electrons)
            ARRAY(I) = LAMBDATE
            I = I +1
          CASE(103)
C           Electronic internal thermal conductivity
            ARRAY(I) = LAMBDAINTE
            I = I +1
          CASE(104)
C           Rotational internal thermal conductivity
            ARRAY(I) = LAMBDAINTR
            I = I +1
          CASE(105)
C           Vibrational internal thermal conductivity
            ARRAY(I) = LAMBDAINTV
            I = I +1
          CASE(106)
C           Total internal thermal conductivity
            ARRAY(I) = LAMBDAINTH
            I = I +1
          CASE(107)
C           Reactive thermal conductivity (SM Magin)
            ARRAY(I) = LAMBDASM
            I = I +1
          CASE(108)
C           Reactive thermal conductivity (SM Ramshaw)
            ARRAY(I) = LAMBDAR
            I = I +1
          CASE(109)
C           Reactive thermal conductivity (SM Kolesnikov)
            ARRAY(I) = LAMBDAK
            I = I +1
          CASE(110)
C           Reactive thermal conductivity (Fick)
            ARRAY(I) = LAMBDAFI
            I = I +1
          CASE(111)
C           Reactive thermal conductivity (Fick-Ramshaw)
            ARRAY(I) = LAMBDAFR
            I = I +1
          CASE(112)
C           Reactive thermal conductivity (Butler and Brokaw)
            ARRAY(I) = LAMBDAREA
            I = I +1
          CASE(113)
C           Heavy thermal diffusion thermal conductivity
            ARRAY(I) = LAMBDATDH
            I = I +1
          CASE(114)
C           Electron thermal diffusion thermal conductivity
            ARRAY(I) = LAMBDATDE
            I = I +1
          CASE(115)
C           Total thermal diffusion thermal conductivity 
            ARRAY(I) = LAMBDATD
            I = I +1
          CASE(116)
C           Total thermal conductivity 
            ARRAY(I) = LAMBDATOT
            I = I +1
          CASE(117)
C           Ambipolar electric field (SM Magin)
            ARRAY(I) = EAMB
            I = I +1
          CASE(118)
C           Ambipolar electric field (SM Ramshaw) 
            ARRAY(I) = EAMBR
            I = I +1
          CASE(119)
C           Ambipolar electric field (SM Kolesnikov)
            ARRAY(I) = EAMBK
            I = I +1
          CASE(120)
C           Electrical conductivity (electrons)
            ARRAY(I) = SIGMA
            I = I +1
          CASE(121)
C           Electrical conductivity (ions)
            ARRAY(I) = SIGMAI
            I = I +1
          CASE(122)
C           Mean free path 
            ARRAY(I) = MFP
          CASE DEFAULT
            WRITE(*,*) 'Mixture output entry',MIX(IND),' wrong.'
        END SELECT
      ENDDO

C     Species properties
      DO IND = 1, NSPC
        SELECT CASE (SPC(IND))
          CASE(1) 
C           Equilibrium composition (mole fraction)
            DO S = 1, NS
              ARRAY(I) = X(S)
              I = I +1
            ENDDO
          CASE(2) 
C           Equilibrium composition (mass fraction)
            DO S = 1, NS
              ARRAY(I) = Y(S)
              I = I +1
            ENDDO
          CASE(3) 
C           Equilibrium composition (number density)
            DO S = 1, NS
              ARRAY(I) = X(S) *ND
              I = I +1
            ENDDO
          CASE(4)  
C           Total enthalpy (per unit mole)
            DO S = 1, NS
              ARRAY(I) = H1(S)
              I = I +1
            ENDDO
          CASE(5)  
C           Translational enthalpy (per unit mole)
            DO S = 1, NS
              ARRAY(I) = H2(S)
              I = I +1
            ENDDO
          CASE(6)  
C           Electronic enthalpy (per unit mole)
            DO S = 1, NS
              ARRAY(I) = H3(S)
              I = I +1
            ENDDO
          CASE(7)  
C           Rotational enthalpy (per unit mole)
            DO S = 1, NS
              ARRAY(I) = H4(S)
              I = I +1
            ENDDO
          CASE(8)  
C           Vibrational enthalpy (per unit mole)
            DO S = 1, NS
              ARRAY(I) = H5(S)
              I = I +1
            ENDDO
          CASE(9)  
C           Formation enthalpy (per unit mole)
            DO S = 1, NS
              ARRAY(I) = H6(S)
              I = I +1
            ENDDO
          CASE(10)  
C           Total enthalpy (per unit mass)
            DO S = 1, NS
              ARRAY(I) = H1(S) /MI(S)
              I = I +1
            ENDDO
          CASE(11)  
C           Translational enthalpy (per unit mass)
            DO S = 1, NS
              ARRAY(I) = H2(S) /MI(S)
              I = I +1
            ENDDO
          CASE(12)  
C           Electronic enthalpy (per unit mass)
            DO S = 1, NS
              ARRAY(I) = H3(S) /MI(S)
              I = I +1
            ENDDO
          CASE(13)  
C           Rotational enthalpy (per unit mass)
            DO S = 1, NS
              ARRAY(I) = H4(S) /MI(S)
              I = I +1
            ENDDO
          CASE(14)  
C           Vibrational enthalpy (per unit mass)
            DO S = 1, NS
              ARRAY(I) = H5(S) /MI(S)
              I = I +1
            ENDDO
          CASE(15)  
C           Formation enthalpy (per unit mass)
            DO S = 1, NS
              ARRAY(I) = H6(S) /MI(S)
              I = I +1
            ENDDO
          CASE(16)  
C           Total entropy (per unit mole)
            DO S = 1, NS
              ARRAY(I) = S1(S)
              I = I +1
            ENDDO
          CASE(17)  
C           Translational entropy (per unit mole)
            DO S = 1, NS
              ARRAY(I) = S2(S)
              I = I +1
            ENDDO
          CASE(18)  
C           Electronic entropy (per unit mole)
            DO S = 1, NS
              ARRAY(I) = S3(S)
              I = I +1
            ENDDO
          CASE(19)  
C           Rotational entropy (per unit mole)
            DO S = 1, NS
              ARRAY(I) = S4(S)
              I = I +1
            ENDDO
          CASE(20)  
C           Vibrational entropy (per unit mole)
            DO S = 1, NS
              ARRAY(I) = S5(S)
              I = I +1
            ENDDO
          CASE(21)  
C           Total entropy (per unit mass)
            DO S = 1, NS
              ARRAY(I) = S1(S) /MI(S)
              I = I +1
            ENDDO
          CASE(22)  
C           Translational entropy (per unit mass)
            DO S = 1, NS
              ARRAY(I) = S2(S) /MI(S)
              I = I +1
            ENDDO
          CASE(23)  
C           Electronic entropy (per unit mass)
            DO S = 1, NS
              ARRAY(I) = S3(S) /MI(S)
              I = I +1
            ENDDO
          CASE(24)  
C           Rotational entropy (per unit mass)
            DO S = 1, NS
              ARRAY(I) = S4(S) /MI(S)
              I = I +1
            ENDDO
          CASE(25)  
C           Vibrational entropy (per unit mass)
            DO S = 1, NS
              ARRAY(I) = S5(S) /MI(S)
              I = I +1
            ENDDO
          CASE(26)  
C           Electronic specific heat (per unit mole)
            DO S = 1, NS
              ARRAY(I) = CPE(S)
              I = I +1
            ENDDO
          CASE(27)  
C           Rotational specific heat (per unit mole)
            DO S = 1, NS
              ARRAY(I) = CPR(S)
              I = I +1
            ENDDO
          CASE(28)  
C           Vibrational specific heat (per unit mole)
            DO S = 1, NS
              ARRAY(I) = CPV(S)
              I = I +1
            ENDDO
          CASE(29)  
C           Internal specific heat (per unit mole)
            DO S = 1, NS
              ARRAY(I) = CPINT(S)
              I = I +1
            ENDDO
          CASE(30)  
C           Electronic specific heat (per unit mass)
            DO S = 1, NS
              ARRAY(I) = CPE(S) /MI(S)
              I = I +1
            ENDDO
          CASE(31)  
C           Rotational specific heat (per unit mass)
            DO S = 1, NS
              ARRAY(I) = CPR(S) /MI(S)
              I = I +1
            ENDDO
          CASE(32)  
C           Vibrational specific heat (per unit mass)
            DO S = 1, NS
              ARRAY(I) = CPV(S) /MI(S)
              I = I +1
            ENDDO
          CASE(33)  
C           Internal specific heat (per unit mass)
            DO S = 1, NS
              ARRAY(I) = CPINT(S) /MI(S)
              I = I +1
            ENDDO
          CASE(100)  
C           Driving forces 
            DO S = 1, NS
              ARRAY(I) = DF(S)
              I = I +1
            ENDDO
          CASE(101)  
C           Pure viscosity 
            DO S = 1, NS
              ARRAY(I) = ETAI(S)
              I = I +1
            ENDDO
          CASE(102)  
C           Pure thermal conductivity 
            DO S = 1, NS
              ARRAY(I) = LAMBDAI(S)
              I = I +1
            ENDDO
          CASE(103) 
C           Thermal diffusion ratio (heavy particles)
            DO S = 1, NS
              ARRAY(I) = CHIH(S)
              I = I +1
            ENDDO
          CASE(104)  
C           Thermal diffusion ratio (free electrons)
            DO S = 1, NS
              ARRAY(I) = CHIE(S)
              I = I +1
            ENDDO
          CASE(105)  
C           Mass diffusion flux (SM Magin)
            DO S = 1, NS
              ARRAY(I) = JDIF(S)
              I = I +1
            ENDDO
          CASE(106)  
C           Mass diffusion flux (SM Ramshaw)
            DO S = 1, NS
              ARRAY(I) = JDIFR(S)
              I = I +1
            ENDDO
          CASE(107)  
C           Mass diffusion flux (SM Kolesnikov)
            DO S = 1, NS
              ARRAY(I) = JDIFK(S)
              I = I +1
            ENDDO
          CASE DEFAULT
            WRITE(*,*) 'Species output entry',SPC(IND),' wrong.'
        END SELECT
      ENDDO

      END SUBROUTINE FILLARRAY 
C-----------------------------------------------------------------------
