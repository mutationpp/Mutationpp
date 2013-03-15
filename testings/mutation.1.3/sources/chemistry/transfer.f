C-----------------------------------------------------------------------
      SUBROUTINE OMEGAVTRANSFER (LWR1, WR1, LWR4, WR4, LWI, WI, P, TVEC,
     &                          Y, RHO, H5T, H5, OMEGAVT, OMEGAVE,
     &                          OMEGAVV)
C-----------------------------------------------------------------------
      IMPLICIT NONE 
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      INTEGER LWR1, LWR4, LWI, WI(1:LWI)
      DOUBLE PRECISION WR1(1:LWR1), WR4(1:LWR4), P, TVEC(1:NVIB+2), 
     &                 Y(1:NS), RHO, H5(1:NS), H5T(1:NS), 
     &                 OMEGAVT(1:NVIB+1), OMEGAVE(1:NVIB+1), 
     &                 OMEGAVV(1:NVIB+1)
C-----------------------------------------------------------------------
      INTEGER I, J, IVT, IC, NVT
      DOUBLE PRECISION VT, VE, VV, TAUVT(1:NS*NS)
C-----------------------------------------------------------------------
C     Millikan-White-Park relaxation time for VT transfer
C     ---------------------------------------------------
      NVT = 0
      DO I = 1, NVIB
        DO J = 1, WI(INVIBSPEI+I-1)
          NVT = NVT +1
        ENDDO
      ENDDO

      CALL COMPUTETAUVT (NE, NS, NVT, NVIB, WI(IVIBSPEI),  P, TVEC(1),
     &                   WR1(IMI), WR1(IUPI),   WR1(IUNA), WR1(IUKB),
     &                   WR4(IMWAIJ),WR4(IMWBIJ), TAUVT)   

C     Energy transfer
C     ---------------
      IC = 0
      DO I = 1, NVIB
        OMEGAVT(I) = 0.D0; OMEGAVV(I) = 0.D0; OMEGAVE(I) = 0.D0
        DO J = 1, WI(INVIBSPEI+I-1)
          IC = IC +1
          IVT = WI(IVIBSPEI+IC-1)
          CALL VTTRANSFER (NS, NE, Y, RHO, H5T(IVT), H5(IVT), WR1(IMI),
     &                     IVT, TAUVT(1+(NS-NE)*(IC-1)), VT) 
          IF (NVIB > 1) THEN
            VV = 0.D0
            VE = 0.D0
C          CALL VVTRANSFER 
C          CALL VETRANSFER  
          ELSE
            VV = 0.D0
            VE = 0.D0
          ENDIF
          OMEGAVT(I) = OMEGAVT(I) +VT
          OMEGAVV(I) = OMEGAVV(I) +VV
          OMEGAVE(I) = OMEGAVE(I) +VE
        ENDDO
      ENDDO 

      END SUBROUTINE OMEGAVTRANSFER
C-----------------------------------------------------------------------
C----------------------------------------------------------------------- 
      SUBROUTINE COMPUTETAUVT (NE, NS, NVT, NVIB, VIBSPE, P, T, MASS,
     &                         PI, NA, KB, MWA, MWB, TAUVT)  
C-----------------------------------------------------------------------
C     This subroutine computes the vibrational-translational relaxation 
C     times based on Millikan-White's formula for the relaxation time 
C     and including Park's correction for high temperatures found in
C     C. Park, Review of chemical kinetics problems of future NASA
C     missions, I: Earth entries, Journal of Thermophysics and Heat 
C     Transfer, 1993, 7(3):385.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NE, NS, NVT, NVIB, VIBSPE(1:NS)
      DOUBLE PRECISION  P, T, MASS(1:NS), PI, NA, KB, 
     &                  MWA(1:NVT*(NS-NE)), MWB(1:NVT*(NS-NE)),
     &                  TAUVT(1:NS*NS)
C-----------------------------------------------------------------------
      INTEGER I, IS, IC
      DOUBLE PRECISION MU, MILLIKAN, PARK, SIGMA
C-----------------------------------------------------------------------
C     Limiting cross sections for Park's correction [m^2]
C     ---------------------------------------------------
      SIGMA = 3.D-21 *(50000.D0 /T)**2

C     Relaxation time
C     ---------------
      IC= 0
      DO I = 1, NVT
        DO IS = NE+1, NS
          IC        = IC +1
          MU   = MASS(IS) *MASS(VIBSPE(I)) /( MASS(IS) +MASS(VIBSPE(I)))
          MILLIKAN  = DEXP(MWA(IC) *(T**(-0.3333333333) -MWB(IC)) 
     &                -18.42D0) *101325.D0 /P
          PARK      = DSQRT(PI *MU *KB *T/(8.D0 *NA)) /(SIGMA *P)
          TAUVT(IC) = MILLIKAN +PARK
        ENDDO
      ENDDO

      END SUBROUTINE COMPUTETAUVT
C-----------------------------------------------------------------------
C----------------------------------------------------------------------- 
      SUBROUTINE VTTRANSFER (NS, NE, Y, RHO, ET, ETV, MASS, IVT, TAUVT, 
     &                       OMEGAVT)  
C-----------------------------------------------------------------------
C     This subroutines computes the vibrational-translational energy
C     transfer for the species IVT based on Landau-Teller model. The 
C     kinetic theory predicts a frequency average for the relaxation 
C     time, see for instance:
C     R.N. Schwarz, Z.I. Slawsky, K.F. Herzfeld, Calculation of 
C     vibrational relaxation times in gases, Journal of Chemical 
C     Physics, 1952, 20:1591. 
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NS, NE, IVT
      DOUBLE PRECISION  Y(1:NS), RHO, ET, ETV, MASS(1:NS),
     &                  TAUVT(1:NS), OMEGAVT
C-----------------------------------------------------------------------
      INTEGER IS
      DOUBLE PRECISION TAU, SUM1, SUM2
C-----------------------------------------------------------------------
C     Frequency average over heavy particles
C     --------------------------------------
      SUM1 = 0.D0; SUM2 = 0.D0
      DO IS = NE+1, NS
        SUM1 = SUM1 +Y(IS) /MASS(IS)
        SUM2 = SUM2 +Y(IS) /(MASS(IS) *TAUVT(IS-NE))
      ENDDO
      TAU = SUM1 /SUM2 

C     Landau-Teller's vibrational-translational relaxation model
C     ----------------------------------------------------------
      OMEGAVT =  RHO *Y(IVT) *(ET -ETV) /TAU

      END SUBROUTINE VTTRANSFER
C-----------------------------------------------------------------------
C----------------------------------------------------------------------- 
      SUBROUTINE VVTRANSFER (NS, I, Y, RHO, ET, ETV, MASS, NA, SIGMA,
     &                       R, PI, TH, OMEGAVV) 
C-----------------------------------------------------------------------
C     This subroutines computes the vibrational-vibrational energy
C     transfer for the species IVT based on Landau-Teller model with 
C     Candler's formula.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NS, I
      DOUBLE PRECISION  Y(1:NS), RHO, ET(1:NS), ETV(1:NS), MASS(1:NS),
     &                  NA, SIGMA, R, PI, TH, OMEGAVV
C-----------------------------------------------------------------------
      INTEGER J
      DOUBLE PRECISION SUM1, PIJ, SQMU
C-----------------------------------------------------------------------
C     Exchange probability
C     --------------------
      PIJ = 1.D-2

      SUM1 = 0.D0
      DO J = 1, NS
          SQMU = DSQRT(MASS(I) *MASS(J) /(MASS(I) +MASS(J)))
          SUM1 = SUM1 +SIGMA *NA *Y(J) /(MASS(J) *SQMU) 
     &           *(ET(I) /ET(J) *ETV(J) -ETV(I))
      ENDDO

      OMEGAVV =  RHO**3 *Y(I) *Y(I) /MASS(I) *DSQRT(8.D0 *R *TH /PI) 
     &           *PIJ *SUM1

      END SUBROUTINE VVTRANSFER
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE INITIALIZETRANSFER (PATH, LWR1, WR1, LWR4, WR4, LWI,
     &                               WI, LWC, WC)
C-----------------------------------------------------------------------
C     This subroutine initializes temperature and pressure independent 
C     variables (Millikan-White constants) stored in the work arrays
C     WR4, WI and WC. For the VT transfer, if no values are specified 
C     in the data/transfer/VT, Millikan-White's formula is used.
C-----------------------------------------------------------------------
      IMPLICIT NONE 
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      INTEGER LWR1, LWR4, LWI, WI(1:LWI), LWC
      DOUBLE PRECISION WR1(1:LWR1), WR4(1:LWR4)
      CHARACTER(100) PATH
      CHARACTER WC(1:LWC)
C-----------------------------------------------------------------------
      INTEGER I, J, IC1, IC2, IC3, IS, K, IVT, LPATH, LCHAR, LSPECIES1, 
     &        LSPECIES2, LTOT, IREAD
      DOUBLE PRECISION MASS(1:NS), SQMU, THETA(1:NV), A, B
      CHARACTER(10) SPECIES1, SPECIES2
      CHARACTER(21) NAME
      CHARACTER(80) FULLCOM1
C-----------------------------------------------------------------------
      LPATH = LCHAR (PATH)
      DO IS = 1, NS
        MASS(IS) = WR1(IMI+IS-1)
      ENDDO 
      IC1 = 0
      DO I = 1, NVIB
        DO J = 1, WI(INVIBSPEI+I-1)
          IC1 = IC1 +1
          IC2 = 0
          DO K = 1, WI(IVIBSPEI+IC1-1) -1
            IC2 = IC2 +WI(IVIBI+K-1)
          ENDDO
          THETA(IC1) = WR1(ITVIK+IC2)
        ENDDO
      ENDDO

C     VT transfer parameters
      OPEN(UNIT=INOUT1,FILE=PATH(1:LPATH)//'/data/transfer/VT',
     &     STATUS='OLD')
        IC1 = 0; IC2 = 0
        DO I = 1, NVIB
          DO J = 1, WI(INVIBSPEI+I-1)
            IC1 = IC1 +1
            IVT = WI(IVIBSPEI+IC1-1) 
            IC3 = 0
            DO K = 1, IVT
              LSPECIES1 = WI(ILNAMEI+K-1)
              IC3 = IC3 +LSPECIES1
            ENDDO
            IC3 = IC3 -LSPECIES1
            DO K = 1, LSPECIES1
              SPECIES1(K:K) = WC(INAMEI+IC3+K-1)
            ENDDO
            IC3 = 0; IF (NE /= 0) IC3 = WI(ILNAMEI)
            DO IS = NE+1, NS
              IC2  = IC2 +1
              LSPECIES2 = WI(ILNAMEI+IS-1)
              DO K = 1, LSPECIES2
                SPECIES2(K:K) = WC(INAMEI+IC3+K-1)
              ENDDO
              IC3 = IC3 +LSPECIES2
              LTOT      = LSPECIES1 +LSPECIES2 +3
              NAME(1:LTOT) = '.'//SPECIES1(1:LSPECIES1)//'-'//SPECIES2
     &                       (1:LSPECIES2)//'.'
              FULLCOM1 = '                     '
              REWIND(INOUT1) ; IREAD = 0
              DO WHILE (FULLCOM1(1:4)/= 'STOP')
                READ(INOUT1,*) FULLCOM1
                IF (FULLCOM1(1:LTOT) == NAME(1:LTOT)) THEN
                  IREAD = IREAD +1
                  READ(INOUT1,*) WR4(IMWAIJ+IC2-1), WR4(IMWBIJ+IC2-1)
                ENDIF
              ENDDO
C     Default: Millikan-White
              IF  (IREAD == 0) THEN
                WRITE(*,*) '-> Default VT model for ',
     &                     NAME(2:LTOT-1), ' pair'
                SQMU = DSQRT(1.D3 *MASS(IS) *MASS(IVT) 
     &                 /( MASS(IS) +MASS(IVT)))
                WR4(IMWAIJ+IC2-1) = 1.16D-3 *SQMU 
     &                              *THETA(IC1)**1.333333333
                WR4(IMWBIJ+IC2-1) = 0.015D0 *DSQRT(SQMU)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        CLOSE(INOUT1)

      END SUBROUTINE INITIALIZETRANSFER
C-----------------------------------------------------------------------
