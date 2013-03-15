C-----------------------------------------------------------------------
      SUBROUTINE INITIALIZE (PATH, MIXTURE, REACTION, TRANSFER, WR1,
     &                       LWR1, WR3, LWR3, WR4, LWR4, WI, LWI, WC,
     &                       LWC, IMOD)
C-----------------------------------------------------------------------
C     This subroutine initializes temperature and pressure independent 
C     variables (mass, formation enthalpy, charge, atomicity, linearity, 
C     symmetry, spectroscopic data of the particle and stoichiometric 
C     matrix, gas and wall reactions,...) stored in the work arrays WR1,
C     WR3, WR4, WI and WC. This subroutine must be called by the user
C     only once at the beginning of the program.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      INTEGER LWR1, LWR3, LWR4, LWI, LWC, WI(1:LWI), TLSPECIES
      DOUBLE PRECISION  WR1(1:LWR1), WR3(1:LWR3), WR4(1:LWR4)
      CHARACTER(10) MIXTURE, REACTION, TRANSFER
      CHARACTER(100) PATH 
      CHARACTER WC(1:LWC)
C-----------------------------------------------------------------------
      INTEGER SPECIESNUCLEI(1:NC), ATOMSUM, DUM1
      DOUBLE PRECISION MASS(1:NS), SQMU, THETA(1:NV), NORMVEC, DUM2
      CHARACTER(10) SPECIES, SPECIES1, SPECIES2, NUCLEI, WORD
      CHARACTER(21) NAME1, NAME2 
      CHARACTER(80) COM1, COM2, FULLCOM1, FULLCOM2, LINE 
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      LMIXTURE  = LCHAR (MIXTURE)
      LTRANSFER = LCHAR (TRANSFER)
      LPATH = LCHAR (PATH)
      DO I = 1, LWR1
        WR1(I) = 0.D0
      ENDDO
      DO I = 1, LWI
        WI(I) = 0
      ENDDO
      DO I = 1, LWC
        WC(I) = '-'
      ENDDO

C     Universal constants
      WR1(IUKB) = 1.380658D-23   
      WR1(IUNA) = 6.0221367D23  
      WR1(IUR)  = WR1(IUKB) *WR1(IUNA)
      WR1(IUE)  = 1.602191D-19 
      WR1(IUPI) = 3.14159265D0
      WR1(IUE0) = 8.854188D-12 
      WR1(IUH)  = 6.626075D-34
      WR1(IUQT) = 1.5D0*DLOG(2.D0*WR1(IUPI)/(WR1(IUNA)*WR1(IUH)**2))
     &            +2.5D0*DLOG(WR1(IUKB))

C     Mixture file
      ICR = 0
      OPEN(UNIT=INOUT1,FILE=PATH(1:LPATH)//'/data/mixture/'
     &     //MIXTURE(1:LMIXTURE)//'.mix',STATUS='OLD')
      COM1 = '   '
      IC1 = 0; II = 0
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
                    READ(INOUT2,*) 
                    READ(INOUT2,*) 
                    CALL BLANK(WORD)
                    READ(INOUT2,*) DUM1, DUM2, WORD   
                    LSPECIES = LCHAR(WORD)
                    CALL BLANK(SPECIES)
                    DO JJ = 1, LSPECIES
                      SPECIES(JJ:JJ) = WORD(JJ:JJ)
                    ENDDO                    
                  ENDIF
                ENDDO
                CLOSE(INOUT2)
              ENDIF
            ENDIF
            II = II +1
            WI(ILNAMEI+II-1) = LSPECIES 
            DO JJ = 1, LSPECIES
              WC(INAMEI+IC1+JJ-1) = SPECIES(JJ:JJ)
            ENDDO
            IC1 = IC1 + LSPECIES
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
          DO I = 1, NS-NSTAR
            CALL BLANK(SPECIES)
            READ(INOUT1,1003) SPECIES 
            LSPECIES  = LCHAR (SPECIES)
C 12
            IF (LSPECIES>5) THEN
              J = LSPECIES-4
              IF (SPECIES(J:LSPECIES)=='-star') THEN 
                ICR = 1
C     Excited specific species
                OPEN(UNIT=INOUT2,FILE=PATH(1:LPATH)//'/data/thermo/'//
     &               SPECIES(1:J-1),STATUS='OLD')
                COM2 = '   '
                DO WHILE (COM2(1:4)/= 'STOP')
                  READ(INOUT2,*) FULLCOM2
                  COM2 = FULLCOM2(1:4)
                  IF (COM2(1:4)=='Elec') THEN
                    READ(INOUT2,*) NCR
                    READ(INOUT2,*) 
                    IF (I<=NC) THEN
                      I1 = 2
                      READ(INOUT2,*) 
                    ENDIF
                    DO K = I1, NCR
                      CALL BLANK(WORD)
                      READ(INOUT2,*) DUM1, DUM2, WORD   
                      LSPECIES = LCHAR(WORD)
                      CALL BLANK(SPECIES)
                      DO JJ = 1, LSPECIES
                        SPECIES(JJ:JJ) = WORD(JJ:JJ)
                      ENDDO
                      II = II +1
                      WI(ILNAMEI+II-1) = LSPECIES 
                      DO JJ = 1, LSPECIES
                       WC(INAMEI+IC1+JJ-1) = SPECIES(JJ:JJ)
                      ENDDO
                      IC1 = IC1 + LSPECIES
                    ENDDO                    
                  ENDIF
                ENDDO
                CLOSE(INOUT2)
              ELSE
                IF (I>NC) THEN
                  ICR = 0
                ENDIF
              ENDIF
C 12
            ELSE
              IF (I>NC) THEN
                ICR = 0
              ENDIF
            ENDIF
C     Non electronic specific species
C 12
            IF ((ICR==0).AND.(I>NC)) THEN
              II = II +1
              WI(ILNAMEI+II-1) = LSPECIES 
              DO JJ = 1, LSPECIES
                WC(INAMEI+IC1+JJ-1) = SPECIES(JJ:JJ)
              ENDDO
              IC1 = IC1 + LSPECIES
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      CLOSE(INOUT1)

C     Stoichiometric matrix
      DO I = 1, NS
        DO J = 1, NS
          IF (I==J) THEN
            WI(INUIK+(I-1)*NS+J-1) = 1
          ELSE
            WI(INUIK+(I-1)*NS+J-1) = 0
          ENDIF
        ENDDO
      ENDDO
C     Nuclei
      IC1 = 0
      DO I = 1, NC
        LSPECIES = WI(ILNAMEI+I-1)
        NUCLEI(I:I) = WC(INAMEI+IC1)
        IC1 = IC1 +LSPECIES
      ENDDO
C     Charge
      DO I = 1, NC
        WI(IQI+I-1)    =  0
        WI(IATOMI+I-1) =  1
      ENDDO
      IF (NE/=0) THEN
        WI(IQI)        = -1
        WI(IATOMI)     =  0
      ENDIF
C     Other species
      DO I = NC+1, NS
        LSPECIES = WI(ILNAMEI+I-1)
        CALL BLANK (SPECIES)
        DO J = 1, LSPECIES
          SPECIES(J:J) = WC(INAMEI+IC1+J-1)
        ENDDO
        CALL CUTMOLECULE(SPECIES, LSPECIES, NUCLEI, SPECIESNUCLEI, NC)
        IC1 = IC1 +LSPECIES  
        IF (NE==1) THEN
          J = 1
          WI(INUIK+(J-1)*NS+I-1)= -SPECIESNUCLEI(J)
          WI(INUIK+(I-1)*NS+J-1)= +SPECIESNUCLEI(J)
        ENDIF
        ATOMSUM = 0 
        DO J = NE+1, NC 
          WI(INUIK+(J-1)*NS+I-1)= +SPECIESNUCLEI(J)
          WI(INUIK+(I-1)*NS+J-1)= -SPECIESNUCLEI(J)
          ATOMSUM = ATOMSUM +SPECIESNUCLEI(J)
        ENDDO
        WI(IQI+I-1) = SPECIESNUCLEI(1)
        WI(IATOMI+I-1) = ATOMSUM
      ENDDO
C     Species not counted twice
      DO I = NC+1, NS
        DO J = NC+1, NS
          IF (I/=J) THEN
            NORMVEC = 0.D0
            DO K = 1, NC
              DIFF = WI(INUIK+(I-1)*NS+K-1)-WI(INUIK+(J-1)*NS+K-1)
              NORMVEC = NORMVEC +DIFF*DIFF
            ENDDO
            IF (NORMVEC<1.D-6) THEN
              L1 = WI(ILNAMEI+I-1)
              L2 = WI(ILNAMEI+J-1)
              IF (L1==L2) THEN
                IC1 = 0
                DO L = 1, I-1
                  LSPECIES = WI(ILNAMEI+L-1)
                  IC1 = IC1 +LSPECIES 
                ENDDO
                IC2 = 0
                DO L = 1, J-1
                  LSPECIES = WI(ILNAMEI+L-1)
                  IC2 = IC2 +LSPECIES 
                ENDDO
                IC1 = INAMEI +IC1
                IC2 = INAMEI +IC2
                IDIF = 0
                DO L = 1, L1
                  IF (WC(IC1+L-1)/=WC(IC2+L-1)) THEN
                    IDIF = IDIF+1
                  ENDIF
                ENDDO
                IF (IDIF==0) THEN
                  WRITE(*,*) 'Species ',I,' and', J, ' are identical...'
                  PAUSE
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDDO          
      ENDDO

C     Species files
      IC1 = 0; IC2 = 0; 
      DO I = 1, NS
        ICHECK = 0
        LSPECIES = WI(ILNAMEI+I-1) 
        TLSPECIES = LSPECIES
        DO J = 1, LSPECIES
          SPECIES(J:J) = WC(INAMEI+IC1+J-1) 
        ENDDO
        ISTAR = 0
        DO J = 1, LSPECIES
          IF (SPECIES(J:J)=="(") THEN
            TLSPECIES = J-1
            ISTAR = 1
          ENDIF
        ENDDO
        IC1 = IC1 + LSPECIES
        OPEN(UNIT=INOUT2,FILE=PATH(1:LPATH)//'/data/thermo/'//
     &       SPECIES(1:TLSPECIES),STATUS='OLD')
        COM2 = '   '
        DO WHILE (COM2(1:4)/= 'STOP')
           READ(INOUT2,*) FULLCOM2
           COM2 = FULLCOM2(1:4)
           IF ((COM2(1:4)=='Elec').AND.(SPECIES(1:2)/='em')) THEN
C       Non electronic-specific species
             IF (ISTAR==0) THEN
               READ(INOUT2,1002) WI(IELEI+I-1)
               IF (WI(IELEI+I-1)/=0) THEN
C       The proton Hp has no bound electrons  
C       It is also possible to neglect the electronic levels
                 READ(INOUT2,*)
                 DO J = 1, WI(IELEI+I-1)
                   READ(INOUT2,*) WI(IEGIK+IC2+J-1), WR1(IEEIK+IC2+J-1)
C       Conversion [cm**(-1)]*100*h*c/kb = [K] 
                   WR1(IEEIK+IC2+J-1) = 1.43876866D0 *WR1(IEEIK+IC2+J-1)
                 ENDDO
               ENDIF
C       Electronic-specific species
             ELSE
               WI(IELEI+I-1) = 1
               READ(INOUT2,1002) NCR
               READ(INOUT2,*)
               DO K = 1, NCR
                 CALL BLANK(WORD)
                 READ(INOUT2,*) DUM1, DUM2, WORD
                 L = LCHAR(WORD)
                 IF (L==LSPECIES) THEN   
                   IF (SPECIES(1:L)==WORD(1:L)) THEN
                     WI(IEGIK+IC2)  = DUM1
                     WR1(IEEIK+IC2) = DUM2
C       Conversion [cm**(-1)]*100*h*c/kb = [K] 
                     WR1(IEEIK+IC2) = 1.43876866D0 *WR1(IEEIK+IC2)
                   ENDIF
                 ENDIF
               ENDDO
             ENDIF
             IC2 = IC2 + WI(IELEI+I-1)
             ICHECK = ICHECK +1
           ELSEIF (COM2(1:4) == 'Mola') THEN
             READ(INOUT2,1004) WR1(IMI+I-1)
             ICHECK = ICHECK +1
           ELSEIF (COM2(1:4) == 'Form') THEN
             READ(INOUT2,1004) WR1(IHFORI+I-1)
             ICHECK = ICHECK +1
           ENDIF 
        ENDDO
        CLOSE(INOUT2)
        IF (SPECIES(1:2)=='em') THEN
          ICHECK = ICHECK +1
        ENDIF
        IF (ICHECK /= 3) THEN
           WRITE(*,*) 'Species file for ',SPECIES(1:LSPECIES),
     &                ' not correct.'
           RETURN
        ENDIF
      ENDDO

C     Molecules    
      IC1 = 0; IC2 = 0; IC3 = 0
      DO I = 1, NS
        ATOM = WI(IATOMI+I-1)
        LSPECIES = WI(ILNAMEI+I-1) 
        DO J = 1, LSPECIES
          SPECIES(J:J) = WC(INAMEI+IC1+J-1) 
        ENDDO
        IC1 = IC1 + LSPECIES
        IF ( ATOM > 1 ) THEN
          ICHECK = 0
          OPEN(UNIT=INOUT2,FILE=PATH(1:LPATH)//'/data/thermo/'//
     &         SPECIES(1:LSPECIES),STATUS='OLD')
          COM2 = '   '
          DO WHILE (COM2(1:4)/= 'STOP')
             READ(INOUT2,*) FULLCOM2
             COM2 = FULLCOM2(1:4)
             IF (COM2(1:4) == 'Line') THEN
               READ(INOUT2,1002) LIN 
               WI(ILINI+I-1) = LIN
               ICHECK = ICHECK +1
             ENDIF
          ENDDO
          CLOSE(INOUT2)
          IF ( LIN == 1 ) THEN
            NVIBMODE = 3 *ATOM -5
          ELSE
            NVIBMODE = 3 *ATOM -6 
          ENDIF
C     The number of rotational temperatures is always equal to 1.
          NROT = 1
          WI(IVIBI+I-1) = NVIBMODE
          WI(IROTI+I-1) = NROT
          OPEN(UNIT=INOUT2,FILE=PATH(1:LPATH)//'/data/thermo/'//
     &         SPECIES(1:LSPECIES),STATUS='OLD')
          COM2 = '   '
          DO WHILE (COM2(1:4)/= 'STOP')
             READ(INOUT2,*) FULLCOM2
             COM2 = FULLCOM2(1:4)
             IF ( COM2(1:4) == 'Vibr' ) THEN
               DO J = 1, NVIBMODE
                 READ(INOUT2,1004) WR1(ITVIK+IC2+J-1)
               ENDDO
               IC2 = IC2 + NVIBMODE
               ICHECK = ICHECK +1
             ELSEIF ( COM2(1:4) == 'Rota' ) THEN
               READ(INOUT2,1004) WR1(ITRIK+IC3)
               IC3 = IC3 + 1
               ICHECK = ICHECK +1
             ELSEIF ( COM2(1:4) == 'Ster' ) THEN
               READ(INOUT2,1002) WI(ISYMI+I-1)
               ICHECK = ICHECK +1
             ENDIF
          ENDDO
          CLOSE(INOUT2)
          IF (ICHECK /= 4) THEN
             WRITE(*,*) 'Species file for molecule ',
     &                  SPECIES(1:LSPECIES),' not correct.'
             RETURN
          ENDIF
        ENDIF
      ENDDO

C     Chemistry file
      OPEN(UNIT=INOUT1,FILE=PATH(1:LPATH)//'/data/chemistry/ceqreact/'
     &//MIXTURE(1:LMIXTURE)//'.ceq',STATUS='OLD')
      COM1 = '   '
      DO WHILE (COM1(1:4)/= 'STOP')
        READ(INOUT1,*) FULLCOM1
        COM1 = FULLCOM1(1:4)
        IF (COM1(1:4) == 'Defa') THEN
          DO I = 1 ,NC
            READ(INOUT1,*) WR1(IXN+I-1)
          ENDDO
        ELSEIF (COM1(1:4) == 'Temp') THEN
            READ(INOUT1,*) WR1(ITMIN)
        ENDIF
      ENDDO
      CLOSE(INOUT1)  
      
C     Check that both block diagonal of stoichiometric matrix are 
C     diagonal matrices
      DO I = 1, NC
        DO J = 1, NC
          IF ((J /= I).AND.(WI(INUIK-1+(I-1)*NS+J)/=0)) THEN
            WRITE(*,*) 'Block diagonal of stoichiometric matrix are not'
     &                  ,' diagonal matrices.'
          ENDIF
        ENDDO
      ENDDO
      DO I = 1, NR
        DO J =1, NR
          IF ((J /= I).AND.(WI(INUIK-1+(NC+I-1)*NS+NC+J)/=0)) THEN
            WRITE(*,*) 'Block diagonal of stoichiometric matrix are not'
     &                 ,' diagonal matrices.'
          ENDIF
        ENDDO
      ENDDO

C     Matrix for Schur complement (equilibrium composition)
      DO I = 1, NR
        DO J = 1, NC
          WR1(IJ22M1J21-1+(I-1)*(NC+1)+J) = WI(INUIK-1+(NC+I-1)*NS+J)
     &                            /WI(INUIK-1+(NC+I-1)*NS+NC+I)
        ENDDO 
        J = NC+1
        WR1(IJ22M1J21-1+(I-1)*(NC+1)+J) = 0.D0
      ENDDO

C     Collision integral fit coefficients : neutral-neutral 
C                                           and neutral-ion
      IF (IMOD == 1) THEN
        OPEN(UNIT=INOUT1,FILE=PATH(1:LPATH)//'/data/transport/heavy',
     &       STATUS='OLD')
        IC1 = 0; IF (NE /= 0) IC1 = WI(ILNAMEI); IC3 = 0
        DO I = NE +1, NE +NN
           LSPECIES1 = WI(ILNAMEI+I-1)
           DO K = 1, LSPECIES1
             SPECIES1(K:K) = WC(INAMEI+IC1+K-1)
           ENDDO
           IC2 = IC1; IC1 = IC1+ LSPECIES1
           DO J = I, NS
             LSPECIES2 = WI(ILNAMEI+J-1)
             DO K = 1, LSPECIES2
               SPECIES2(K:K) = WC(INAMEI+IC2+K-1)
             ENDDO
             IC2 = IC2 + LSPECIES2
             LTOT = 3 +LSPECIES1 +LSPECIES2
             NAME1(1:LTOT) = '.'//SPECIES1(1:LSPECIES1)//'-'//SPECIES2
     &                       (1:LSPECIES2)//'.'
             NAME2(1:LTOT) = '.'//SPECIES2(1:LSPECIES2)//'-'//SPECIES1
     &                       (1:LSPECIES1)//'.'
             FULLCOM1 = '                     '
             REWIND(INOUT1) ; IREAD = 0
             DO WHILE (FULLCOM1(1:4)/= 'STOP')
               READ(INOUT1,*) FULLCOM1
               IF ( (FULLCOM1(1:LTOT) == NAME1(1:LTOT) ) .OR. 
     &              (FULLCOM1(1:LTOT) == NAME2(1:LTOT) ) ) THEN
                 READ(INOUT1,*) ( WR1(IFITQ11IJ+IC3+K-1), K = 1, 4 )
                 READ(INOUT1,*) ( WR1(IFITAIJ+IC3+K-1), K = 1, 4 )
                 READ(INOUT1,*) ( WR1(IFITBIJ+IC3+K-1), K = 1, 4 )
                 READ(INOUT1,*) ( WR1(IFITCIJ+IC3+K-1), K = 1, 4 )
                 IREAD = IREAD +1
               ENDIF
             ENDDO
             IF  ( IREAD == 0) THEN
               WRITE(*,*) 'Interaction ',NAME1(1:LTOT),' not found'
               STOP
             ENDIF
             IF ( IREAD > 1) THEN
               WRITE(*,*) 'Problem with interaction ',NAME1(1:LTOT),
     &                    ' or ',NAME2(1:LTOT),' written ',IREAD, 
     &                    ' times'
               STOP
             ENDIF
             IC3 = IC3 +4
           ENDDO
        ENDDO
        CLOSE(INOUT1) 

C     Collision integral fit coefficients : electron-neutral
        IF (NE /= 0) THEN
          OPEN(UNIT=INOUT1,FILE=PATH(1:LPATH)//'/data/transport/electron
     &         ',STATUS='OLD')
          IC3 = 0
          I = 1
          LSPECIES1 = WI(ILNAMEI+I-1)
          DO K = 1, LSPECIES1
            SPECIES1(K:K) = WC(INAMEI+K-1)
          ENDDO
          IC2 = LSPECIES1;
          DO J = NE +1, NE +NN 
            LSPECIES2 = WI(ILNAMEI+J-1)
            DO K = 1, LSPECIES2
              SPECIES2(K:K) = WC(INAMEI+IC2+K-1)
            ENDDO
            IC2 = IC2 + LSPECIES2           
            LTOT = 3 +LSPECIES1 +LSPECIES2
            NAME1(1:LTOT) = '.'//SPECIES1(1:LSPECIES1)//'-'//SPECIES2
     &                      (1:LSPECIES2)//'.'
            FULLCOM1 = '                     '
            REWIND(INOUT1) ; IREAD = 0
            DO WHILE (FULLCOM1(1:4)/= 'STOP')
              READ(INOUT1,*) FULLCOM1
              IF ( FULLCOM1(1:LTOT) == NAME1(1:LTOT) ) THEN 
                 READ(INOUT1,*) ( WR1(IFITQ11EI+IC3+K-1), K = 1, 4 )
                 READ(INOUT1,*) ( WR1(IFITQ12EI+IC3+K-1), K = 1, 4 )
                 READ(INOUT1,*) ( WR1(IFITQ13EI+IC3+K-1), K = 1, 4 )
                 READ(INOUT1,*) ( WR1(IFITQ14EI+IC3+K-1), K = 1, 4 )
                 READ(INOUT1,*) ( WR1(IFITQ15EI+IC3+K-1), K = 1, 4 )
                 IREAD = IREAD +1
              ENDIF
            ENDDO
            IF ( IREAD == 0) THEN
              WRITE(*,*) 'Interaction ',NAME1(1:LTOT),' not found'
              STOP
            ENDIF
            IF ( IREAD > 1) THEN
              WRITE(*,*) 'Problem with interaction ',NAME1(1:LTOT),
     &                   ' written ',IREAD, ' times'
              STOP
            ENDIF
            IC3 = IC3 +4
          ENDDO
          CLOSE(INOUT1) 
        ENDIF
      ENDIF

C     Reaction rates
      IF (REACTION(1:5) /= 'empty') THEN
        CALL READRATES (PATH, REACTION, WI, LWI, WC, LWC, WR1, LWR1, 
     &                  WR3, LWR3)
      ENDIF

C     Species considered in the energy transfer equations and 
C     temperatures.
      IF (TRANSFER(1:5) /= 'empty') THEN
        OPEN(UNIT=INOUT1,FILE=PATH(1:LPATH)//'/data/transfer/'
     &       //TRANSFER(1:LTRANSFER),STATUS='OLD')
        COM1 = '   '
        DO WHILE (COM1(1:4)/= 'STOP')
          CALL BLANK (FULLCOM1)
          READ(INOUT1,*) FULLCOM1
          COM1 = FULLCOM1(1:4)
          IF (COM1(1:4) == 'Vibr') THEN
            READ(INOUT1,*)
            READ(INOUT1,*)
            CALL BLANK (SPECIES)
            IC = 0; IV = 0
            DO IVIB = 1, NVIB
              CALL BLANK (FULLCOM1)
              READ(INOUT1,'(A)') FULLCOM1
              CALL RMBLANK (FULLCOM1, LINE, ILEN)
              WI(INVIBSPEI+IVIB-1) = 0
              DO I = 1, ILEN
                IC = IC +1
                IF ((LINE(I:I)=='/').OR.(I==ILEN)) THEN
                  IF (I==ILEN) THEN
                    SPECIES(IC:IC) = LINE(I:I)                
                  ENDIF
                  CALL TESTSPECIES (SPECIES, WI, LWI, WC, LWC, ISPECIES)
                  IF (ISPECIES==0) THEN
                    WRITE(*,*) 'Species ', SPECIES,
     &                         ' does not exist!', LINE
                    PAUSE
                  ELSE
                    IF (WI(IATOMI + ISPECIES-1) > 1) THEN
                      WI(INVIBSPEI+IVIB-1) = WI(INVIBSPEI+IVIB-1) +1
                      IV = IV +1
                      WI(IVIBSPEI+IV-1)  = ISPECIES
                    ELSE
                      WRITE(*,*) 'Species ', SPECIES,
     &                           ' is not a molecule'
                      PAUSE
                    ENDIF 
                  ENDIF
                  IC = 0
                  CALL BLANK (SPECIES)
                ELSE
                  SPECIES(IC:IC) = LINE(I:I)
                ENDIF
              ENDDO      
            ENDDO
          ENDIF
        ENDDO
        CLOSE(INOUT1)
C     Check whether the molecules are entered more than once
          IC1 = 0
          DO IV = 1, NVIB
            DO I = 1, WI(INVIBSPEI+IV-1)
              IC1 = IC1 +1
              IC2 = 0
              DO IV2 = 1, NVIB
                DO I2 = 1, WI(INVIBSPEI+IV2-1)
                  IC2 = IC2 +1 
                    IF ( (WI(IVIBSPEI+IC1-1) == WI(IVIBSPEI+IC2-1))
     &                    .AND. (IC1 /= IC2) ) THEN
                      WRITE(*,*) 'Molecule entered twice in the ',
     &                           ' tranfer file'
                      PAUSE
                    ENDIF
                ENDDO 
              ENDDO             
            ENDDO 
          ENDDO 
C     If the species IS is not a molecule WI(IVIBTEMPI+IS-1) = 1
        DO IS = 1, NS
          WI(IVIBTEMPI+IS-1) = 1
        ENDDO
        IF (NVIB >= 1) THEN
C     Default vibrational temperature of the first vibrational equation
          DO IS = 1, NS
            IF (WI(IATOMI+IS-1) > 1) THEN
              WI(IVIBTEMPI+IS-1) = 2
            ENDIF
          ENDDO 
C     Vibrational temperature specified in the transfer data file
          IC1 = 0
          DO IV = 1, NVIB
            DO I = 1, WI(INVIBSPEI+IV-1)
              IC1 = IC1 +1
              WI(IVIBTEMPI+WI(IVIBSPEI+IC1-1)-1) = 1+IV
            ENDDO 
          ENDDO
        ENDIF
C     VT relaxation time
        CALL INITIALIZETRANSFER (PATH, LWR1, WR1, LWR4, WR4, LWI, WI, 
     &                           LWC, WC)
      ENDIF

1002  FORMAT (I5) 
1003  FORMAT (A) 
1004  FORMAT (f20.6) 

      END SUBROUTINE INITIALIZE
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE CUTMOLECULE(SPECIES, LSPECIES, NUCLEI, SPECIESNUCLEI,
     &                       NC)
C-----------------------------------------------------------------------
C     This subroutines determines the number of atoms or charges
C     corresponding to a given nucleus for SPECIES (maximum 9 atoms of 
C     a given nucleus). The charges are denoted by the symbols "p" 
C     (positive charge) and "m" (negative charge), these symbols are 
C     repeated when molecules exhibit multiple charges. 
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER LSPECIES, NC, SPECIESNUCLEI(1:NC)
      CHARACTER(10) SPECIES, NUCLEI
C-----------------------------------------------------------------------
      INTEGER I, J, JMEM, INUCLEUS, NNUCLEI, TLSPECIES
C-----------------------------------------------------------------------
      TLSPECIES = LSPECIES
      DO I = 1, LSPECIES
        IF (SPECIES(I:I)=="(") THEN
          TLSPECIES = I-1
        ENDIF
      ENDDO
      DO I = 1, NC
        SPECIESNUCLEI(I) = 0
      ENDDO
      DO I = 1, TLSPECIES
        INUCLEUS = 0
        DO J = 1, NC
          IF (SPECIES(I:I)==NUCLEI(J:J)) THEN
            SPECIESNUCLEI(J) = 1
            JMEM = J
            INUCLEUS = 1
          ENDIF
        ENDDO
        IF (INUCLEUS == 0) THEN
          IF (SPECIES(I:I)=='p') THEN
            SPECIESNUCLEI(1) = SPECIESNUCLEI(1) +1
          ELSEIF (SPECIES(I:I)=='m') THEN
            SPECIESNUCLEI(1) = SPECIESNUCLEI(1) -1
          ELSE
            CALL ITRANSLATE(SPECIES(I:I), NNUCLEI)
            IF ((NNUCLEI>0).AND.(NNUCLEI<10)) THEN
              SPECIESNUCLEI(JMEM) = NNUCLEI
            ELSE
              WRITE(*,*) 'Unknown ', SPECIES(1:LSPECIES), ' molecule...'
              PAUSE
            ENDIF
          ENDIF
        ENDIF
      ENDDO

      END SUBROUTINE CUTMOLECULE
C-----------------------------------------------------------------------
