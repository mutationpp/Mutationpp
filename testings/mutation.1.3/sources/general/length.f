C-----------------------------------------------------------------------
      SUBROUTINE LENGTH (PATH, MIXTURE, REACTION, TRANSFER, LWR1, LWR2,
     &                   LWR3, LWR4, LWI, LWC, NNS, NNE, NNC, NNREA, 
     &                   NNV, NMAX, NNVIB, NNELEQ)
C-----------------------------------------------------------------------
C     This subroutine computes the length of the work arrays
C     WR1, WR2, WR3, WI and WC and initializes their pointers. 
C     This subroutine should be called by the user once at the
C     beginning of the program.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER LWR1, LWR2, LWR3, LWR4, LWI, LWC, NNS, NNE, NNC, NNREA, 
     &        LPATH, NNV, NMAX, NNVIB, NNELEQ
      CHARACTER(10) MIXTURE, REACTION, TRANSFER
      CHARACTER(100) PATH 
C-----------------------------------------------------------------------
      INTEGER S, LSPECIES, LCHAR, CHARGE, ATOM, LIN, NLEVEL, 
     &        NNSYM, NSSYM, NHSYM, NA, NAN, NAI, NP, NPN, NPI,
     &        LMIXTURE, LREACTION, LTRANSFER, I, ILEN, IC, IM, IVIB, I1,
     &        SPECIESNUCLEI(1:9), J, NCR, LWORD, DUM1
      DOUBLE PRECISION COEFFICIENT, DUM2

      CHARACTER(4) COM1, COM2 
      CHARACTER(10) WORD, NUCLEI, SPECIES
      CHARACTER(80) FULLCOM1, FULLCOM2, LINE
      INCLUDE '../general/memory.cmn'
C-----------------------------------------------------------------------
      LWR1 = 1; LWR2 = 1; LWR3 = 1; LWR4 = 1; LWI = 1; LWC = 1
      LPATH      = LCHAR (PATH)
      LMIXTURE   = LCHAR (MIXTURE)
      LREACTION  = LCHAR (REACTION)
      LTRANSFER = LCHAR (TRANSFER)

C     Input-output
      INOUT1 = 20; INOUT2 = 21; INOUT3 = 22; INOUT4 = 23;  
 
C     Universal constants
      IUKB = LWR1; LWR1 = LWR1 + 1
      IUNA = LWR1; LWR1 = LWR1 + 1
      IUR  = LWR1; LWR1 = LWR1 + 1
      IUE  = LWR1; LWR1 = LWR1 + 1
      IUPI = LWR1; LWR1 = LWR1 + 1
      IUE0 = LWR1; LWR1 = LWR1 + 1
      IUH  = LWR1; LWR1 = LWR1 + 1
      IUQT = LWR1; LWR1 = LWR1 + 1

C     Number of species      
      OPEN(UNIT=INOUT1,FILE=PATH(1:LPATH)//'/data/mixture/'
     &     //MIXTURE(1:LMIXTURE)//'.mix',STATUS='OLD')
      COM1 = '   '
      DO WHILE (COM1(1:4)/= 'STOP')
        READ(INOUT1,*) FULLCOM1
        COM1 = FULLCOM1(1:4)
        IF (COM1(1:4) == 'Numb') THEN
          READ(INOUT1,1002) NS
        ENDIF 
      ENDDO
      NNS = NS
      REWIND (UNIT=INOUT1) 
      COM1 = '   '
      NSTAR = 0
      DO WHILE (COM1(1:4)/= 'STOP')
        CALL BLANK (FULLCOM1)
        READ(INOUT1,*) FULLCOM1
        COM1 = FULLCOM1(1:4)
        IF (COM1(1:4) == 'Spec') THEN
          DO I = 1, NNS
            CALL BLANK (SPECIES)
            READ(INOUT1,1003) SPECIES 
            LSPECIES  = LCHAR (SPECIES)
            IF (LSPECIES>5) THEN
              J = LSPECIES-4
              IF (SPECIES(J:LSPECIES)=='-star') THEN 
                OPEN(UNIT=INOUT2,FILE=PATH(1:LPATH)//'/data/thermo/'//
     &               SPECIES(1:J-1),STATUS='OLD')
                COM2 = '   '
                DO WHILE (COM2(1:4)/= 'STOP')
                  READ(INOUT2,*) FULLCOM2
                  COM2 = FULLCOM2(1:4)
                  IF (COM2(1:4)=='Elec') THEN
                    READ(INOUT2,*) NCR
                    NSTAR = NSTAR +NCR-1
                    NS  = NS +NCR -1
                  ENDIF
                ENDDO
                CLOSE(INOUT2)
              ENDIF
            ENDIF           
          ENDDO
        ENDIF 
      ENDDO
      CLOSE(INOUT1) 
      NNS = NS

C     Length of the species name
      ILNAMEI = LWI; LWI = LWI + NS 

C     Number of independent chemical reactions
      OPEN(UNIT=INOUT1,FILE=PATH(1:LPATH)//'/data/chemistry/ceqreact/'
     &     //MIXTURE(1:LMIXTURE)//'.ceq',STATUS='OLD')
      COM1 = '   '
      DO WHILE (COM1(1:4)/= 'STOP')
        READ(INOUT1,*) FULLCOM1
        COM1 = FULLCOM1(1:4)
        IF (COM1(1:4) == 'Numb') THEN
          READ(INOUT1,1002) NC
          NNC = NC
        ENDIF
      ENDDO
      CLOSE(INOUT1)
      IF (NC > 9) THEN
        WRITE(*,*) 'Number of nuclei:', NC, ' > 9, ',
     &             'modify the CUTMOLECULE routine'
        PAUSE
      ENDIF    
      NR  = NS  -NC

C     Species name
      INAMEI = LWC
      OPEN(UNIT=INOUT1,FILE=PATH(1:LPATH)//'/data/mixture/'
     &     //MIXTURE(1:LMIXTURE)//'.mix',STATUS='OLD')
      COM1 = '   '
      DO WHILE (COM1(1:4)/= 'STOP')
        READ(INOUT1,*) FULLCOM1
        COM1 = FULLCOM1(1:4)
        IF (COM1(1:4) == 'Spec') THEN
          DO S = 1, NC
            CALL BLANK (SPECIES)
            READ(INOUT1,1003) SPECIES 
            LSPECIES  = LCHAR (SPECIES)
            NUCLEI(S:S) = SPECIES(1:1)
            IF (LSPECIES==1) THEN
              LWC = LWC +LSPECIES
            ELSEIF ((LSPECIES==2).AND.(SPECIES(1:2)=='em')) THEN
              LWC = LWC +LSPECIES
              NE = 1
            ELSEIF ((LSPECIES==6).AND.(SPECIES(2:6)=='-star')) THEN
              OPEN(UNIT=INOUT2,FILE=PATH(1:LPATH)//'/data/thermo/'//
     &             SPECIES(1:1),STATUS='OLD')
              COM2 = '   '
              DO WHILE (COM2(1:4)/= 'STOP')
                READ(INOUT2,*) FULLCOM2
                COM2 = FULLCOM2(1:4)
                IF (COM2(1:4)=='Elec') THEN
                  READ(INOUT2,*) NCR
                  READ(INOUT2,*) 
                  CALL BLANK(WORD)
                  READ(INOUT2,*) DUM1, DUM2, WORD   
                  LWORD = LCHAR(WORD)
                  LWC = LWC + LWORD
                ENDIF
              ENDDO
              CLOSE(INOUT2)
            ELSE
                WRITE(*,*) 'Nucleus not found'
            ENDIF
          ENDDO
          DO I = 1, NC
            DO J = 1, NC
              IF ((NUCLEI(I:I)==NUCLEI(J:J)).AND.(I/=J)) THEN
                WRITE(*,*)  I, J, NUCLEI(I:I), NUCLEI(J:J)
                WRITE(*,*) 'Twice the same nucleus'
                PAUSE
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      REWIND(UNIT=INOUT1)
      COM1 = '    '
      DO WHILE (COM1(1:4)/= 'STOP')
        CALL BLANK (FULLCOM1)
        READ(INOUT1,*) FULLCOM1
        COM1 = FULLCOM1(1:4)
        IF (COM1(1:4) == 'Spec') THEN
          DO S = 1, NS-NSTAR
            CALL BLANK (SPECIES)
            READ(INOUT1,1003) SPECIES 
            LSPECIES  = LCHAR (SPECIES)
C 12
            IF (LSPECIES>5) THEN
              J = LSPECIES-4
              IF (SPECIES(J:LSPECIES)=='-star') THEN 
                OPEN(UNIT=INOUT2,FILE=PATH(1:LPATH)//'/data/thermo/'//
     &               SPECIES(1:J-1),STATUS='OLD')
                COM2 = '   '
                DO WHILE (COM2(1:4)/= 'STOP')
                  READ(INOUT2,*) FULLCOM2
                  COM2 = FULLCOM2(1:4)
                  IF (COM2(1:4)=='Elec') THEN
                    READ(INOUT2,*) NCR
                    READ(INOUT2,*)
                    IF (S<=NC) THEN
C     Ground state nucleus already taken into account
                      READ(INOUT2,*)
                      I1 = 2
                    ELSE
C     Cround state is not a nucleus
                      I1 = 1
                    ENDIF
                    DO I = I1, NCR 
                      CALL BLANK(WORD)
                      READ(INOUT2,*) DUM1, DUM2, WORD   
                      LWORD = LCHAR(WORD)
                      LWC = LWC + LWORD
                    ENDDO
                  ENDIF
                ENDDO
                CLOSE(INOUT2)
              ELSE
C     Species with a name length > 5 but not electronic specific
                LWC = LWC +LSPECIES
              ENDIF
C 12
            ELSE
C     Species with a name length <= 5
              IF (S > NC) THEN
                LWC = LWC +LSPECIES
              ENDIF
            ENDIF                       
          ENDDO
        ENDIF 
      ENDDO
      CLOSE(INOUT1) 

C     Mass
      IMI = LWR1; LWR1 = LWR1 +NS 
      
C     Charge
      IQI = LWI; LWI = LWI +NS 

C     Number of constitutive atoms: ATOM =  0, electron
C                                        =  1, atom
C                                        = >1, polyatomic 
      IATOMI = LWI; LWI = LWI +NS 

C     Formation enthalpy 
      IHFORI = LWR1; LWR1 = LWR1 +NS

C     Number of electronic levels
      IELEI = LWI; LWI = LWI +NS

C     Electronic level energy
      IEEIK = LWR1 

C     Electronic level degeneracy
      IEGIK = LWI
      OPEN(UNIT=INOUT1,FILE=PATH(1:LPATH)//'/data/mixture/'
     &     //MIXTURE(1:LMIXTURE)//'.mix',STATUS='OLD')
      COM1 = '   '
      DO WHILE (COM1(1:4)/= 'STOP')
        READ(INOUT1,*) FULLCOM1
        COM1 = FULLCOM1(1:4)
        IF (COM1(1:4) == 'Spec') THEN
          DO S = 1, NS-NSTAR
            READ(INOUT1,1003) SPECIES 
            LSPECIES  = LCHAR (SPECIES)
            IF (LSPECIES>5) THEN
              J = LSPECIES-4
              IF (SPECIES(J:LSPECIES)=='-star') THEN
              ELSE
                J = LSPECIES +1
              ENDIF
            ELSE
              J = LSPECIES +1
            ENDIF   
            OPEN(UNIT=INOUT2,FILE=PATH(1:LPATH)//'/data/thermo/'//
     &           SPECIES(1:J-1),STATUS='OLD')
            COM2 = '   '; NLEVEL = 0
            DO WHILE (COM2(1:4)/= 'STOP')
              CALL BLANK (FULLCOM2)
              READ(INOUT2,*) FULLCOM2
              COM2 = FULLCOM2(1:4)
              IF (COM2(1:4) == 'Elec') THEN
                READ(INOUT2,1002) NLEVEL
              ENDIF 
            ENDDO
            LWR1 = LWR1 +NLEVEL
            LWI = LWI +NLEVEL 
            CLOSE(INOUT2)
          ENDDO
        ENDIF
      ENDDO
      CLOSE(INOUT1)

C     Linearity of polyatomic molecules: LIN = 1, linear
C                                              0, non linear
      ILINI = LWI; LWI = LWI +NS 

C     Symmetry number: SYM = 1, non symmetric molecule
C                          = 2, symmetric molecule
      ISYMI = LWI; LWI = LWI +NS 

C     Number of vibrational temperature
      IVIBI = LWI; LWI = LWI +NS

C     Vibrational temperature(s)
      NV = 0 ; ITVIK = LWR1 
      NAN = 0; NAI = 0; NPN = 0; NPI = 0
      OPEN(UNIT=INOUT1,FILE=PATH(1:LPATH)//'/data/mixture/'
     &     //MIXTURE(1:LMIXTURE)//'.mix',STATUS='OLD')
      COM1 = '   '
      DO WHILE (COM1(1:4)/= 'STOP')
        READ(INOUT1,*) FULLCOM1
        COM1 = FULLCOM1(1:4)
        IF (COM1(1:4) == 'Spec') THEN
          DO S = 1, NS-NSTAR
            CALL BLANK (SPECIES)
            READ(INOUT1,1003) SPECIES 
            LSPECIES  = LCHAR (SPECIES)
C     Nucleus (including excited states)
            IF (S<=NC)THEN
              IF ((S==1).AND.(NE==1)) THEN
                J      =  2
                CHARGE = -1
                ATOM   =  0
              ELSE
                J      =  1
                CHARGE =  0
                ATOM   =  1  
                IF (SPECIES(2:LSPECIES)=='-star') THEN
                  OPEN(UNIT=INOUT2,FILE=PATH(1:LPATH)//'/data/thermo/'//
     &                 SPECIES(1:J),STATUS='OLD')
                  COM2 = '   '
                  DO WHILE (COM2(1:4)/= 'STOP')
                    READ(INOUT2,*) FULLCOM2
                    COM2 = FULLCOM2(1:4)
                    IF (COM2(1:4)=='Elec') THEN
                      READ(INOUT2,*) NCR
                    ENDIF
                  ENDDO 
                  CLOSE(INOUT2)
                ELSE
                  NCR = 1
                ENDIF
              ENDIF
C     Other species 
            ELSE
              IF (LSPECIES>5) THEN
                J = LSPECIES-4
                IF (SPECIES(J:LSPECIES)=='-star') THEN
                  OPEN(UNIT=INOUT2,FILE=PATH(1:LPATH)//'/data/thermo/'//
     &                 SPECIES(1:J-1),STATUS='OLD')
                  COM2 = '   '
                  DO WHILE (COM2(1:4)/= 'STOP')
                    READ(INOUT2,*) FULLCOM2
                    COM2 = FULLCOM2(1:4)
                    IF (COM2(1:4)=='Elec') THEN
                      READ(INOUT2,*) NCR
                    ENDIF
                  ENDDO 
                  CLOSE(INOUT2)
                  J = LSPECIES-5
                ELSE
                  J   = LSPECIES
                  NCR = 1
                ENDIF
              ELSE
                J   = LSPECIES
                NCR = 1
              ENDIF  
              CALL CUTMOLECULE(SPECIES, J, NUCLEI, SPECIESNUCLEI,
     &                         NC)
              ATOM = 0 
              DO I = NE+1, NC 
                ATOM = ATOM +SPECIESNUCLEI(I)
              ENDDO
              IF (NE.EQ.0) THEN
                CHARGE = 0
              ELSE    
                CHARGE = SPECIESNUCLEI(1)
              ENDIF
            ENDIF
            OPEN(UNIT=INOUT2,FILE=PATH(1:LPATH)//'/data/thermo/'
     &          //SPECIES(1:J),STATUS='OLD')
            COM2 = '   '
            DO WHILE (COM2(1:4)/= 'STOP')
              READ(INOUT2,*) FULLCOM2
              COM2 = FULLCOM2(1:4)
              IF (COM2(1:4) == 'Line') THEN
                READ(INOUT2,1002) LIN  
              ENDIF 
            ENDDO
            IF ( ATOM > 1 ) THEN
              IF ( LIN == 1 ) THEN
                LWR1 = LWR1 +3 *ATOM -5
                NV = NV +3 *ATOM -5
              ELSE
                LWR1 = LWR1 +3 *ATOM -6
                NV = NV +3 *ATOM -6
              ENDIF
            ENDIF
            CLOSE(INOUT2)
            SELECT CASE(ATOM)
               CASE(0)                   !Electrons
               CASE(1) 
                 IF (CHARGE == 0 ) THEN
                   NAN = NAN +NCR          !Neutral atoms
                 ELSE
                   NAI = NAI +NCR          !Ion atoms
                 ENDIF
               CASE DEFAULT
                 IF (CHARGE == 0 ) THEN
                   NPN = NPN +NCR          !Neutral polyatomic
                 ELSE
                   NPI = NPI +NCR          !Ion polyatomic
                 ENDIF
             END SELECT
             CLOSE(INOUT2)
          ENDDO
        ENDIF
      ENDDO
      CLOSE(INOUT1)
      NNV = NV

C     Number of rotational temperatures
      IROTI = LWI; LWI = LWI +NS

C     Rotational temperature
      ITRIK = LWR1
      OPEN(UNIT=INOUT1,FILE=PATH(1:LPATH)//'/data/mixture/'
     &     //MIXTURE(1:LMIXTURE)//'.mix',STATUS='OLD')
      COM1 = '   '
      DO WHILE (COM1(1:4)/= 'STOP')
        READ(INOUT1,*) FULLCOM1
        COM1 = FULLCOM1(1:4)
        IF (COM1(1:4) == 'Spec') THEN
          DO S = 1, NS-NSTAR
            READ(INOUT1,1003) SPECIES 
            LSPECIES  = LCHAR (SPECIES)
            IF (LSPECIES>5) THEN
              J = LSPECIES-4
              IF (SPECIES(J:LSPECIES)=='-star') THEN
                J = LSPECIES-5
              ELSE
                J   = LSPECIES
              ENDIF
            ELSE
              J   = LSPECIES
            ENDIF  
            OPEN(UNIT=INOUT2,FILE=PATH(1:LPATH)//'/data/thermo/'//
     &           SPECIES(1:J),STATUS='OLD')
            COM2 = '   '
            DO WHILE (COM2(1:4)/= 'STOP')
              READ(INOUT2,*) FULLCOM2
              COM2 = FULLCOM2(1:4)
              IF (COM2(1:4) == 'Cons') THEN
                READ(INOUT2,1002) ATOM
              ELSEIF (COM2(1:4) == 'Line') THEN
                READ(INOUT2,1002) LIN  
              ENDIF 
            ENDDO
            IF ( ATOM > 1 ) THEN
              LWR1 = LWR1 +1
            ENDIF
            CLOSE(INOUT2)
          ENDDO
        ENDIF
      ENDDO
      CLOSE(INOUT1)

C     Number of atoms
      NA  = NAN +NAI

C     Number of polyatomic molecules
      NP  = NPN +NPI

C     Number of heavy particles
      NH  = NA  +NP

C     Number of neutral particles
      NN = NAN +NPN

C     Number of ions (charged atoms and polyatomic molecules)
      NI = NAI +NPI

      NNE = NE
      NNSYM = NN *(NN +1) /2
      NHSYM = NH *(NH +1) /2
      NSSYM = NS *(NS +1) /2     

C     Temperature treshold for a robust computation of composition
      ITMIN = LWR1; LWR1 = LWR1 +1

C     Stoichiometric matrix
      INUIK  = LWI; LWI = LWI +NS *NS 

C     Matrix for Schur complement (equilibrium composition)
      IJ22M1J21 = LWR1; LWR1 = LWR1 +NR*(NC+1)

C     Default nuclear fraction for chemical equilibrium composition
      IXN = LWR1; LWR1 = LWR1 +NC

C     Curve fit coefficients for the collision integrals
      IFITQ11IJ = LWR1; LWR1 = LWR1 +4 *(NNSYM +NN *NI) 
      IFITAIJ   = LWR1; LWR1 = LWR1 +4 *(NNSYM +NN *NI)
      IFITBIJ   = LWR1; LWR1 = LWR1 +4 *(NNSYM +NN *NI)
      IFITCIJ   = LWR1; LWR1 = LWR1 +4 *(NNSYM +NN *NI)
      IFITQ11EI = LWR1; LWR1 = LWR1 +4 *NN
      IFITQ12EI = LWR1; LWR1 = LWR1 +4 *NN
      IFITQ13EI = LWR1; LWR1 = LWR1 +4 *NN
      IFITQ14EI = LWR1; LWR1 = LWR1 +4 *NN
      IFITQ15EI = LWR1; LWR1 = LWR1 +4 *NN
   
C     Species viscosity
      IETAI  = LWR2; LWR2 = LWR2 +NH
      IQEE22 = LWR2; LWR2 = LWR2 +1

C     Modified binary diffusion coefficient
      IBINIJ = LWR2; LWR2 = LWR2 +NSSYM

C     Collision integral ratios
      IAIJ      = LWR2; LWR2 = LWR2 +NHSYM    
      IBIJ      = LWR2; LWR2 = LWR2 +NHSYM    
      ICIJ      = LWR2; LWR2 = LWR2 +NHSYM   

C     Devoto collision integrals for electron gas
      IDQ00EE   = LWR2; LWR2 = LWR2 +NE
      IDQ10EE   = LWR2; LWR2 = LWR2 +NE
      IDQ20EE   = LWR2; LWR2 = LWR2 +NE
      IDQ11EE   = LWR2; LWR2 = LWR2 +NE
      IDQ12EE   = LWR2; LWR2 = LWR2 +NE
      IDQ22EE   = LWR2; LWR2 = LWR2 +NE
      IDQ10EI   = LWR2; LWR2 = LWR2 +NH
      IDQ20EI   = LWR2; LWR2 = LWR2 +NH

C     Number of reactions (chemical nonequilibrium, Arrhenius's law)
      NREA = 0
      IF (REACTION(1:5) /= 'empty') THEN
        OPEN(UNIT=INOUT3,FILE=PATH(1:LPATH)//'/data/chemistry/gasreact/'
     &       //REACTION(1:LREACTION),STATUS='OLD')
        CALL BLANK (FULLCOM1)
        READ(INOUT3,'(A)') FULLCOM1
        CALL RMBLANK (FULLCOM1, LINE, ILEN)
        COM1 = LINE(1:4)
        DO WHILE (COM1(1:4) /= 'STOP')
          IC = 0
          IM = 0
          I = 1
          DOWHILE (LINE(I:I)/='/')
            I = I +1
            IC = IC +1
            IF ((LINE(I:I)=='+').OR.(LINE(I:I)=='=')
     &          .OR.(LINE(I:I)==':')) THEN
              CALL CUTWORD (WORD, COEFFICIENT, SPECIES)
              IF (SPECIES(1:2)=='M ') THEN
                IM = IM +1
              ENDIF
              CALL BLANK(SPECIES)
              CALL BLANK (WORD) 
              IC = 0
            ENDIF
            IF (IC/=0) THEN
              WORD(IC:IC) = LINE(I:I)
            ENDIF
          ENDDO
          CALL BLANK (FULLCOM1)
          READ(INOUT3,'(A)') FULLCOM1
          CALL RMBLANK (FULLCOM1, LINE, ILEN)
          COM1 = FULLCOM1(1:4)
          IF (IM>0) THEN
C     => M = some specific species
            IF ((LINE(1:1)=='M').AND.(LINE(2:2)=='=')) THEN
              NREA = NREA +1
              CALL BLANK (FULLCOM1)
              READ(INOUT3,'(A)') FULLCOM1
              CALL RMBLANK (FULLCOM1, LINE, ILEN)
              COM1 = FULLCOM1(1:4)
C     => M = all species
            ELSE
              NREA = NREA +1
            ENDIF
            IM = 0 
C     => Noncatalytic reaction
          ELSE
            NREA = NREA +1
          ENDIF
        ENDDO
        CLOSE(INOUT3)
      ENDIF
      NNREA = NREA

C     Reaction rates
      IRATE = LWR3;   LWR3 = LWR3 +NREA *4
  
C     Reaction thermal nonequilibrium mechanism
      ITNEQ = LWI;    LWI  = LWI  +NREA

C     Molecules dissociated by atomic or molecular impact
      IDISSOCIATION = LWI;    LWI  = LWI  +NREA

C     Stoichiometric coefficients
      ISTOIR = LWR3;  LWR3 = LWR3 +NREA *NS
      ISTOIP = LWR3;  LWR3 = LWR3 +NREA *NS

C     Catalytic reactions
      ICATA  = LWI;   LWI  = LWI +NREA
      IALPHA = LWR3;  LWR3 = LWR3 +NREA *NS

C     
      IF (TRANSFER(1:5) == 'empty') THEN
        NNVIB  = 0
        NNELEQ = 0
        IC     = 0
      ELSE
        OPEN(UNIT=INOUT1,FILE=PATH(1:LPATH)//'/data/transfer/'
     &       //TRANSFER(1:LTRANSFER),STATUS='OLD')
        COM1 = '   '
        DO WHILE (COM1(1:4)/= 'STOP')
          READ(INOUT1,*) FULLCOM1
          COM1 = FULLCOM1(1:4)
          IF (COM1(1:4) == 'Vibr') THEN
            READ(INOUT1,1002) NNVIB
            IC = 0
            READ(INOUT1,*)
            DO IVIB = 1, NNVIB
              CALL BLANK (FULLCOM1)
              READ(INOUT1,'(A)') FULLCOM1
              CALL RMBLANK (FULLCOM1, LINE, ILEN)
              IC = IC +1
              DO I = 1, ILEN
                IF (LINE(I:I)=='/') THEN
                  IC = IC +1
                ENDIF
              ENDDO
            ENDDO
            INVIBSPEI  = LWI; LWI = LWI +NVIB
            IVIBSPEI   = LWI; LWI = LWI +IC
            IVIBTEMPI  = LWI; LWI = LWI +NS
          ELSEIF (COM1(1:4) == 'Elec') THEN
            READ(INOUT1,1002) NNELEQ
            NELEQ = NNELEQ
          ENDIF 
        ENDDO
        IF (NNELEQ /= 0) THEN
          WRITE(*,*) 'Electron temperature not implemented yet'
          PAUSE
        ENDIF
        CLOSE(INOUT1)
      ENDIF
      NVIB = NNVIB
     
      IMWAIJ = LWR4; LWR4 = LWR4 +IC*(NS-NE)
      IMWBIJ = LWR4; LWR4 = LWR4 +IC*(NS-NE)

      LWR1 = LWR1 -1; LWR2 = LWR2 -1; LWR3 = LWR3 -1
      LWR4 = LWR4 -1; LWI  = LWI  -1; LWC  = LWC  -1  

1002  FORMAT (I5)
1003  FORMAT (A)

      END SUBROUTINE LENGTH
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION LCHAR (SYMBOL)
C-----------------------------------------------------------------------
C     This function determines the effective length (non blank) 
C     of a character 
C-----------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(*) SYMBOL
      INTEGER LCHAR
C-----------------------------------------------------------------------
      INTEGER C, I 
C-----------------------------------------------------------------------
      C = 0 
      DO I = 1, LEN(SYMBOL)
         IF ( SYMBOL(I:I) /= " " ) THEN
            C = C +1
         ENDIF
      ENDDO
      
      LCHAR = C 

      END FUNCTION LCHAR
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE RMBLANK (SYMBOLIN, SYMBOLOUT, ILEN)
C-----------------------------------------------------------------------
C     This subroutine  selects the non blank characters of SYMBOLIN 
C     (number = ILEN) and put them in SYMBOLOUT.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER ILEN
      CHARACTER(*) SYMBOLIN, SYMBOLOUT
C-----------------------------------------------------------------------
      INTEGER I
C-----------------------------------------------------------------------
      ILEN = 0
      DO I = 1, LEN(SYMBOLIN)
        IF (SYMBOLIN(I:I) /= " ") THEN
          ILEN = ILEN +1
          SYMBOLOUT(ILEN:ILEN) = SYMBOLIN(I:I)
        ENDIF
      ENDDO 
      DO I =  ILEN +1, LEN(SYMBOLIN)
        SYMBOLOUT(I:I) = " "
      ENDDO

      END SUBROUTINE RMBLANK
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE BLANK (SYMBOL)
C-----------------------------------------------------------------------
C     This subroutine fill SYMBOL with blanks.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(*) SYMBOL
C-----------------------------------------------------------------------
       INTEGER I
C-----------------------------------------------------------------------
      DO I = 1, LEN (SYMBOL)
        SYMBOL(I:I) = " "
      ENDDO

      END SUBROUTINE BLANK
C-----------------------------------------------------------------------
