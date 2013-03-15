C-----------------------------------------------------------------------
      SUBROUTINE PRESENTATION (I) 
C-----------------------------------------------------------------------
C     Presentation screen
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER I, SCREEN 
C-----------------------------------------------------------------------
      SCREEN = I
      SELECT CASE(SCREEN)
        CASE(1)
          WRITE(*,*) ''
          WRITE(*,*) '    MUTATION 1.0 is running in thermo-chemical ',
     &                    'equilibrium...'
          WRITE(*,*) ''
        CASE(2)
          WRITE(*,*) ''
          WRITE(*,*) '    Good bye! Many thanks!'
          WRITE(*,*) ''
        CASE(3)
          WRITE(*,*) ''
          WRITE(*,*) '    SHOCKING is running...'
      END SELECT 
C-----------------------------------------------------------------------
      END SUBROUTINE PRESENTATION
C-----------------------------------------------------------------------

