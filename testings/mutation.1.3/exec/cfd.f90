      module mutvar
      implicit none
 
      integer :: LWR1, LWR2, LWR3, LWR4, LWI, LWC, NS, NE, NC, NREA, NV,&
&                 NVIB, NELEQ, IMOD, NMAX
      real*8, dimension(:), allocatable :: WR1, WR2, WR3, WR4
      integer, dimension(:), allocatable :: WI
      character, dimension(:), allocatable :: WC 
      CHARACTER(10) MIXTURE, REACTION, TRANSFER
      CHARACTER(100) PATH

      end module mutvar


      program cfd
      use mutvar
      implicit none 
      integer :: i 

!     Files and path
!      -mixture file located in directory /data/mixture
      MIXTURE = 'air5'
!      -reaction file located in directory /data/chemistry/gasreact
!       for nonreacting mixture, use 'empty'
      REACTION = 'air5'
!      -energy transfer file located in directory /data/transfer
      TRANSFER = 'empty'
!      -path for the library location
      PATH = '..'

!     Pseudo-dynamic allocation (mandatory routine)
      CALL LENGTH (PATH, MIXTURE, REACTION, TRANSFER, LWR1, LWR2, LWR3, &
     &             LWR4, LWI, LWC, NS, NE, NC, NREA, NV, NMAX, NVIB,    &
     &             NELEQ)


      end program cfd
