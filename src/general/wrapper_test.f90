program main
    use mutationpp
    implicit none
    character(len=10) :: mixture
    character(len=12) :: species
    character(len=10) :: element
    integer :: i, j, ne, ns, nr
    
    real :: T, P, rho, n, cp, cv, h, mw, e
    real, dimension(:), allocatable :: element_x
    real, dimension(:), allocatable :: species_x
    real, dimension(:), allocatable :: species_y
    real, dimension(:), allocatable :: conc
    real, dimension(:), allocatable :: wdot
    real, dimension(:), allocatable :: mwi
    
    mixture = "air11"
    
    call mpp_initialize(mixture)
    ne = mpp_nelements()
    ns = mpp_nspecies()
    nr = mpp_nreactions()
    
    write(*,*) "ne: ", ne, "ns: ", ns, "nr: ", nr
    
    allocate(element_x(ne))
    allocate(species_x(ns))
    allocate(species_y(ns))
    allocate(conc(ns))
    allocate(wdot(ns))
    allocate(mwi(ns))
    
    element = "N";  element_x(mpp_element_index(element)) = 0.79;
    element = "O";  element_x(mpp_element_index(element)) = 0.21;
    element = "e-"; element_x(mpp_element_index(element)) = 0.0;
    
    P = 101325.0
    call mpp_species_mw(mwi)
    
    write(*,'(A12)',advance='no') "T(K)        "
    do j = 1,ns
        call mpp_species_name(j, species)
        write(*,'(A12)',advance='no') "X_"//species
    end do
    write(*,'(A12)',advance='no') "P(Pa)       "
    write(*,'(A12)',advance='no') "rho(kg/m^3) "
    write(*,'(A12)',advance='no') "n(1/m^3)    "
    write(*,'(A12)',advance='no') "Mw(kg/mol)  "
    write(*,'(A12)',advance='no') "Cp(J/kg-K)  "
    write(*,'(A12)',advance='no') "Cv(J/kg-K)  "
    write(*,'(A12)',advance='no') "h(J/kg)     "
    write(*,'(A12)',advance='no') "e(J/kg)     "
    write(*,*)
    
    ! Loop over temperature and compute equilibrium properties
    do i = 1,295
        T = dble(i-1)*50.0 + 300.0
        
        call mpp_equilibrate_mole(T, P, element_x, species_x)
            
        write(*,'(E12.4)',advance='no') T
        do j = 1,ns
            write(*,'(E12.4)',advance='no') species_x(j)
        end do
        
        call mpp_convert_x_to_y(species_x, species_y)        
        call mpp_number_density(T, P, n)
        call mpp_density(T, P, species_x, rho)
        call mpp_mixture_mw_mole(species_x, mw)
        call mpp_mixture_frozen_cp_mass(T, species_y, cp)
        call mpp_mixture_frozen_cv_mass(T, species_y, cv)
        call mpp_mixture_h_mass(T, species_y, h)
        call mpp_mixture_e_mass(T, rho, species_y, e)
        
        !conc = species_y * rho / mwi
        !call mpp_net_production_rates(T, conc, wdot)
        
        
        write(*,'(8E12.4)',advance='no') P, rho, n, mw, cp, cv, h, e
        write(*,*)
    end do
    
    ! Clean up the memory stored in the mutation++ library
    call mpp_destroy()

end program main
