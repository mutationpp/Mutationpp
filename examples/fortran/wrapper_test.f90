!!
!! Copyright 2014 von Karman Institute for Fluid Dynamics (VKI)
!!
!! This file is part of MUlticomponent Thermodynamic And Transport
!! properties for IONized gases in C++ (Mutation++) software package.
!!
!! Mutation++ is free software: you can redistribute it and/or modify
!! it under the terms of the GNU Lesser General Public License as
!! published by the Free Software Foundation, either version 3 of the
!! License, or (at your option) any later version.
!!
!! Mutation++ is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU Lesser General Public License for more details.
!!
!! You should have received a copy of the GNU Lesser General Public
!! License along with Mutation++.  If not, see
!! <http://www.gnu.org/licenses/>.
!!

!!
!> @file wrapper_test.f90
!!
!! Example Fortran90 program which makes use of the \ref FortranWrapper to 
!! mutation++.  Computes equilibrium properties for the air11 mixture.
!!
program main
    use mutationpp
    implicit none
    character(len=10) :: mixture, state_model
    character(len=12) :: species
    integer :: i, j, ne, ns, nr, var
    
    real(8) :: T, P, rho, n, cp, cv, h, mw, e
    real(8), dimension(:), allocatable :: species_x, wdot, species_y
    real(8), dimension(:), allocatable :: mwi
    
    mixture     = "air_11"
    state_model = "Equil"
    var = 1
    
    call mpp_initialize(mixture, state_model)
    ne = mpp_nelements()
    ns = mpp_nspecies()
    nr = mpp_nreactions()
    
    write(*,*) "ne: ", ne, "ns: ", ns, "nr: ", nr
    
    allocate(species_x(ns))
    allocate(species_y(ns))
    allocate(wdot(ns))
    allocate(mwi(ns))
    
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
        
        call mpp_set_state(T, P, var)
        call mpp_x(species_x)
            
        write(*,'(E12.4)',advance='no') T
        do j = 1,ns
            write(*,'(E12.4)',advance='no') species_x(j)
        end do
             
        n   = mpp_number_density()
        rho = mpp_density()
        mw  = mpp_mixture_mw()
        cp  = mpp_mixture_frozen_cp_mass()
        cv  = mpp_mixture_frozen_cv_mass()
        h   = mpp_mixture_h_mass()
        e   = mpp_mixture_e_mass()
        
        write(*,'(8E12.4)') P, rho, n, mw, cp, cv, h, e
    end do
    
    ! Clean up the memory stored in the mutation++ library
    call mpp_destroy()

end program main
