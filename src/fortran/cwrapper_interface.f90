!!> @file cwrapper_interface.f90

!!
!! Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
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
!! Provides an explicit interface for all Mutation++ functions which have a 
!! non-void return value which is need to compile user fortran codes.
!!
module mutationpp

    interface
        integer function mpp_nelements()
        end function
    
        integer function mpp_nspecies()
        end function
        
        integer function mpp_nreactions()
        end function
        
        integer function mpp_n_mass_eqns()
        end function
        
        integer function mpp_n_energy_eqns()
        end function
        
        integer function mpp_element_index(element)
            character(len=*) :: element
        end function
        
        integer function mpp_species_index(species)
            character(len=*) :: species
        end function
        
        real(kind=8) function mpp_number_density()
        end function
        
        real(kind=8) function mpp_density()
        end function

        !real(kind=8) function mpp_density_tpx(T, P, X)
        !    real(kind=8), intent(in) :: T
        !    real(kind=8), intent(in) :: P
        !    real(kind=8), dimension(:) :: X
        !end function

        real(kind=8) function mpp_pressure()
        end function
        
        real(kind=8) function mpp_mixture_mw()
        end function

        real(kind=8) function mpp_mixture_t()
        end function
        
        real(kind=8) function mpp_mixture_frozen_cp_mass()
        end function
        
        real(kind=8) function mpp_mixture_frozen_cv_mass()
        end function
        
        real(kind=8) function mpp_mixture_frozen_gamma()
        end function
        
        real(kind=8) function mpp_mixture_frozen_sound_speed()
        end function
        
        real(kind=8) function mpp_mixture_h_mass()
        end function
        
        real(kind=8) function mpp_mixture_e_mass()
        end function
        
        integer function mpp_ncollision_pairs()
        end function
        
        real(kind=8) function mpp_viscosity()
        end function

        real(kind=8) function mpp_equilibrium_thermal_conductivity()
        end function
        
        real(kind=8) function mpp_heavy_thermal_conductivity()
        end function

        real(kind=8) function mpp_electron_thermal_conductivity()
        end function

        real(kind=8) function mpp_internal_thermal_conductivity()
        end function

        real(kind=8) function mpp_reactive_thermal_conductivity()
        end function
   
        real(kind=8) function mpp_sigma()
        end function
        
    end interface

end module mutationpp
