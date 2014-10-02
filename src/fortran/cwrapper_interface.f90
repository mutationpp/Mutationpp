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
        
        real(kind=8) function mpp_mixture_frozen_cv_mass();
        end function
        
        real(kind=8) function mpp_mixture_frozen_gamma();
        end function
        
        real(kind=8) function mpp_mixture_frozen_sound_speed();
        end function
        
        real(kind=8) function mpp_mixture_h_mass()
        end function
        
        real(kind=8) function mpp_mixture_e_mass()
        end function
        
        integer function mpp_ncollision_pairs()
        end function
        
        real(kind=8) function mpp_viscosity()
        end function

        real(kind=8) function mpp_frozen_thermal_conductivity()
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
