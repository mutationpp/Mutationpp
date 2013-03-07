module mutationpp

    interface
        integer function mpp_nelements()
        end function
    
        integer function mpp_nspecies()
        end function
        
        integer function mpp_nreactions()
        end function
        
        integer function mpp_element_index(element)
            character(len=*) :: element
        end function
        
        integer function mpp_species_index(species)
            character(len=*) :: species
        end function
        
        real function mpp_number_density()
        end function
        
        real function mpp_density()
        end function
        
        real function mpp_pressure()
        end function
        
        real function mpp_mixture_mw()
        end function
        
        real function mpp_mixture_frozen_cp_mass()
        end function
        
        real function mpp_mixture_frozen_cv_mass();
        end function
        
        real function mpp_mixture_h_mass()
        end function
        
        real function mpp_mixture_e_mass()
        end function
    end interface

end module mutationpp
