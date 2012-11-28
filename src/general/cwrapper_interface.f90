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
    end interface

end module mutationpp
