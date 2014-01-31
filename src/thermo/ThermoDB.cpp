
#include "ThermoDB.h"
#include "Utilities.h"

#include <algorithm>
#include <cassert>

namespace Mutation {
    namespace Thermodynamics {

//==============================================================================

ThermoDB::ThermoDB(double sst, double ssp)
    : m_sst(sst), m_ssp(ssp)
{
    assert(m_sst >= 0.0);
    assert(m_ssp >  0.0);
}

//==============================================================================

bool ThermoDB::load(const SpeciesListDescriptor& descriptor)
{
    // It is possible that this isn't the first call to load so just make sure
    // species and elements are cleared
    m_species.clear();
    m_elements.clear();
    
    // First we need to load the entire element database for use in constructing
    // our species list
    Utilities::IO::XmlDocument element_doc(
        Utilities::getEnvironmentVariable("MPP_DATA_DIRECTORY") +
        "/thermo/elements.xml");
    Utilities::IO::XmlElement::const_iterator element_iter =
        element_doc.root().begin();

    std::vector<Element> elements;
    for ( ; element_iter != element_doc.root().end(); ++element_iter)
        elements.push_back(Element(*element_iter));
    
    // Load all possible species from the concrete database type
    std::list<Species> species_list;
    loadAvailableSpecies(species_list, elements);
    
    // Check for duplicate species in the database
    std::list<Species>::iterator iter1 = species_list.begin();
    std::list<Species>::iterator iter2;
    
    while (iter1 != species_list.end()) {
        (iter2 = iter1)++;
        while (iter2 != species_list.end()) {
            if (iter1->name() == iter2->name()) {
                std::cout << "Warning, species \"" << iter1->name()
                          << "\" is defined more than once in the thermodynamic"
                          << " database.  I will ignore recurrences..."
                          << std::endl;
                iter2 = species_list.erase(iter2);
                while (iter2 != species_list.end()) {
                    if (iter1->name() == iter2->name())
                        iter2 = species_list.erase(iter2);
                    else
                        iter2++;
                }
            } else
                iter2++;
        }
        iter1++;
    }
    
    // Use the species list descriptor to remove unwanted species
    iter1 = species_list.begin();
    while (iter1 != species_list.end()) {
        if (descriptor.matches(*iter1))
            iter1++;
        else
            iter1 = species_list.erase(iter1);
    }
    
//    std::cout << "species before ordering" << std::endl;
//    iter = species_list.begin();
//    while (iter != species_list.end()) {
//        std::cout << iter->name() << std::endl;
//        iter++;
//    }
    
    // Now we have all of the species that we want but possibly in the wrong
    // order so use the descriptor to tell us the correct order
    std::vector<std::string> missing;
    descriptor.order(species_list, m_species, missing);
    
    if (missing.size() > 0) {
        std::cout << "Missing the following species!" << std::endl;
        for (int i = 0; i < missing.size(); ++i)
            std::cout << "  " << missing[i] << std::endl;
        return false;
    }
    
//    std::cout << "species after ordering" << std::endl;
//    for (int i = 0; i < m_species.size(); ++i)
//        std::cout << m_species[i].name() << std::endl;
    
    // Finally fill our elements vector with only the elements that are required
    // for the species list
    std::set<std::string> element_names;
    for (int i = 0; i < m_species.size(); ++i) {
        Species::StoichList::const_iterator iter =
            m_species[i].stoichiometry().begin();
        for ( ; iter != m_species[i].stoichiometry().end(); ++iter)
            element_names.insert(iter->first);
    }
    
    for (int i = 0; i < elements.size(); ++i)
        if (element_names.count(elements[i].name()) > 0)
            m_elements.push_back(elements[i]);
    
    // Tell the concrete class to load the necessary thermodynamic data based on
    // the final species list
    loadThermodynamicData();
    return true;
}

//==============================================================================

// Used in ThermoDB::cv
struct MinusOne {
    void operator () (double& x) const { x -= 1.0; }
} MinusOne;

void ThermoDB::cv(
    double Th, double Te, double Tr, double Tv, double Tel, 
    double* const cv = NULL, double* const cvt = NULL, double* const cvr = NULL,
    double* const cvv = NULL, double* const cvel = NULL)
{
    // Compute Cp/Ru
    cp(Th, Te, Tr, Tv, Tel, cv, cvt, cvr, cvv, cvel);
    
    const size_t ns = m_species.size();

    // Cv/Ru = Cp/Ru - 1
    if (cv   != NULL) std::for_each(  cv+0,   cv+ns, MinusOne);
    if (cvt  != NULL) std::for_each( cvt+0,  cvt+ns, MinusOne);
    if (cvr  != NULL) std::for_each( cvr+0,  cvr+ns, MinusOne);
    if (cvv  != NULL) std::for_each( cvv+0,  cvv+ns, MinusOne);
    if (cvel != NULL) std::for_each(cvel+0, cvel+ns, MinusOne);
}

//==============================================================================

    } // namespace Thermodynamics
} // namespace Mutation

