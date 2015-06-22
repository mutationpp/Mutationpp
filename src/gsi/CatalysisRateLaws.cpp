#include <cassert>
#include <stdlib.h>
#include "math.h"

#include "CatalysisRateLaws.h"
#include "Constants.h"

using namespace Mutation::Utilities;

namespace Mutation{
    namespace gsi{

inline double CatalysisRateLaw::AverageThermalSpeed(const int& i_species, const double* const p_rhoie) const {
    return sqrt( (8 * RU * p_rhoie[0]) / (PI * m_thermo.speciesMw(i_species))) ;
}

inline double CatalysisRateLaw::Average2DThermalSpeed(const int&  i_species, const double* const p_rhoie) const {
//    for (int i = 0; i < nSpecies(); i++){
//        double c = sqrt(PI * RU /* * mp_Tw */ / (PI * MolarMass));
//    }
    
    return 0.E0; //sqrt(RU /* * mp_Tw */ / (PI * m_thermo.speciesMw(i_species))) ;
}

//==============================================================================
//============================ GammaModelConst =================================
//==============================================================================
      
/** SEEMS TO BE WORKING. TO BE REIMPLEMENTED BECAUSE IT IS NOT VERY GOOD **/

GammaModelConst::GammaModelConst(const Mutation::Utilities::IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics& thermo, std::vector<int> reactants)
                                : CatalysisRateLaw(thermo),
                                  m_reactants(reactants)
{
    assert( node.tag() == "gamma_const" );

    int l_diff_reac = 1;
    for (int i_reac = 0; i_reac < reactants.size()-1; ++i_reac ){
        if (reactants[i_reac] != reactants[i_reac + 1]) l_diff_reac++;
    }

    /** @todo give an error somewhere if you have more than 2 reactants. **/

    std::vector<std::string> tokens;
    String::tokenize(node.text(), tokens, ":, ");
    int index_ref = -2;
    for (int i = 0; i < tokens.size(); i+=2) {
        int index = thermo.speciesIndex(tokens[i]);
        if(index_ref > index){
            m_gamma.insert(m_gamma.begin(), atof(tokens[i+1].c_str()));
        } else {
            m_gamma.insert(m_gamma.end(), atof(tokens[i+1].c_str()));
        }
        index_ref = index;
    } /** @todo: MOST UGLY THING IN MY LIFE. DEBUG IT WITH THE EXTREME CASES **/

    if (l_diff_reac != m_gamma.size()){
        std::cerr << "ERROR @ GammaModelConst @ CatalysisRateLaws.cpp. #gammas should be = #different reactants." << std::endl; /** @todo cout better error */
        exit(1);
    }
    /**
     * @todo Add Controls here? Or in the rate manager
     * @todo What should I add here?
     */

}

//==============================================================================

double GammaModelConst::forwardReactionRate(const double * const p_rhoi, const double* const p_rhoie) const {
    /** @todo Do it with pointer arithmetics if it called this way... */

    double l_test_flux;
    double l_min_output_flux, l_test_output_flux;
    int index = m_reactants[0];
    int increase_index = 1;
    while (index == m_reactants[increase_index] && increase_index < m_reactants.size() ) increase_index++;

    l_min_output_flux = AverageThermalSpeed(index, p_rhoie) / (4.E0) * p_rhoi[index]
                                                                     / m_thermo.speciesMw(index) ;
    double l_minim_flux = l_min_output_flux / increase_index;
    l_min_output_flux = l_minim_flux * m_gamma[0];

    for(int i_diff_species = 1; i_diff_species < m_gamma.size(); ++i_diff_species){
        int i_increase = 1;
        index = m_reactants[increase_index];

        while (index == m_reactants[increase_index + 1] && increase_index < m_reactants.size() ) {increase_index++ ; i_increase++;}

        l_test_output_flux = AverageThermalSpeed(index, p_rhoie) / (4.E0 ) * p_rhoi[index]
                                                                     / m_thermo.speciesMw(index) ;

        l_test_flux = l_test_output_flux / i_increase ;

        l_test_output_flux = l_test_flux * m_gamma[i_diff_species];
        if ( l_minim_flux > l_test_flux ){
            l_min_output_flux = l_test_output_flux;
            l_minim_flux = l_test_flux;
        }
    }

    return l_min_output_flux; // It can be done with less variables. DEBUG!
}

//==============================================================================
// 
// const double* GammaModelConst::getgammaCoefficient() const
// {
//     return &m_gamma[0];
// }
// 
//==============================================================================
/*
GammaModelT::GammaModelT(const Mutation::Utilities::IO::XmlElement& node){assert( node.tag() == "gamma_T" );}
GammaModelTP::GammaModelTP(const Mutation::Utilities::IO::XmlElement& node){assert( node.tag() == "gamma_const" );}
    */
//==============================================================================
      
//==============================================================================
//============================ Physisorption ===================================
//==============================================================================
      
Physisorption::Physisorption(const Mutation::Utilities::IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics& thermo)
    : CatalysisRateLaw(thermo),
      m_S_coef(1.E0),
      m_beta(0.E0),
      m_E_act(0.E0)
{
    assert( node.tag() == "physisorption" );
     
    node.getAttribute("sticking_coef", m_S_coef, m_S_coef);
    node.getAttribute("b", m_beta, m_beta);
    node.getAttribute("E_activation", m_E_act, "In a physisorption reaction the activation energy should be provided");
     
     /**
      * @todo Error parser if sb tries to add another kind of input!
      */
}
      
double Physisorption::forwardReactionRate(const double* const p_rhoi, const double* const p_rhoie) const
{
     return( m_S_coef  * pow(p_rhoie[0], m_beta) * exp(-(m_E_act)/(RU * p_rhoie[0])));

     /** @todo Multiply by the (1 - phi) */

}

//==============================================================================
//=========================== ThermalDesoption =================================
//==============================================================================
      
ThermalDesorption::ThermalDesorption(const Utilities::IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics& thermo)
    : CatalysisRateLaw(thermo),
      m_steric_factor(1.E0),
      m_beta(0.E0),
      m_vib_perp(1.E0),
      m_E_des(0.E0)
{
    assert( node.tag() == "thermal_desorption");
    
    node.getAttribute("steric_factor", m_steric_factor, m_steric_factor);
    node.getAttribute("b", m_beta, m_beta);
    node.getAttribute("freq_vibration_perp", m_vib_perp, "In a thermal desorption reaction the vibrational frequency perpendicular to the surface should be provide");
    node.getAttribute("E_desorption", m_E_des, "In a thermal desorption reaction the activation energy should be provide");
    
}

double ThermalDesorption::forwardReactionRate(const double* const p_rhoi, const double* const p_rhoie) const
{
    return( m_steric_factor * m_vib_perp * pow(p_rhoie[0], m_beta) * exp(-(m_E_des)/(RU * p_rhoie[0] )));
}


//==============================================================================
//============================ Chemisorption ===================================
//==============================================================================

Chemisorption::Chemisorption(const Utilities::IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics& thermo)
    : CatalysisRateLaw(thermo),
      m_S_coef(1.E0),
      m_beta(0.E0),
      m_E_act(0.E0)
{
    assert( node.tag() == "chemisorption" );
     
    node.getAttribute("sticking_coef", m_S_coef, m_S_coef);
    node.getAttribute("b", m_beta, m_beta);
    node.getAttribute("E_activation", m_E_act, "In a chemisorption reaction the activation energy should be provided");

}

double Chemisorption::forwardReactionRate(const double* const p_rhoi, const double* const p_rhoie) const
{
    return( m_S_coef  * pow(p_rhoie[0], m_beta) * exp(-(m_E_act)/(RU  * p_rhoie[0])));
}

//==============================================================================
//=========================== E-R Recombination ================================
//==============================================================================

ERRecombination::ERRecombination(const Utilities::IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics& thermo)
    : CatalysisRateLaw(thermo),
      m_steric_factor(1.E0),
      m_beta(0.E0),
      m_E_act_recomb(0.E0)
{
    assert( node.tag() == "E-R" );
    
    node.getAttribute("steric_factor", m_steric_factor, m_steric_factor);
    node.getAttribute("b", m_beta, m_beta);
    node.getAttribute("E_activ_recomb", m_E_act_recomb, "In an Eley-Rideal reaction the activation energy should be provided");
}

double ERRecombination::forwardReactionRate(const double* const p_rhoi, const double* const p_rhoie) const
{
    return( m_steric_factor * pow(p_rhoie[0], m_beta) * exp(-(m_E_act_recomb)/(RU  * p_rhoie[0] )));
}

//==============================================================================
//===================== Physisorption to Chemisorption =========================
//==============================================================================

PhysisorptiontoChemisorption::PhysisorptiontoChemisorption(const Utilities::IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics& thermo)
    : CatalysisRateLaw(thermo),
      m_steric_factor(1.0E0),
      m_vib_par(1.0E0),
      m_E_diff(0.E0)
{
    assert( node.tag() == "phys_to_chem" );
  
    node.getAttribute("steric_factor", m_steric_factor, m_steric_factor);
    node.getAttribute("freq_vibration_parallel", m_vib_par, "In a diffusion reaction from a physisorption to a chemisorption sitet he vibrational frequency parallel to the surface should be provide");
    node.getAttribute("E_diffusion_barrier", m_E_diff, "In a diffusion reaction from a physisorption to a chemisorption site the energy of the diffusion barrier should be provide");
}


//==============================================================================
//========================= L-H Recombination ==================================
//==============================================================================

LHRecombination::LHRecombination(const Utilities::IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics& thermo)
    : CatalysisRateLaw(thermo) //CHECK WHICH ONE NEEDS THERMO COPY
{
    assert( node.tag() == "L-H" );
    
    /**
     * @todo L-H for guerra needs no data. This is not the case for Marschall's approach
     */
    
}

//==============================================================================
//==============================================================================
//==============================================================================

    } // namespace gsi
} // namespace Mutation
