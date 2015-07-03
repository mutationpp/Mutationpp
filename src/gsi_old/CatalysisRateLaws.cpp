#include <cassert>
#include <stdlib.h>
#include "math.h"

#include "CatalysisRateLaws.h"
#include "Constants.h"

using namespace Mutation::Utilities;

namespace Mutation{
    namespace gsi{

inline double CatalysisRateLaw::computeAverageThermalSpeed( const int& l_sp, const double* const lp_Twall ) const {
    return sqrt( ( 8 * Mutation::RU * lp_Twall[0] ) / ( Mutation::PI * m_thermo.speciesMw( l_sp ) ) ) ;
}

inline double CatalysisRateLaw::computeAverage2DThermalSpeed(const int&  l_sp, const double* const p_Twall) const {
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
                                  v_reactants(reactants)
{
    assert( node.tag() == "gamma_const" );

    int l_diff_reac = 1;
    for ( int i_reac = 0; i_reac < reactants.size() - 1; ++i_reac ){
        if (reactants[i_reac] != reactants[i_reac + 1]) l_diff_reac++;
    }

    /** @todo give an error somewhere if you have more than 2 reactants. **/

    std::vector<std::string> tokens;
    String::tokenize(node.text(), tokens, ":, ");
    int index_ref = -2;
    for (int i = 0; i < tokens.size(); i+=2) {
        int index = thermo.speciesIndex(tokens[i]);
        if(index_ref > index){
            v_gamma.insert(v_gamma.begin(), atof(tokens[i+1].c_str()));
        } else {
            v_gamma.insert(v_gamma.end(), atof(tokens[i+1].c_str()));
        }
        index_ref = index;
    } /** @todo: MOST UGLY THING IN MY LIFE. DEBUG IT WITH THE EXTREME CASES **/

    if (l_diff_reac != v_gamma.size()){
        std::cerr << "ERROR @ GammaModelConst @ CatalysisRateLaws.cpp. #gammas should be = #different reactants." << std::endl; /** @todo cout better error */
        exit(1);
    }
    /**
     * @todo Add Controls here? Or in the rate manager
     * @todo What should I add here?
     */

    v_output_impinging_flux.resize(v_gamma.size());
    v_impinging_flux_per_stoichiometric_coefficient.resize(v_gamma.size());

}

//==============================================================================

double GammaModelConst::forwardReactionRate(const double * const lp_rhoi, const double* const lp_Twall) const {

    l_index_v_reactants = 0;

    for( int i_gammas = 0 ; i_gammas < v_gamma.size() ; i_gammas++ ){

        getSpeciesIndexandStoichiometricCoefficient( l_index_v_reactants, l_species_index, l_stoichiometric_coefficient );
        v_output_impinging_flux[ i_gammas ] = computeWallImpingingMassFlux( l_species_index, lp_rhoi, lp_Twall );

        v_impinging_flux_per_stoichiometric_coefficient[ i_gammas ] = v_output_impinging_flux[ i_gammas ] / l_stoichiometric_coefficient;
        v_output_impinging_flux[ i_gammas ] = v_impinging_flux_per_stoichiometric_coefficient[ i_gammas ] * v_gamma[ i_gammas ];

        l_index_v_reactants += l_stoichiometric_coefficient;

    }

    return getLimitingImpingingMassFlux( );

}

//==============================================================================

inline void GammaModelConst::getSpeciesIndexandStoichiometricCoefficient( int l_index_v_reactants, int& l_species_index, int& l_stoichiometric_coefficient ) const {

    l_species_index = v_reactants[l_index_v_reactants];
    l_stoichiometric_coefficient = 1;
    l_index_v_reactants++;
    
    while( l_index_v_reactants < v_reactants.size() ){
        if ( l_species_index != v_reactants[ l_index_v_reactants ] ) {
            break;
        } 
        l_stoichiometric_coefficient++;
        l_index_v_reactants++;
    }

}

//==============================================================================

inline double GammaModelConst::computeWallImpingingMassFlux( const int& l_sp, const double* const lp_rhoi, const double* const lp_Twall ) const {

    return computeAverageThermalSpeed( l_sp, lp_Twall ) / ( 4.E0 ) * lp_rhoi[ l_sp ] / m_thermo.speciesMw( l_sp );

}

//==============================================================================

inline double GammaModelConst::getLimitingImpingingMassFlux() const {

    return v_output_impinging_flux[ min_element( v_impinging_flux_per_stoichiometric_coefficient.begin(), v_impinging_flux_per_stoichiometric_coefficient.end()) - v_impinging_flux_per_stoichiometric_coefficient.begin() ];

}

//==============================================================================
//============================ Physisorption ===================================
//==============================================================================
      
Physisorption::Physisorption( const Mutation::Utilities::IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics& thermo )
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

LHRecombination::LHRecombination( const Utilities::IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics& thermo )
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
