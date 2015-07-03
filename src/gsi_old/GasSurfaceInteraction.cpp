#include "GasSurfaceInteraction.h"
#include "Utilities.h"
#include "Vector.h"

using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities;

namespace Mutation {
    namespace gsi {
      
//==============================================================================
      
GasSurfaceInteraction::GasSurfaceInteraction( Mutation::Thermodynamics::Thermodynamics& thermo, Mutation::Transport::Transport& transport, std::string gsi_catalysis_mechanism_file )
                                            : m_thermo( thermo ),
                                              m_transport( transport ),
                                              m_wall_state_set( false ),
                                              mp_wall_solver( NULL ){

    if ( gsi_catalysis_mechanism_file == "none" ){ return; }
    
    gsi_catalysis_mechanism_file =
    getEnvironmentVariable( "MPP_DATA_DIRECTORY" ) + "/gsi/" +
    gsi_catalysis_mechanism_file + ".xml";

    IO::XmlDocument doc(gsi_catalysis_mechanism_file);
    IO::XmlElement root = doc.root();

    if ( root.tag() != "gsi_mechanism" ) {
        std::cout << "Root element in gsi_mechanism file " << gsi_catalysis_mechanism_file
                  << " is not of 'gsi_mechanism' type!" << std::endl;
        exit(1);
    }

    root.getAttribute( "category", m_category, m_category );
    if (m_category != "catalysis") {
        std::cout << "The Gas-Surface interaction category in the "
           << gsi_catalysis_mechanism_file << " file is not of 'catalysis' type!" << std::endl;
        exit(1);
    }

    root.getAttribute( "model", m_catalytic_model, m_catalytic_model );

//    mp_catalysis_rate = Config::Factory<CatalysisRateManager>::create( m_catalytic_model, thermo, mp_surface_properties, v_catalytic_reactions );

    if ( m_catalytic_model == "gamma" ){
        mp_surface_properties = NULL;

        IO::XmlElement::const_iterator iter = root.begin();
        for ( ; iter != root.end(); ++iter) {
            if (iter->tag() == "reaction")
                addCatalyticReaction(CatalysisReaction(*iter, thermo, m_catalytic_model, mp_surface_properties));
        }

        // Setup the rate manager
        mp_catalysis_rates = new CatalysisRateManagerGamma( m_thermo, v_catalytic_reactions );
        mp_wall_state = new WallState( m_thermo );
    }
    else if (m_catalytic_model == "finite_rate_chemistry"){
        root.getAttribute("surface", m_surface, m_surface);
        mp_surface_properties = new CatalysisSurfaceProperties(m_surface, m_thermo);

        IO::XmlElement::const_iterator iter = root.begin();
        for ( ; iter!=root.end(); ++iter) {
            if (iter->tag() == "reaction")
                addCatalyticReaction(CatalysisReaction(*iter, thermo, m_catalytic_model, mp_surface_properties));
        }
        mp_catalysis_rates = new CatalysisRateManagerFiniteRateChemistry( thermo, mp_surface_properties, v_catalytic_reactions );
        mp_wall_state = new WallStateFRC( *mp_surface_properties, m_thermo );
    }
    else{
        std::cout << "The catalytic model " << m_catalytic_model
        << " has not been implemented yet!";
        exit(1);
    }
        
    //Closing the GSI reaction mechanisms
    closeGSIReactions(true);
}

//==============================================================================

GasSurfaceInteraction::~GasSurfaceInteraction(){

    if ( mp_catalysis_rates != NULL ) delete mp_catalysis_rates;
    if ( mp_wall_state != NULL ) delete mp_wall_state;
    if ( mp_surface_properties != NULL ) delete mp_surface_properties;
    if ( mp_wall_solver != NULL ) delete mp_wall_solver;
}


//==============================================================================

void GasSurfaceInteraction::netGSIProductionRates(double * const p_wdot){

    if (!m_wall_state_set){
       std::cerr << "Gas-Surface Interaction Error: The state of the wall should be set before calling for the gas-surface interaction production rates." << endl;
       std::cerr << "Gas-Surface Interaction Error: Call the setWallState function with the partial densities and the wall temperatures" << endl;
       exit(1);
    }
  
    for (int i = 0; i < m_thermo.nSpecies(); ++i)
        p_wdot[i] = 0.;
    
    if (nCatalyticReactions() != 0) {
        mp_catalysis_rates->computeRate(*mp_wall_state, p_wdot);
    }

}

//==============================================================================

void GasSurfaceInteraction::solveWallMassBalance(double * const p_rhoi_wall, const double * const p_rhoi_edge, const double * const p_Twall, const double& dx_gradient_distance_wall_edge){

    int m_ns = m_thermo.nSpecies();

    if (mp_wall_solver == NULL){
        mp_wall_solver = new WallSolver(m_thermo, m_transport, mp_wall_state, mp_catalysis_rates);
    }

    mp_wall_solver->setDataSolver(p_rhoi_edge, p_Twall, dx_gradient_distance_wall_edge); 
    Mutation::Numerics::RealVector v_solution(m_thermo.nSpecies());
    for (int i_ns = 0 ; i_ns < m_ns ; i_ns++){
        v_solution(i_ns) = p_rhoi_wall[i_ns];
    }

    // This is something that can be done internally I guess! Maybe method of mp_wall_solver-> solveUsingPartialDensities(double * const p_rhoi_wall);
    mp_wall_solver->initializeSolutionwithPartialDensities( v_solution );

    mp_wall_solver->retrieveInitialMolarFractions( v_solution );
    v_solution = mp_wall_solver->solve( v_solution );

    mp_wall_solver->returnSolutioninPartialDensities( v_solution );

    for (int i_ns = 0; i_ns < m_ns ; i_ns++){
        p_rhoi_wall[i_ns] = v_solution(i_ns);
    }

}

//==============================================================================

void GasSurfaceInteraction::addCatalyticReaction(const CatalysisReaction& catalytic_reaction){

    v_catalytic_reactions.push_back(catalytic_reaction);

}

//==============================================================================

void GasSurfaceInteraction::closeGSIReactions(const bool validate_gsi_mechanism)
{
   const size_t ns = m_thermo.nSpecies();
   
   // GSI Mechanism validation
   if (validate_gsi_mechanism) {
   
       bool is_valid = true;
       
/**
 * @todo For this in the GSIReaction the Big Three should be implemented
 *  http://stackoverflow.com/questions/4172722/what-is-the-rule-of-three
 */
       
       //Check for duplicate gsi reactions
       //for (size_t i = 0; i < nCatalyticReactions()-1; i++) {
       //    for (size_t j = i+1; j < nGSIReactions(); ++j) {
       //        if (m_gsi_reactions[i] == m_gsi_reactions)
       //    }
       //}
       
       // Check for elemental mass and charge conservation
//       for (size_t i = 0; i < nGSIReactions(); i++) {
//           if (!m_gsi_reactions[i].conservesChargeAndMass()){
//              std::cerr << "GSI reaction " << i+1 << " \"" << m_gsi_reactions[i].formula()
//                        << " \" does not conserve charge or mass." << endl;
//              is_valid = false;
//           }
//       }
       
//       // Check for the existance of active sites
//       for (size_t i = 0; i < nGSIReactions(); i++) {
//           if (!m_gsi_reactions[i].hasActiveSites()){
//              std::cerr << "GSI reaction " << i+1 << " \"" << m_gsi_reactions[i].formula()
//                        << " \" has no active sites or surface species." << endl;
//              is_valid = false;
//           }
//       }
       
       // Check for the conservation of active sites
//       for (size_t i = 0; i < nGSIReactions(); i++) {
//           if (!m_gsi_reactions[i].conservesActiveSites()){
//              std::cerr << "GSI reaction " << i+1 << " \"" << m_gsi_reactions[i].formula()
//                        << " \" does not conserve active sites" << endl;
//              is_valid = false;
//           }
//       }
       
       if (!is_valid) {
           cout << "Validation check of the gas surface interaction reactions failed!" << endl;
           exit(1);
       }
       
   }
   
}

//==============================================================================
      
void GasSurfaceInteraction::setWallState(const double* const p_rhoi, const double* const p_rhoie)
{
    mp_wall_state->setWallState(p_rhoi, p_rhoie);
    m_wall_state_set = true;
}
      
    }// namespace GSI
} // namespace Mutation

//==============================================================================
