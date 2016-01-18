#include "GasSurfaceInteraction.h"
#include "Utilities.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//======================================================================================

GasSurfaceInteraction::GasSurfaceInteraction( Mutation::Thermodynamics::Thermodynamics& l_thermo, Mutation::Transport::Transport& l_transport, std::string l_gsi_input_file )
                                            : m_thermo( l_thermo ),
                                              m_transport( l_transport ),
                                              mp_surf_solver( NULL ),
                                              mp_surf_props( NULL ),
                                              mp_wall_state( NULL ),
                                              v_mass_prod_rate( m_thermo.nSpecies() ),
                                              v_mole_frac_edge( m_thermo.nSpecies() ) {

    using namespace Mutation::Utilities::IO;

    if ( l_gsi_input_file == "none" ){ return; }

    locateGSIInputFile( l_gsi_input_file );
    
    XmlDocument l_xml_doc( l_gsi_input_file );
    XmlElement l_root_element = l_xml_doc.root();

    errorWrongTypeofGSIFile( l_root_element.tag() );

    l_root_element.getAttribute( "gsi_mechanism", m_gsi_mechanism, "none" );

    XmlElement::const_iterator xml_position_surf_props = l_root_element.findTag("surface_properties");
    XmlElement::const_iterator xml_position_diff_model = l_root_element.findTag("diffusion_model");
    XmlElement::const_iterator xml_position_prod_terms = l_root_element.findTag("production_terms");

    DataSurfaceProperties l_data_surface_properties = { m_thermo, *xml_position_surf_props };
    mp_surf_props = Mutation::Utilities::Config::Factory<SurfaceProperties>::create( m_gsi_mechanism, l_data_surface_properties );
    mp_wall_state = new WallState( m_thermo, *mp_surf_props );

    // Creating the SurfaceBalanceSolver class
    DataSurfaceBalanceSolver l_data_surface_balance_solver = { m_thermo, m_transport, m_gsi_mechanism, *xml_position_diff_model, 
                                                               *xml_position_prod_terms, *mp_surf_props, *mp_wall_state };
    mp_surf_solver = Mutation::Utilities::Config::Factory<SurfaceBalanceSolver>::create( m_gsi_mechanism, l_data_surface_balance_solver );

}

//======================================================================================

GasSurfaceInteraction::~GasSurfaceInteraction(){
    if ( mp_surf_props != NULL ) { delete mp_surf_props; }
    if ( mp_wall_state != NULL ) { delete mp_wall_state; }
    if ( mp_surf_solver != NULL ) { delete mp_surf_solver; }
}

//======================================================================================

void GasSurfaceInteraction::setWallState( const double* const l_mass, const double* const l_energy, const int state_variable ){
    mp_wall_state->setWallState( l_mass, l_energy, state_variable );
}

//======================================================================================

void GasSurfaceInteraction::getWallState( double* const l_mass, double* const l_energy, const int state_variable ){
    mp_wall_state->getWallState( l_mass, l_energy, state_variable );
}

//======================================================================================

void GasSurfaceInteraction::surfaceProductionRates( double* const lp_mass_prod_rate ){

    mp_surf_solver->computeGSIProductionRate( v_mass_prod_rate );

    for ( int i_sp = 0 ; i_sp < v_mass_prod_rate.size() ; ++i_sp ){
        lp_mass_prod_rate[ i_sp ] = v_mass_prod_rate( i_sp );
    }

}

//======================================================================================

void GasSurfaceInteraction::setDiffusionModel( const double* const l_mole_frac_edge, const double& l_dx ){

    for ( int i_ns = 0 ; i_ns < m_thermo.nSpecies() ; ++i_ns ){
        v_mole_frac_edge(i_ns) = *(l_mole_frac_edge + i_ns);
    }

    mp_surf_solver->setDiffusionModel( v_mole_frac_edge, l_dx );

}

//======================================================================================

void  GasSurfaceInteraction::solveSurfaceBalance(){

//    mp_surf_solver->solveSurfaceBalance( mp_surf_descr->getWallRhoi(), mp_surf_descr->getWallT() );
    mp_surf_solver->solveSurfaceBalance( mp_wall_state->getWallRhoi(), mp_wall_state->getWallT() );

}

//======================================================================================

inline void GasSurfaceInteraction::locateGSIInputFile( std::string& l_gsi_input_file ){

    // Check if the file is in the current directory
    l_gsi_input_file = l_gsi_input_file + ".xml";
    std::ifstream file( l_gsi_input_file.c_str() , std::ios::in );
    
    // If it is not, check in MPP_DATA_DIRECTORY/gsi
    if ( !file.is_open() ){
        l_gsi_input_file = Mutation::Utilities::getEnvironmentVariable("MPP_DATA_DIRECTORY") + "/gsi/" + l_gsi_input_file;
    }

    /** @todo FIX IT HERE */
    // Give and error

}

//======================================================================================

inline void GasSurfaceInteraction::errorWrongTypeofGSIFile( const std::string& l_gsi_root_tag ){

    if ( l_gsi_root_tag != "gsi" ) {
        std::cerr << "Root element in Gas Surface Interaction input file " << l_gsi_root_tag << " is not of 'gassurfaceinteraction' type!" << std::endl; // @todo FIX ERROR not l_gsi_root_tag. Instead name of file...
        exit(1);
    }

}

//======================================================================================

inline void GasSurfaceInteraction::errorInvalidGSIFileProperties( const std::string& l_gsi_option ) {

    std::cerr << l_gsi_option << " is not a valid gas surface interaction file option!" << std::endl;
    exit(1);

}

//======================================================================================

    } // namespace GasSurfaceInteraction 
} // namespace Mutation
