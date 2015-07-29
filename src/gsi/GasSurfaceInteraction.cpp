#include "GasSurfaceInteraction.h"
#include "Utilities.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//======================================================================================

GasSurfaceInteraction::GasSurfaceInteraction( Mutation::Thermodynamics::Thermodynamics& l_thermo, Mutation::Transport::Transport& l_transport, std::string l_gsi_input_file )
                                            : m_thermo( l_thermo ),
                                              m_transport( l_transport ),
                                              mp_surf_descr( NULL ),
                                              mp_surf_solver( NULL ),
                                              v_mass_prod_rate( m_thermo.nSpecies() ),
                                              v_mole_frac_edge( m_thermo.nSpecies() ) {

    if ( l_gsi_input_file == "none" ){ return; }

    locateGSIInputFile( l_gsi_input_file );
    
    Mutation::Utilities::IO::XmlDocument l_xml_doc( l_gsi_input_file );
    Mutation::Utilities::IO::XmlElement l_root_element = l_xml_doc.root();

    errorWrongTypeofGSIFile( l_root_element.tag() );

    l_root_element.getAttribute( "gsi_mechanism", m_gsi_mechanism, "none" );

    Mutation::Utilities::IO::XmlElement::const_iterator xml_position_surf_descr;
    Mutation::Utilities::IO::XmlElement::const_iterator xml_position_diff_model;
    Mutation::Utilities::IO::XmlElement::const_iterator xml_position_prod_terms;
    
    getXmlPositionPointerSurfPropsDiffModelProdTerm( xml_position_surf_descr, xml_position_diff_model, xml_position_prod_terms, l_root_element );

    mp_surf_descr = new SurfaceDescription( m_thermo, m_gsi_mechanism, *xml_position_surf_descr );
    mp_surf_solver = new SurfaceBalanceSolver( m_thermo, m_transport, m_gsi_mechanism, *xml_position_diff_model, *xml_position_prod_terms, *mp_surf_descr );

}

//======================================================================================

GasSurfaceInteraction::~GasSurfaceInteraction(){

    if ( mp_surf_descr != NULL ) { delete mp_surf_descr; }
    if ( mp_surf_solver != NULL ) { delete mp_surf_solver; }

}

//======================================================================================

void GasSurfaceInteraction::setWallState( const double* const l_mass, const double* const l_energy, const int state_variable ){

    mp_surf_solver->setWallState( l_mass, l_energy, state_variable );

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

//    mp_surf_solver->

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

}

//======================================================================================

inline void GasSurfaceInteraction::errorWrongTypeofGSIFile( const std::string& l_gsi_root_tag ){

    if ( l_gsi_root_tag != "gsi" ) {
        std::cerr << "Root element in Gas Surface Interaction input file " << l_gsi_root_tag << " is not of 'gassurfaceinteraction' type!" << std::endl;
        exit(1);
    }

}

//======================================================================================

inline void GasSurfaceInteraction::getXmlPositionPointerSurfPropsDiffModelProdTerm( Mutation::Utilities::IO::XmlElement::const_iterator& l_index_surf_descr, Mutation::Utilities::IO::XmlElement::const_iterator& l_index_diff_model, Mutation::Utilities::IO::XmlElement::const_iterator& l_index_prod_terms, const Mutation::Utilities::IO::XmlElement& root ){

// Replace it with getTag();

        Mutation::Utilities::IO::XmlElement::const_iterator iter = root.begin();
        for ( ; iter != root.end(); ++iter) {
            if ( iter->tag() == "surface_properties" ) {
                l_index_surf_descr = iter;
            } else if ( iter->tag() == "diffusion_model" ) {
                l_index_diff_model = iter;
            } else if ( iter->tag() == "production_terms" ) {
                l_index_prod_terms = iter;
            } else {
                inline void errorInvalidGSIFileProperties( );
            }
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
