#include <cstddef>
#include <iostream>
#include <stdlib.h>

#include "InputFileParser.h"

InputFileParser::InputFileParser( const std::string& l_input_file_name ) 
                                : m_input_file_name( l_input_file_name ),
                                  mp_mixture( NULL ), 
                                  mp_mixture_options( NULL ),
                                  mp_mesh_options( NULL ) {
    
    locateInputFile();

    // Parsing input file for Mixture Options
    while( getline( m_input_file, line ) ) {
        if( line[0] == '-' ) continue;
        if ( line.compare( "Name of the mixture:" ) == 0 ) {
            getline( m_input_file, m_mixture );
        }
        if ( line.compare("State Model:") == 0 ) {
            getline( m_input_file, m_state_model );
        }
        if ( line.compare("Thermodynamic Database:") == 0 ) {
            getline( m_input_file, m_thermo_db);
        }
        if( line.compare("End Mixture Options") == 0 ){ break; }
    }

    // Initializing Mutation++ Library
    mp_mixture_options = new Mutation::MixtureOptions( m_mixture.c_str() );
    mp_mixture_options->setStateModel( m_state_model.c_str() );
    mp_mixture_options->setThermodynamicDatabase( m_thermo_db.c_str() );
    mp_mixture = new Mutation::Mixture( *mp_mixture_options );

    // Input Data
    mp_input_data = new InputData( mp_mixture, m_input_file );

    // Mesh Options
    while( getline( m_input_file, line ) ) {
        if (line.compare("Mesh integration options:") == 0 ) { 
//            mp_mesh_options = new MeshOptions( /* m_input_file */ );
        }
    }
    
    // Closing File
    m_input_file.close();
       
}

//=============================================================================================================

InputFileParser::~InputFileParser(){
    if ( mp_mixture != NULL ) delete mp_mixture;
    if ( mp_mixture_options != NULL ) delete mp_mixture_options;
    if ( mp_mesh_options != NULL ) delete mp_mesh_options;
}

//=============================================================================================================

inline void InputFileParser::locateInputFile(){

    // Check if the file is in the current directory
    m_input_file_name = m_input_file_name + ".in";
    m_input_file.open( m_input_file_name.c_str() , std::ios::in );
    
    // If it is not, check in ../input/
    if ( !m_input_file.is_open() ){
        m_input_file_name = "../input/" + m_input_file_name;
    }
    m_input_file.open( m_input_file_name.c_str() , std::ios::in );

    if ( !m_input_file.is_open() ){
        std::cerr << "Cannot Locate Input File: "<< m_input_file_name << std::endl;
        exit(1);
    }

}

//=============================================================================================================

