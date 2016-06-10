#include <cstddef>

#include "InputData.h"

InputData::InputData( Mutation::Mixture* const p_mixture, std::ifstream& l_input_file )
                    : m_ns( p_mixture->nSpecies() ),
                      m_nEnergyEqns( p_mixture->nEnergyEqns() ),
                      v_rhoi( m_ns, 0.E0 ),
                      v_X( m_ns, 0.E0 ),
                      v_temp( m_nEnergyEqns, 0.E0 ),
                      set_pressure( 0 ),
                      set_temperature( 0 ),
                      set_velocity( 0 ),
                      set_state_rhoi_T( 1 ) {

    while( getline( l_input_file, line ) ){
        // Parsing Comments
        if( line[0] == '-' ) continue;
        // Parsing Pressure
        if ( line.compare( "Pressure [Pa] (Pre-shock Conditions):" ) == 0 ) {
            getline( l_input_file, line );
            std::stringstream ss(line);
            ss >> m_P;
            set_pressure = 1;
        }
        // Parsing Temperature
        if (line.compare( "Temperature [K] (Pre-shock Conditions):" ) == 0 ) {
            getline( l_input_file, line );
            std::stringstream ss(line);
            ss >> v_temp[0];
            set_temperature = 1;
        }
        if (line.compare( "Shock speed [m/s]:" ) == 0 ) {
            getline( l_input_file, line );
            std::stringstream ss(line);
            ss >> m_Vs;
            set_velocity = 1;
        }
        if (line.compare( "End Pre Shock Conditions" ) == 0 ){ break; }
        // @todo GIVE ERROR!
    }

    // @todo ERROR PUT ME IN A FUNCTION
    if ( !set_velocity || !set_temperature || !set_velocity ){
        std::cout << "ERROR SOMETHING NOT SET" << std::endl;
        exit(1);
    }

    for ( int i_nEn = 1; i_nEn < m_nEnergyEqns; i_nEn++ ){
        v_temp[i_nEn] = v_temp[0];
    }

    // Equilibrium Composition
    p_mixture->equilibriumComposition( v_temp[0], m_P, &v_X[0] );
    m_rho = p_mixture->density( v_temp[0], m_P, &v_X[0] );
    p_mixture->convert<Mutation::Thermodynamics::X_TO_Y>( &v_X[0], &v_rhoi[0] );
    
    for ( int i_ns = 0; i_ns < m_ns; ++i_ns){
        v_rhoi[ i_ns ] *= m_rho;
    }

}
