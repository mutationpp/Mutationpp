#include <sstream>

#include "MeshOptions.h"

MeshOptions::MeshOptions( std::ifstream& l_input_file ){

    // Parsing for Mesh Options
    while( getline( l_input_file, line ) ){
        if( line[0] == '-' ) continue;
        if ( line.compare( "x_0" ) == 0 ) {
            getline( l_input_file, line );
            std::stringstream ss(line);
            ss >> m_x_init;
        }
        if ( line.compare( "x_end" ) == 0 ) {
            getline( l_input_file, line );
            std::stringstream ss(line);
            ss >> m_x_end;
        }
        if ( line.compare( "dx" ) == 0 ) {
            getline( l_input_file, line );
            std::stringstream ss(line);
            ss >> m_dx;
        }
        if (line.compare( "End Mesh Options" ) == 0 ){ break; }
        /** @todo Parse Error */

    }

    createMesh();

}

void MeshOptions::createMesh(){

    double x = m_x_init;
    double x_tol = 1.E-16;

    while ( x < m_x_end + x_tol ){
        m_mesh.push_back( x );
        x += m_dx;
    }

    n_mesh_points = m_mesh.size();

}
