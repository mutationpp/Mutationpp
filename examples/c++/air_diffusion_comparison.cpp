#include "mutation++.h"
#include <vector>
#include <iostream>

#include <Eigen/Dense>

using namespace Mutation;
using namespace Mutation::Numerics;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Transport;

int main()
{
    // Generate the default options for the air11 mixture
    MixtureOptions opts("air5");
    opts.setStateModel("EquilTP");
    
    // Change from default thermodynamic database (RRHO) to NASA 9-coefficient
    // polynomial database
    opts.setThermodynamicDatabase("RRHO");
    
    // Load the mixture with the new options
    Mixture mix(opts);
    int ns = mix.nSpecies();

    Eigen::MatrixXd m_Dij(ns, ns);
    Eigen::MatrixXd m_ram_Dij(ns, ns);
    Eigen::VectorXd v_Vd_ram(ns);
    Eigen::VectorXd v_Vd_sm(ns);
    Eigen::VectorXd v_Vd(ns);
    Eigen::VectorXd v_b(ns);
    Eigen::VectorXd v_rhoi(ns);

    double E_field = 0.E0;

    double ** const p_Dij = new double* [ns] ;
    for ( int i = 0; i < ns; ++i ){
        p_Dij[i] = new double [ns];
    }

    // Loop over range of temperatures and compute equilibrium values at 1 atm    
    double P = ONEATM;
    double T;
    double rho;

    // Defining Driving Forces Respecting the constraint that they should add equal to 1.
    v_b.setConstant(-1.0);
    v_b(0) = 1.0 * ( ns - 1 );

    for (int k = 0; k < 250; ++k) {
        // Compute the temperature
        T = 200.0 + static_cast<double>(k) * 200.0;
        
        // Set the mixture state equal to the equilibrium state for the given
        // temperature and pressure
        mix.setState(&T, &P);
        rho = mix.density();
        for ( int i = 0 ; i < ns ; ++i ){
            v_rhoi(i) = rho * mix.Y()[i];
        }

        // Exact Diffusion Coefficients
        mix.exactDiffusionMatrix( p_Dij );

        v_Vd.setZero();
        
        for ( int i = 0; i < ns; ++i ){
            for ( int j = 0; j < ns; ++j ){
                v_Vd(i) -= p_Dij[i][j] * v_b(j);
            }
        }

        // Ramsaw Diffusion Coefficients
        m_ram_Dij = mix.diffusionMatrix();
        v_Vd_ram = -m_ram_Dij * v_b;

        // Comparing with Stefan Maxwell
        mix.stefanMaxwell( &v_b(0), &v_Vd_sm(0), E_field );

        // Output
        std::cout << "Sum of the Fluxes for T = " << T << "K" << std::endl;
        std::cout << "Exact Diffusion Matrix = " << v_rhoi.dot( v_Vd ) << std::endl;
        std::cout << "Ramsaw Diffusion Matrix = " << v_rhoi.dot( v_Vd_ram ) << std::endl;
        std::cout << "Stefan-Maxwell = " << v_rhoi.dot( v_Vd_sm ) << std::endl;
        std::cout << "Velocities L2 Norm error Between Exact matrix and Ramsaw = " << ( v_Vd_ram - v_Vd ).norm() << std::endl;
        std::cout << "Velocities L2 Norm error Between Exact matrix and Stefan-Maxwell = " << ( v_Vd_sm -v_Vd ).norm() << std::endl;
        std::cout << std::endl;

    } // Loop over temperatures
 
    // Cleaning
    for ( int i = 0 ; i < ns; ++i ){
        delete [] p_Dij[i];
    }
    delete [] p_Dij;

    return 0;
}
