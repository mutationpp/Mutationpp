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

    // Initializing
    //std::vector<double> v_rhoi( ns, 0.0 );

    //std::vector<double> v_Vd( ns, 0.0 );
    //std::vector<double> v_Vd_ram( ns, 0.0 );
    //std::vector<double> v_Vd_sm( ns, 0.0 );
    //std::vector<double> v_b( ns, 0.0);    

    Eigen::MatrixXd m_ram_Dij(ns, ns);
    Eigen::VectorXd v_Vd_ram(ns);
    Eigen::VectorXd v_Vd_sm(ns);
    Eigen::VectorXd v_Vd(ns);
    Eigen::VectorXd v_b(ns);
    Eigen::VectorXd v_rhoi(ns);

    double E_field = 0.E0;

    double sum;
    double sum_ram;
    double sum_sm;

    double ** const p_Dij = new double* [ns] ;
    for ( int i = 0; i < ns; ++i ){
        p_Dij[i] = new double [ns];
    }

    // Loop over range of temperatures and compute equilibrium values at 1 atm    
    double P = ONEATM;

    // L = 1
    v_b[0] = 1.0 * ( ns - 1 );
    for ( int i = 1 ; i < ns ; i++ ){
        v_b[i] = - 1.0;
    }

    for (int k = 0; k < 6; ++k) {
        // Compute the temperature
        double T = 500.0 + static_cast<double>(k) * 5000.0;
        
        // Set the mixture state equal to the equilibrium state for the given
        // temperature and pressure
        mix.setState(&T, &P);
        double rho = mix.density();
        for ( int i = 0 ; i < ns ; ++i ){
            v_rhoi[i] = rho * mix.Y()[i];
        }

        // Exact Diffusion Coefficients
        mix.exactDiffusionMatrix( p_Dij );

        for ( int i = 0 ; i < ns ; ++i ) {
            v_Vd[i] = 0.0;
            v_Vd_ram[i] = 0.0;
        }
        for ( int i = 0; i < ns; ++i ){
            for ( int j = 0; j < ns; ++j ){
                v_Vd[i] -= p_Dij[i][j] * v_b[j];
            }
        }

        // Ramsaw Diffusion Coefficients
        m_ram_Dij = mix.diffusionMatrix();
        v_Vd_ram = -m_ram_Dij*v_b;
        //for ( int i = 0; i < ns; ++i ){
        //    for ( int j = 0; j < ns; ++j ){
        //        v_Vd_ram[i] -= m_ram_Dij(i, j) * v_b[j];
        //    }
       // }

        // Comparing with Stefan Maxwell
        mix.stefanMaxwell( &v_b[0], &v_Vd_sm[0], E_field );
        std::cout << "Stefan Maxwell E-field = " << E_field << std::endl;

        // Sums
        sum = 0.0;
        sum_ram = 0.0;
        sum_sm = 0.0;
        for( int i = 0 ; i < ns ; i++ ){
            sum += v_rhoi[i] * v_Vd[i];
            sum_ram += v_rhoi[i] * v_Vd_ram[i];
            sum_sm += v_rhoi[i] * v_Vd_sm[i];
        }

        // Output
        std::cout << "Temperature = " << T << std::endl;
        std::cout << "Fluxes Sum = " << sum << " Sum Ram = " << sum_ram << " Sum SM = " << sum_sm << std::endl;

        std::cout << "Diffusion Matrix" << std::endl;
        for ( int i = 0; i < ns; ++i ){
            for ( int j = 0; j < ns; ++j ){
                std::cout << std::setw(13) << p_Dij[i][j] ;
            }
            std::cout << std::endl;
        }
            std::cout << std::endl;
        std::cout << "Ramsaw Diffusion Matrix" << std::endl;
        for ( int i = 0; i < ns; ++i ){
            for ( int j = 0; j < ns; ++j ){
                std::cout << std::setw(13) << m_ram_Dij(i,j) ;
            }
            std::cout << std::endl;
        }

        std::cout << std::endl;
   
        std::cout << "Velocity Error against Ramsaw" << std::endl;;
        for ( int i = 0; i < ns; ++i ){
            std::cout << std::abs((v_Vd_ram[i] - v_Vd[i])/v_Vd_ram[i]) << " ";
        }
        std::cout << std::endl;
        std::cout << "Velocity Error against SM" << std::endl;
        for ( int i = 0; i < ns; ++i ){
            std::cout << std::abs((v_Vd_sm[i] - v_Vd[i])/v_Vd_sm[i]) << " ";
        }
        
        std::cout << "\nRamshaw Velocities: \n" << v_Vd_ram << std::endl;
        std::cout << "\nSM Velocities: \n" << v_Vd_sm << std::endl;
        std::cout << "\nExact Matrix Velocities: \n" << v_Vd << std::endl;
        
        std::cout << "\nRamshaw Fluxes: \n" << v_Vd_ram.array()*v_rhoi.array() << std::endl;
        std::cout << "\nSM Fluxes: \n" << v_Vd_sm.array()*v_rhoi.array() << std::endl;
        std::cout << "\nExact Matrix Fluxes: \n" << v_Vd.array()*v_rhoi.array() << std::endl;
	
	    std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << " " << std::endl;


    } // Loop over temperatures
 
    // Cleaning
    for ( int i = 0 ; i < ns; ++i ){
        delete [] p_Dij[i];
    }

    delete [] p_Dij;

    return 0;
}
