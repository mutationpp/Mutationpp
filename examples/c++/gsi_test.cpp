/**
 * @file gsi_test.cpp
 */

// Must include this header file to use the Mutation++ library
#include "mutation++.h"
#include <iostream>
#include <fstream>
#include <cmath>

//using namespace std;

using std::ofstream;
using std::setw;
using std::cout;
using std::endl;

using namespace Mutation;
using namespace Mutation::Numerics;

int main()
{
//================================================================================================================
    // Initializing the library
    Mutation::MixtureOptions opts("air5gsi");
    opts.setStateModel("ChemNonEq1T");
    opts.setThermodynamicDatabase("RRHO");
    Mutation::Mixture mix(opts);

//================================================================================================================

    const int ns   = mix.nSpecies();
    double* p_rhoi = new double [ns];
    double* p_output = new double [ns];
    double* p_rhoeint = new double [1];

//================================================================================================================
    // Useful to know the species indices
    const int iN   = mix.speciesIndex("N");
    const int iO   = mix.speciesIndex("O");
    const int iNO  = mix.speciesIndex("NO");
    const int iN2  = mix.speciesIndex("N2");
    const int iO2  = mix.speciesIndex("O2");

//============================================================================
    // Conditions T = 10000K and p = 1000Pa
    // Initial partial densities

    p_rhoi[iN] = 0.000133077;
    p_rhoi[iO] = 4.04114e-05;
    p_rhoi[iNO] = 3.79216e-10;
    p_rhoi[iN2] = 1.35089e-08;
    p_rhoi[iO2] = 5.56891e-12;

    double rho = p_rhoi[iN] + p_rhoi[iO] + p_rhoi[iNO] + p_rhoi[iN2] + p_rhoi[iO2];
    
    *p_rhoeint = 10000.0e0;

//===========================================================================
    // Initializing Mutation++
    
    int state = 1;
    //mix.setState(p_rhoi,p_rhoeint, state); 
    mix.setWallState(p_rhoi, p_rhoeint);
    //mix.netGSIProductionRates(p_output);
    
    std::cout << "Output = " << std::endl;
    //std::cout << "Output = " << p_output[0] << endl;
 
//================================================================================================================ 
    // Clean up
    delete [] p_rhoi;
    delete [] p_rhoeint;
    delete [] p_output;
}
