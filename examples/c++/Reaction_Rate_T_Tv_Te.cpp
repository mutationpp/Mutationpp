/**
 * @file Reaction_Rate_T_Tv_Te.cpp
 *
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
    Mutation::MixtureOptions opts("test_mixture");
    opts.setStateModel("TTvTeRhoi");
    opts.setThermodynamicDatabase("RRHO");
    Mutation::Mixture mix(opts);

//================================================================================================================
    // Useful to know the species indices
    const int ns   = mix.nSpecies();
    const int ie   = mix.speciesIndex("e-");
    const int iN   = mix.speciesIndex("N");
    const int iNp  = mix.speciesIndex("N+");
    const int iO   = mix.speciesIndex("O");
    const int iOp  = mix.speciesIndex("O+");
    const int iNO  = mix.speciesIndex("NO");
    const int iN2  = mix.speciesIndex("N2");
    const int iN2p = mix.speciesIndex("N2+");
    const int iO2  = mix.speciesIndex("O2");
    const int iO2p = mix.speciesIndex("O2+");
    const int iNOp = mix.speciesIndex("NO+");
    const int iOm  = mix.speciesIndex("O-");

//================================================================================================================
    // Initial mass fractions. They are not needed though
    double* pY = new double [ns];
    pY[ie] = 0.0;
    pY[iN] = 0.0;
    pY[iNp] = 0.0;
    pY[iO] = 0.0;
    pY[iOp] = 0.0;
    pY[iNO] = 0.0;
    pY[iN2] = 0.7;
    pY[iN2p] = 0.0;
    pY[iO2] = 0.3;
    pY[iO2p] = 0.0;
    pY[iNOp] = 0.0;
    pY[iOm] = 0.0;

//================================================================================================================
    // Initial densities
    double rho = 2.85e-4;  // Give the density of the mixture in kg/m^3 
    double rho_is [ns];
    for (int ic_1 = 0; ic_1 < ns; ic_1++){
      rho_is[ic_1] = pY[ic_1] * rho;
    }

//================================================================================================================
    // Loop for the temperature matrix

    int n_temp_size = 100;

    double TTvTe_min [3];
    TTvTe_min[0] = 100.0;   // Minimum Translational temperature in K
    TTvTe_min[1] = 1310.0;   // Minimum Vibrational temperature in K
    TTvTe_min[2] = 2520.0;   // Minimum Electronic temperature in K

    double TTvTe_max [3];
    TTvTe_max[0] = 50000.0;   // Maximum Translational temperature in K
    TTvTe_max[1] = 51310.0;   // Maximum Vibrational temperature in K
    TTvTe_max[2] = 52520.0;   // Maximum Electronic temperature in K

    double TTvTe [3]; 
    double* reaction_values = new double [mix.nReactions()];

    ofstream f_reaction_rates;
    f_reaction_rates.open ("Reaction_Rates.dat", std::ios::out);

    f_reaction_rates << "Forward Reaction Rates" << endl;

    f_reaction_rates << setw(16) << "T_trans [K]";
    f_reaction_rates << setw(16) << "T_vib [K]";
    f_reaction_rates << setw(16) << "T_ele [K]";
    f_reaction_rates << setw(16) << "T_park [K]";

    for (int ic_4 = 0; ic_4 < mix.nReactions(); ic_4++)
    {
    f_reaction_rates << setw(12) << "Reaction" << setw(4) << ic_4;
    }
    f_reaction_rates << endl;


    for (int ic_3 = 1; ic_3 < n_temp_size; ic_3++)
    {
      TTvTe[0] = (TTvTe_max[0]-TTvTe_min[0])/(n_temp_size - 1) * (ic_3 - 1) + TTvTe_min[0];
      TTvTe[1] = (TTvTe_max[1]-TTvTe_min[1])/(n_temp_size - 1) * (ic_3 - 1) + TTvTe_min[1];
      TTvTe[2] = (TTvTe_max[2]-TTvTe_min[2])/(n_temp_size - 1) * (ic_3 - 1) + TTvTe_min[2];
    
      f_reaction_rates << std::fixed << setw(16) << TTvTe[0];
      f_reaction_rates << setw(16) << TTvTe[1];
      f_reaction_rates << setw(16) << TTvTe[2];
      f_reaction_rates << setw(16) << sqrt(TTvTe[0]*TTvTe[1]);

      mix.setState(TTvTe, rho_is);

      mix.forwardRateCoefficients(reaction_values);

      for (int ic_2 = 0; ic_2 < mix.nReactions(); ic_2++){
      f_reaction_rates << std::scientific << setw(16) << reaction_values[ic_2];
      }
      f_reaction_rates << endl;
   }

//================================================================================================================ 
//Backward reaction rates  

    f_reaction_rates << "Backward Reaction Rates" << endl;

    f_reaction_rates << setw(16) << "T_trans [K]";
    f_reaction_rates << setw(16) << "T_vib [K]";
    f_reaction_rates << setw(16) << "T_ele [K]";
    f_reaction_rates << setw(16) << "T_park [K]";

    for (int ic_5 = 0; ic_5 < mix.nReactions(); ic_5++)
    {
    f_reaction_rates << setw(12) << "Reaction" << setw(4) << ic_5;
    }
    f_reaction_rates << endl;

    for (int ic_6 = 1; ic_6 < n_temp_size; ic_6++)
    {
      TTvTe[0] = (TTvTe_max[0]-TTvTe_min[0])/(n_temp_size - 1) * (ic_6 - 1) + TTvTe_min[0];
      TTvTe[1] = (TTvTe_max[1]-TTvTe_min[1])/(n_temp_size - 1) * (ic_6 - 1) + TTvTe_min[1];
      TTvTe[2] = (TTvTe_max[2]-TTvTe_min[2])/(n_temp_size - 1) * (ic_6 - 1) + TTvTe_min[2];
    
      f_reaction_rates << std::fixed << setw(16) << TTvTe[0];
      f_reaction_rates << setw(16) << TTvTe[1];
      f_reaction_rates << setw(16) << TTvTe[2];
      f_reaction_rates << setw(16) << sqrt(TTvTe[0]*TTvTe[1]);

      mix.setState(TTvTe, rho_is);

      mix.backwardRateCoefficients(reaction_values);

      for (int ic_7 = 0; ic_7 < mix.nReactions(); ic_7++){
      f_reaction_rates << std::scientific << setw(16) << reaction_values[ic_7];
      }
      f_reaction_rates << endl;
   }


    f_reaction_rates.close();

    // !!! Comments:
    // How do I get the number of temperatures rigorously???

//================================================================================================================ 
    // Clean up
    delete [] pY;
    delete [] reaction_values;
}
