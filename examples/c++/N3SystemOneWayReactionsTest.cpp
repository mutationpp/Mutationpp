#include "mutation++.h"

//========================== GLOBAL CONDITIONS =================================
const static double Temperature = 300;     // K
const static double Density     = 1.e-5;   // kg/m^3

const static double dT          = 10.e-7;  // s

static void initializeComposition( std::vector<double>& l_comp ){
    double l_tolerance = 1.e-20;
    double l_initial_seeding = 1.e-1; // 10 per mass%
    size_t pos_N          = 0;
	size_t pos_N2_ground  = 1;

	for ( std::vector<double>::iterator iter = l_comp.begin(); iter != l_comp.end() ; ++iter )
	    *iter = l_tolerance;

    l_comp[ pos_N ]         = l_initial_seeding * Density;
    l_comp[ pos_N2_ground ] = (1 - l_initial_seeding ) * Density;
}

//===================================================================================

class DataHandler{
public:
	DataHandler( Mutation::Mixture& l_mix )
               : m_ns ( l_mix.nSpecies() ),
				 m_T ( Temperature ),
                 m_rho( Density ),
				 mv_rhoi( m_ns ) {
		initializeComposition( mv_rhoi );
    }

	~DataHandler(){ }

	const double& getTemperature() const { return m_T; }
	const std::vector<double>& getPartialDensities() const { return mv_rhoi; }

private:
	size_t m_ns;

	double m_T;
	double m_rho;
	std::vector<double> mv_rhoi;
};

//===================================================================================

class ChemistryRateSolver{
public:
	ChemistryRateSolver()
                       : m_dT( dT ) {}
    ~ChemistryRateSolver(){}

    void solve(){


    	// Printing Results


    }

private:
	double m_dT;

	printResults();

};

//===================================================================================

int main(){

	// Initialize Mutation++
    Mutation::MixtureOptions opts("O2");
    opts.setStateModel("ChemNonEq1T");
    opts.setThermodynamicDatabase("RRHO");
    Mutation::Mixture mix(opts);

    // Initializing Conditions
    DataHandler m_data_handler( mix );
    ChemistryRateSolver m_chem_rate_solver();

    m_chem_rate_solver.solve();
	// solve dw_i^(n+1) = F^n * dt with F^n = law of mass action


	// mixture.getProductionRates()

    return 0;

//    // Setup arrays
//    double rhoi [mix.nSpecies()];
//    double wdot [mix.nSpecies()];
//
//    // Set state of mixture to the initial conditions
//    rhoi[mix.speciesIndex("O")]  = 0.0;
//    rhoi[mix.speciesIndex("O2")] = rho_init;
//    mix.setState(rhoi, &T_init, 1); // set state using {rho_i, T}
//
//    // Get the mixture energy which is conserved during the simulation
//    double etot = mix.mixtureEnergyMass() * mix.density(); // J/m^3
//
//    // Write the results header and initial conditions
//    printHeader(mix);
//    printResults(0.0, mix);
//
//    // Integrate the species mass conservation equations forward in time
//    double dt   = 0.0;
//    double time = 0.0;
//    double tout = 0.0001;
//
//    while (time < 100.0) {
//        // Get the species production rates
//        mix.netProductionRates(wdot);
//
//        // Compute "stable" time-step based on maximum allowed change in
//        // species densities
//        dt = 0.0;
//        for (int i = 0; i < mix.nSpecies(); ++i)
//            dt += 1.0e-7 * std::abs(rho_init * mix.Y()[i] / wdot[i]);
//        dt = std::min(dt / mix.nSpecies(), tout-time);
//
//        // Integrate in time
//        for (int i = 0; i < mix.nSpecies(); ++i)
//            rhoi[i] += wdot[i]*dt;
//        time += dt;
//
//        // Set the new state of the mixture (using conserved variables)
//        mix.setState(rhoi, &etot);
//
//        // Print results at specified time intervals
//        if (std::abs(tout-time) < 1.0e-10) {
//           printResults(time, mix);
//
//           if (tout > 9.9)          tout += 10.0;
//           else if (tout > 0.99)    tout += 1.0;
//           else if (tout > 0.099)   tout += 0.1;
//           else if (tout > 0.0099)  tout += 0.01;
//           else if (tout > 0.00099) tout += 0.001;
//           else                     tout += 0.0001;
//        }
//    }
//

}

//===================================================================================

// // Prints the header of the output table
// void printHeader(const Mutation::Mixture& mix)
// {
//     std::cout << std::setw(8) << "Time[s]";
//     std::cout << std::setw(12) << "T[K]";
//     std::cout << std::setw(12) << "P[Pa]";
//     std::cout << std::setw(12) << "rho[kg/m^3]";
//     std::cout << std::setw(12) << "e[J/kg]";
//     for (int i = 0; i < mix.nSpecies(); ++i)
//         std::cout << std::setw(12) << "Y_" + mix.speciesName(i);
//     std::cout << std::endl;
// }
//
// // Prints the mixture properties at the given time
// void printResults(double time, const Mutation::Mixture& mix)
// {
//     std::cout << std::setw(8) << time;
//     std::cout << std::setw(12) << mix.T();
//     std::cout << std::setw(12) << mix.P();
//     std::cout << std::setw(12) << mix.density();
//     std::cout << std::setw(12) << mix.mixtureEnergyMass();
//     for (int i = 0; i < mix.nSpecies(); ++i)
//         std::cout << std::setw(12) << mix.Y()[i];
//     std::cout << std::endl;
// }
