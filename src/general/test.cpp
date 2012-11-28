#include <cstdlib>

#include "mutation++.h"
#include <fenv.h>

using namespace std;


/*template <typename Thermo, typename Equil, typename Trans, typename Kin>
void makeCEQInput(Mixture<Thermo, Equil, Trans, Kin> *p_mixture, double *p_c)
{
    const int ns = p_mixture->nSpecies();
    const int ne = p_mixture->nElements();
    
    // Write the mixture size
    cout << ne << endl;
    cout << ns << endl;
    
    // Write out elemental compositions
    for (int i = 0; i < ne; ++i)
        cout << p_c[i] << endl;
    
    // Write out the species names
    for (int i = 0; i < ns; ++i)
        cout << p_mixture->speciesName(i) << endl;

    // Write out the element matrix
    for (int i = 0; i < ne; ++i) {
        for (int j = 0; j < ns; ++j)
            cout << setw(4) << p_mixture->elementMatrix()(j,i);
        cout << endl;
    }
    
    cout.precision(7);
    cout.setf(ios::scientific | ios::right);
    
    // Write out the species molecular weights
    for (int i = 0; i < ne; ++i)
        cout << setw(15) << p_mixture->atomicMass(i);
    cout << endl;
    
    // Next write out the thermodynamic polynomials in the order that CEQ wants
    double *p_temp = new double [ns];
    
    std::vector< std::pair<int, Nasa7Polynomial> >::const_iterator iter = 
        p_mixture->mp_thermodb->m_nasa_7_polynomials.begin();
    std::vector< std::pair<int, Nasa7Polynomial> >::const_iterator end  =
        p_mixture->mp_thermodb->m_nasa_7_polynomials.end();
    
    // upper temperature on lower temp-range
    for ( ; iter != end; ++iter)
        p_temp[iter->first] = iter->second.m_mid_temp;
    
    for (int i = 0; i < ns; ++i)
        cout << setw(15) << p_temp[i];
    cout << endl;
    
    // lower temp-range polynomial
    for (int k = 0; k < 7; ++k) {
        iter = p_mixture->mp_thermodb->m_nasa_7_polynomials.begin();    
        for ( ; iter != end; ++iter)
            p_temp[iter->first] = iter->second.mp_coefficients[0][k];
        
        for (int i = 0; i < ns; ++i)
            cout << setw(15) << p_temp[i];
        
        cout << endl;
    }
    
    // upper temp-range polynomial
    for (int k = 0; k < 7; ++k) {
        iter = p_mixture->mp_thermodb->m_nasa_7_polynomials.begin();    
        for ( ; iter != end; ++iter)
            p_temp[iter->first] = iter->second.mp_coefficients[1][k];
        
        for (int i = 0; i < ns; ++i)
            cout << setw(15) << p_temp[i];
        
        cout << endl;
    }
    exit(1); 
}*/


int main(int argc, char** argv)
{
    feenableexcept(FE_OVERFLOW|FE_DIVBYZERO|FE_INVALID);
    //feenableexcept(FE_OVERFLOW);
    //feenableexcept(FE_DIVBYZERO);

    if (argc <= 1) exit(1);
    
    string name(argv[argc-1]);
    cout << name << endl;
    
    Mixture mix(name);
        
    const int ne = mix.nElements();
    const int ns = mix.nSpecies();
    const int nr = mix.nReactions();
    
    double* p_c = new double [ne];
    double* p_x = new double [ns];
    double* p_y = new double [ns];
    
    for (int i = 0; i < ne; ++i)
        p_c[i] = mix.getDefaultComposition(i);
    
    // Write out a header
    const int width = 13;
    cout << setw(width) << "T(K)";
    cout << setw(width) << "rho(kg/m^3)";    
    cout << setw(width) << "Mw(kg/mol)";
    cout << setw(width) << "Cp(J/kg-K)";
    cout << setw(width) << "Cv(J/kg-K)";
    cout << setw(width) << "gamma";
    cout << setw(width) << "H(J/kg)";
    cout << setw(width) << "S(J/kg-K)";
    //cout << setw(width) << "mu(Pa-s)";
    //cout << setw(width) << "Lambda";
    
    for (int i = 0; i < ns; ++i)
        cout << setw(width) << "X_" + mix.speciesName(i);
    cout << endl;
    
    for (int i = 0; i < 199; ++i) {
        // Temperature and pressure
        double T  = static_cast<double>(i) * 100.0 + 200.0;
        double P  = 101325.0 * 0.1;
        
        // Get equilibrium mole and mass fractions
        mix.equilibrate(T, P, p_c, p_x);
        mix.convert<X_TO_Y>(p_x, p_y);
        
        // Write out some mixture quantities
        cout << setw(width) << T;
        cout << setw(width) << mix.density();
        cout << setw(width) << mix.mixtureMw();
        cout << setw(width) << mix.mixtureEquilibriumCpMass();
        cout << setw(width) << mix.mixtureEquilibriumCvMass();
        cout << setw(width) << mix.mixtureEquilibriumGamma();
        cout << setw(width) << mix.mixtureHMass();
        cout << setw(width) << mix.mixtureSMass();
        //cout << setw(width) << mix.eta();
        //cout << setw(width) << mix.lambda();
        
        for (int i = 0; i < ns; ++i)
            cout << setw(width) << p_x[i];
        cout << endl;
    }
    
    delete [] p_c;
    delete [] p_x;
    delete [] p_y;

    return 0;
}
