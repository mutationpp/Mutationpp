
#include "mutation++.h"
#include "MultiPhaseEquilSolver.h"

#ifdef _GNU_SOURCE
#include <fenv.h>
#endif

using namespace Mutation::Thermodynamics;
using namespace std;
using Mutation::ONEATM;

void bishnuTable3()
{
    Mutation::MixtureOptions opts("water8");
    opts.setThermodynamicDatabase("NASA-9");
    Mutation::Mixture mix(opts);
    
    double c[mix.nElements()];
    double x[mix.nSpecies()];
    
    for (int i = 0; i < mix.nElements(); i++)
        c[i] = mix.getDefaultComposition(i);
    
    cout << endl;
    cout << "Bisnu et. al. Table 3:" << endl;
    cout << setw(5) << "Point" << setw(10) << "T (K)" << setw(10) << "P (atm)";
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(14) << "X_" + mix.speciesName(i);
    cout << endl;
    
    
    mix.equilibrate(373.15, ONEATM, c, x);
    cout << setw(5) << "BP" << setw(10) << 373.15 << setw(10) << 1;
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(14) << x[i];
    cout << endl;
    
    mix.equilibrate(473.15, ONEATM, c, x);
    cout << setw(5) << "G1" << setw(10) << 473.15 << setw(10) << 1;
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(14) << x[i];
    cout << endl;
    
    mix.equilibrate(541.15, ONEATM, c, x);
    cout << setw(5) << "G2" << setw(10) << 541.15 << setw(10) << 1;
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(14) << x[i];
    cout << endl;
    
    mix.equilibrate(1500.0, ONEATM, c, x);
    cout << setw(5) << "G3" << setw(10) << 1500.0 << setw(10) << 1;
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(14) << x[i];
    cout << endl;
    
    mix.equilibrate(647.25, 218.0*ONEATM, c, x);
    cout << setw(5) << "CP" << setw(10) << 647.25 << setw(10) << 218;
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(14) << x[i];
    cout << endl << endl;
}

void bishnuTable5()
{
    Mutation::Mixture mix("water8");
    
    const int ne = mix.nElements();
    
    double c[mix.nElements()+1];
    double x[mix.nSpecies()];
    
    for (int i = 0; i < mix.nSpecies(); ++i)
        x[i] = static_cast<double>(1);
    mix.addEquilibriumConstraint(x);
    
    c[mix.elementIndex("H")] = double(4);
    c[mix.elementIndex("O")] = double(2);
    
    cout << endl;
    cout << "Bisnu et. al. Table 5:" << endl;
    cout << setw(10) << "M" << setw(10) << "T (K)" << setw(10) << "P (atm)";
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(14) << "X_" + mix.speciesName(i);
    cout << endl;
    
    for (double m = 2.0; m < 6.0; m += 0.01) {
        c[ne] = m;
        mix.equilibrate(1500.0, ONEATM, c, x);
        cout << setw(10) << m << setw(10) << 1500 << setw(10) << 1;
        for (int i = 0; i < mix.nSpecies(); ++i)
            cout << setw(14) << x[i];
        cout << endl;
    }
    
    
    /*c[ne] = double(2);
    equil.equilibrate(1500.0, ONEATM, c, x);
    cout << setw(5) << "2" << setw(10) << 1500 << setw(10) << 1;
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(14) << x[i];
    cout << endl;
    
    c[ne] = 4.0;
    equil.equilibrate(1500.0, ONEATM, c, x);
    cout << setw(5) << "4" << setw(10) << 1500 << setw(10) << 1;
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(14) << x[i];
    cout << endl;
    
    c[ne] = 5.9999999999;
    equil.equilibrate(1500.0, ONEATM, c, x);
    cout << setw(5) << "6" << setw(10) << 1500 << setw(10) << 1;
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(14) << x[i];
    cout << endl << endl;*/
    
}

void dXdT()
{
    Mutation::Mixture mix("water8");
    
    const int ne = mix.nElements();
    const int ns = mix.nSpecies();
    
    double c[ne];
    double x[ns];
    
    c[mix.elementIndex("H")] = double(4);
    c[mix.elementIndex("O")] = double(2);
    
    cout << setw(10) << "T(K)";
    for (int i = 0; i < ns; ++i)
        cout << setw(14) << "X_" + mix.speciesName(i);
    for (int i = 0; i < ns; ++i)
        cout << setw(14) << "dX/dT_" + mix.speciesName(i);
    cout << setw(14) << "sum(dX/dT)";
    for (int i = 0; i < ns; ++i)
        cout << setw(14) << "fddX/dT_" + mix.speciesName(i);
    cout << setw(14) << "fdsum(dX/dT)" << endl;
    
    for (int i = 0; i < 115; ++i) {
        double T = 300.0 + 50.0*double(i);
        mix.equilibrate(T, ONEATM, c, x);
        
        cout << setw(10) << T;
        for (int k = 0; k < ns; ++k)
            cout << setw(14) << x[k];
        
        mix.dXidT(x);
        double sum = 0.0;
        for (int k = 0; k < ns; ++k) {
            sum += x[k];
            cout << setw(14) << x[k];
        }
        
        cout << setw(14) << sum;
        
        double Teps = T*1.0E-8;
        mix.equilibrate(T+Teps, ONEATM, c, x, false);
        sum = 0.0;
        for (int k = 0; k < ns; ++k) {
            double dxdt = (x[k]-mix.X()[k])/Teps;
            sum += dxdt;
            cout << setw(14) << dxdt;
        }
        
        cout << setw(14) << sum << endl;
    }
    
    
}

int main()
{
#ifdef _GNU_SOURCE
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

    bishnuTable3();
    bishnuTable5();
    
    cout << endl;
    dXdT();
}
