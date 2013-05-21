
#include "mutation++.h"
#include "MultiPhaseEquilSolver.h"

//#ifdef __GNU__
#include <fenv.h>
//#endif

using namespace Mutation::Thermodynamics;
using namespace std;

void bishnuTable3()
{
    Mutation::MixtureOptions opts("water8");
    opts.setThermodynamicDatabase("NASA-9");
    Mutation::Mixture mix(opts);
    MultiPhaseEquilSolver equil(mix);
    
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
    
    
    equil.equilibrate(373.15, ONEATM, c, x);
    cout << setw(5) << "BP" << setw(10) << 373.15 << setw(10) << 1;
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(14) << x[i];
    cout << endl;
    
    equil.equilibrate(473.15, ONEATM, c, x);
    cout << setw(5) << "G1" << setw(10) << 473.15 << setw(10) << 1;
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(14) << x[i];
    cout << endl;
    
    equil.equilibrate(541.15, ONEATM, c, x);
    cout << setw(5) << "G2" << setw(10) << 541.15 << setw(10) << 1;
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(14) << x[i];
    cout << endl;
    
    equil.equilibrate(1500.0, ONEATM, c, x);
    cout << setw(5) << "G3" << setw(10) << 1500.0 << setw(10) << 1;
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(14) << x[i];
    cout << endl;
    
    equil.equilibrate(647.25, 218.0*ONEATM, c, x);
    cout << setw(5) << "CP" << setw(10) << 647.25 << setw(10) << 218;
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(14) << x[i];
    cout << endl << endl;
}

void bishnuTable5()
{
    Mutation::Mixture mix("water8");
    MultiPhaseEquilSolver equil(mix);
    
    const int ne = mix.nElements();
    
    double c[mix.nElements()+1];
    double x[mix.nSpecies()];
    
    for (int i = 0; i < mix.nSpecies(); ++i)
        x[i] = 1.0;
    equil.addConstraint(x);
    
    c[mix.elementIndex("H")] = 4.0;
    c[mix.elementIndex("O")] = 2.0;
    
    cout << endl;
    cout << "Bisnu et. al. Table 5:" << endl;
    cout << setw(5) << "M" << setw(10) << "T (K)" << setw(10) << "P (atm)";
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(12) << "X_" + mix.speciesName(i);
    cout << endl;
    
    
    c[ne] = 2.0;
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
    
    c[ne] = 6.0;
    equil.equilibrate(1500.0, ONEATM, c, x);
    cout << setw(5) << "6" << setw(10) << 1500 << setw(10) << 1;
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(14) << x[i];
    cout << endl << endl;
    
}

int main()
{
//#ifdef __GNU__
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
//#endif

    bishnuTable3();
    bishnuTable5();
}
