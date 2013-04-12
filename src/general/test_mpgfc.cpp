
#include "mutation++.h"
#include "MultiPhaseEquilSolver.h"

using namespace Mutation::Thermodynamics;
using namespace std;

void bishnuTable3()
{
    Mutation::Mixture mix("water8");
    MultiPhaseEquilSolver equil(mix);
    
    double c[mix.nElements()];
    double x[mix.nSpecies()];
    
    for (int i = 0; i < mix.nElements(); i++)
        c[i] = mix.getDefaultComposition(i);
    
    cout << endl;
    cout << "Bisnu et. al. Table 3:" << endl;
    cout << setw(5) << "Point" << setw(10) << "T (K)" << setw(10) << "P (atm)";
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(12) << "X_" + mix.speciesName(i);
    cout << endl;
    
    
    equil.equilibrate(373.15, ONEATM, c, x);
    cout << setw(5) << "BP" << setw(10) << 373.15 << setw(10) << 1;
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(12) << x[i];
    cout << endl;
    
    equil.equilibrate(473.15, ONEATM, c, x);
    cout << setw(5) << "G1" << setw(10) << 473.15 << setw(10) << 1;
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(12) << x[i];
    cout << endl;
    
    equil.equilibrate(541.15, ONEATM, c, x);
    cout << setw(5) << "G2" << setw(10) << 541.15 << setw(10) << 1;
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(12) << x[i];
    cout << endl;
    
    equil.equilibrate(1500.0, ONEATM, c, x);
    cout << setw(5) << "G3" << setw(10) << 1500.0 << setw(10) << 1;
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(12) << x[i];
    cout << endl;
    
    equil.equilibrate(647.25, 218.0*ONEATM, c, x);
    cout << setw(5) << "CP" << setw(10) << 647.25 << setw(10) << 218;
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(12) << x[i];
    cout << endl << endl;
}

void bishnuTable5()
{
    Mutation::Mixture mix("water8");
    MultiPhaseEquilSolver equil(mix);
    
    double c[mix.nElements()];
    double x[mix.nSpecies()];
    
    c[mix.elementIndex("H")] = 4.0;
    c[mix.elementIndex("O")] = 2.0;
    
    cout << endl;
    cout << "Bisnu et. al. Table 5:" << endl;
    cout << setw(5) << "M" << setw(10) << "T (K)" << setw(10) << "P (atm)";
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(12) << "X_" + mix.speciesName(i);
    cout << endl;
    
    
    for (int i = 0; i < mix.nSpecies(); ++i)
        x[i] = 1.0;
    equil.addConstraint(x, 2.0);
    equil.equilibrate(1500.0, ONEATM, c, x);
    cout << setw(5) << "2" << setw(10) << 1500 << setw(10) << 1;
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(12) << x[i];
    cout << endl;
    
    for (int i = 0; i < mix.nSpecies(); ++i)
        x[i] = 1.0;
    equil.clearConstraints();
    equil.addConstraint(x, 4.0);
    equil.equilibrate(1500.0, ONEATM, c, x);
    cout << setw(5) << "4" << setw(10) << 1500 << setw(10) << 1;
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(12) << x[i];
    cout << endl;
    
    for (int i = 0; i < mix.nSpecies(); ++i)
        x[i] = 1.0;
    equil.clearConstraints();
    equil.addConstraint(x, 6.0);
    equil.equilibrate(1500.0, ONEATM, c, x);
    cout << setw(5) << "6" << setw(10) << 1500 << setw(10) << 1;
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(12) << x[i];
    cout << endl << endl;
    
}

int main()
{
    bishnuTable3();
    bishnuTable5();
}
