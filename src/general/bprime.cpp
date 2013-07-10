
#include "mutation++.h"
#include <fenv.h>
#include <iostream>

using namespace std;
using namespace Mutation;



int main(int argc, char* argv[])
{
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);

    Mixture mix(argv[6]);
    
    const int ne = mix.nElements();
    const int ns = mix.nSpecies();
    
    double* p_Yke = new double [ne];
    double* p_Ykg = new double [ne];
    double* p_Xw  = new double [ns];
    
    // Run conditions
    double T1 = atof(argv[1]);
    double T2 = atof(argv[2]);
    double dt = atof(argv[3]);
    double P  = atof(argv[4])*ONEATM;
    double Bg = atof(argv[5]);
    double Bc, hw;
    
    p_Yke[mix.elementIndex("C")] = 0.0;
    p_Yke[mix.elementIndex("H")] = 0.0;
    p_Yke[mix.elementIndex("O")] = 0.21;
    p_Yke[mix.elementIndex("N")] = 0.79;
    mix.convert<Thermodynamics::XE_TO_YE>(p_Yke, p_Yke);
    
    p_Ykg[mix.elementIndex("C")] = 0.206;
    p_Ykg[mix.elementIndex("H")] = 0.679;
    p_Ykg[mix.elementIndex("O")] = 0.115;
    p_Ykg[mix.elementIndex("N")] = 0.0;
    mix.convert<Thermodynamics::XE_TO_YE>(p_Ykg, p_Ykg);
    
    cout << endl << setw(10) << "Tw[K]" << setw(15) << "B'c" << setw(15) << "hw[MJ/kg]";
    for (int i = 0; i < ns; ++i)
        cout << setw(15) << mix.speciesName(i);
    cout << endl;
    
    for (double T = T1; T < T2 + 1.0e-6; T += dt) {
        mix.surfaceMassBalance(p_Yke, p_Ykg, T, P, Bg, Bc, hw, p_Xw);
        cout << setw(10) << T << setw(15) << Bc << setw(15) << hw / 1.0e6;
        for (int i = 0; i < ns; ++i)
            cout << setw(15) << p_Xw[i] / (1.0 - p_Xw[1]);
        cout << endl;
    }
    
    delete [] p_Yke;
    delete [] p_Ykg;
    delete [] p_Xw;
}
