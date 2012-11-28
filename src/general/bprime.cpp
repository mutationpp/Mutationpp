
#include "mutation++.h"
#include "BPrimeProvider.h"

#include <fenv.h>

int main()
{
    feenableexcept(FE_INVALID|FE_DIVBYZERO);

    Mixture mix("phenolJDM");
    BPrimeProvider bprime(mix);
    
    const int ne = mix.nElements();
    const int ns = mix.nSpecies();
    
    double* p_Yke = new double [ne];
    double* p_Ykg = new double [ne];
    double* p_X   = new double [ns];
    
    double Bg, Bc, T, P;
    Bg = 10.0;
    T  = 3000.0;
    P  = 101325.0;
    
    p_Yke[mix.elementIndex("C")] = 0.0;
    p_Yke[mix.elementIndex("H")] = 0.0;
    p_Yke[mix.elementIndex("O")] = 0.29;
    p_Yke[mix.elementIndex("N")] = 0.71;
    
    p_Ykg[mix.elementIndex("C")] = 0.375;
    p_Ykg[mix.elementIndex("H")] = 0.162;
    p_Ykg[mix.elementIndex("O")] = 0.463;
    p_Ykg[mix.elementIndex("N")] = 0.0;
    
    bprime.massBalance(p_Yke, p_Ykg, Bg, T, P, Bc, p_X);
    
    cout << "Bc: " << Bc << endl;
    
    delete [] p_Yke;
    delete [] p_Ykg;
    delete [] p_X;
}
