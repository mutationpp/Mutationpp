
#include "mutation++.h"
#include <fenv.h>
#include <iostream>

using namespace std;
using namespace Mutation;
using namespace Mutation::Numerics;

double* p_Yke;
double* p_Ykg;
double* p_Ykc;
double Bg, T, P;

void computeF(Mixture& mix, RealVector& Xew, RealVector& F)
{
    const int ne = mix.nElements();
    const int ns = mix.nSpecies();
    
    double* p_Ykw = new double [ne];
    double* p_X   = new double [ns];
    double* p_Xe  = new double [ne];

    // Compute the equilibrium composition
    asVector(p_Xe, ne) = Xew(0,ne);
    mix.equilibrate(T, P, p_Xe, p_X);
    
    //cout << "Mole Fractions: " << endl;
    //for (int i = 0; i < ns; ++i)
    //    cout << setw(5) << mix.speciesName(i) << " " << p_X[i] << endl;
    //cout << "Sum: " << asVector(p_X, ns).sum() << endl;
    
    // Now compute the wall mass fractions
    double mwg = 0.0;
    
    for (int i = 0; i < ne; ++i)
        p_Ykw[i] = 0.0;;
    
    for (int j = 0; j < ns; ++j) {
        if (mix.species(j).phase() == Thermodynamics::GAS) {
            mwg += mix.speciesMw(j) * p_X[j];
            for (int i = 0; i < ne; ++i)
                p_Ykw[i] += mix.elementMatrix()(j,i) * p_X[j];
        }
    }
    
    for (int i = 0; i < ne; ++i)
        p_Ykw[i] *= mix.atomicMass(i) / mwg;
        
    cout << "sum(Ykw) = " << asVector(p_Ykw, ne).sum() << endl;
    //for (int i = 0; i < ne; ++i)
    //    cout << setw(5) << mix.elementName(i) << " " << p_Ykw[i] << endl;

    F(ne) = -1.0;
    double Bc = Xew(ne);
    for (int i = 0; i < ne; ++i) {
        F(i) = p_Ykw[i]*(1.0 + Bg + Bc) - (p_Yke[i] + Bg*p_Ykg[i] + Bc*p_Ykc[i]);
        F(ne) += Xew(i);
    }
    
    delete [] p_Ykw;
    delete [] p_X;
    delete [] p_Xe;
}

void computeJwFD(Mixture& mix, RealVector& Xew, RealMatrix& J)
{
    const int ne = mix.nElements();
    
    RealVector F0(ne+1);
    RealVector Fp(ne+1);
    RealVector Xp(ne+1);
    
    double h;
    
    computeF(mix, Xew, F0);
    
    for (int i = 0; i < ne+1; ++i) {
        h = Xew(i)*std::sqrt(NumConst<double>::eps);
        Xp = Xew;
        Xp(i) += h;
        computeF(mix, Xp, Fp);
        for (int j = 0; j < ne+1; ++j) {
            J(j,i) = (Fp(j)-F0(j))/h;
        }
    }
}


int main()
{
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

    Mixture mix("tacot24");
    
    const int ne = mix.nElements();
    const int ns = mix.nSpecies();
    
    p_Yke = new double [ne];
    p_Ykg = new double [ne];
    p_Ykc = new double [ne];
    
    double* p_Yw = new double [ne];
    double* p_Xw = new double [ne];
    double* p_X = new double [ns];
    
    RealVector F(ne+1);
    RealVector Xew(ne+1);
    RealMatrix J(ne+1,ne+1);
    RealVector dX(ne+1);
    
    // Run conditions
    Bg = 0.6;
    //T  = 2000.0;
    P  = 101325.0;
    
    p_Yke[mix.elementIndex("C")] = 0.0;
    p_Yke[mix.elementIndex("H")] = 0.0;
    p_Yke[mix.elementIndex("O")] = 0.21;
    p_Yke[mix.elementIndex("N")] = 0.79;
    mix.convert<Thermodynamics::XE_TO_YE>(p_Yke, p_Yke);
    cout << "Boundary Layer Gas" << endl;
    for (int i = 0; i < ne; ++i)
        cout << setw(10) << "Y_" + mix.elementName(i) << setw(10) << p_Yke[i] << endl;
    
    p_Ykg[mix.elementIndex("C")] = 0.206;
    p_Ykg[mix.elementIndex("H")] = 0.679;
    p_Ykg[mix.elementIndex("O")] = 0.115;
    p_Ykg[mix.elementIndex("N")] = 0.0;
    mix.convert<Thermodynamics::XE_TO_YE>(p_Ykg, p_Ykg);
    cout << endl << "Pyrolysis Gas" << endl;
    for (int i = 0; i < ne; ++i)
        cout << setw(10) << "Y_" + mix.elementName(i) << setw(10) << p_Ykg[i] << endl;
    
    p_Ykc[mix.elementIndex("C")] = 1.0;
    p_Ykc[mix.elementIndex("H")] = 0.0;
    p_Ykc[mix.elementIndex("O")] = 0.0;
    p_Ykc[mix.elementIndex("N")] = 0.0;
    
    cout << endl << setw(10) << "Tw (K)" << setw(15) << "B'c" << setw(15) << "hw (MJ/kg)" << endl;
    
    for (double T = 300.0; T < 3800.1; T += 50.0) {
    
    // Initialize the wall element fractions to be the pyrolysis gas fractions
    double sum = 0.0;
    for (int i = 0; i < ne; ++i) {
        p_Yw[i] = p_Yke[i] + Bg*p_Ykg[i] + 100.0*Bg*p_Ykc[i];
        sum += p_Yw[i];
    }
    
    for (int i = 0; i < ne; ++i)
        p_Yw[i] /= sum;
        
    mix.convert<Thermodynamics::YE_TO_XE>(p_Yw, p_Xw);
    // cout << endl << "Wall Composition" << endl;
    //for (int i = 0; i < ne; ++i)
    //    cout << setw(10) << "Y_" + mix.elementName(i) << setw(10) << p_Yw[i] << endl;
    
    mix.equilibrate(T, P, p_Xw, p_X);
    //cout << "Mole Fractions" << endl;
    //for (int i = 0; i < ns; ++i)
    //    cout << setw(10) << mix.speciesName(i) << setw(15) << p_X[i] << endl;
    
    // compute the Ykw
    double mwg = 0.0;
    
    for (int i = 0; i < ne; ++i)
        p_Yw[i] = 0.0;;
    
    for (int j = 0; j < ns; ++j) {
        if (mix.species(j).phase() == Thermodynamics::GAS) {
            mwg += mix.speciesMw(j) * p_X[j];
            for (int i = 0; i < ne; ++i)
                p_Yw[i] += mix.elementMatrix()(j,i) * p_X[j];
        }
    }
    
    for (int i = 0; i < ne; ++i)
        p_Yw[i] *= mix.atomicMass(i) / mwg;
    
    double Bc = (p_Yke[0] + Bg*p_Ykg[0] - p_Yw[0]*(1.0 + Bg)) / (p_Yw[0] - p_Ykc[0]);
    cout << setw(10) << T << setw(15) << Bc << setw(15) << mix.mixtureHMass() / 1.0e6 << endl;
    
    
    }
    
    // Initialize the B'c value to 1
    /*Xew(ne) = 1.0;
    
    for (int i = 0; i < 10; ++i) {    
        // Now compute the F vector
        cout << "X = " << endl;
        cout << Xew << endl;
        computeF(mix, Xew, F);
        
        // Compute the jacobian
        computeJwFD(mix, Xew, J);
        
        QRP<double> qrp(J);
        F *= -1.0;
        qrp.solve(dX, F);
        
        // Output initial conditions
        cout << "F = " << endl;
        cout << F << endl;
        cout << "norm(F) = " << F.norm2() << endl;
        cout << "J = " << endl;
        cout << J << endl;
        cout << "dX = " << endl;
        cout << dX << endl;
        
        Xew += dX*0.01;
    }*/
    
    delete [] p_Yke;
    delete [] p_Ykc;
    delete [] p_Ykg;
    
    delete [] p_Yw;
    delete [] p_Xw;
    delete [] p_X;
}
