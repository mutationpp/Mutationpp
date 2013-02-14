
#include "Transport.h"
#include "Constants.h"
#include "Numerics.h"

using namespace Mutation::Numerics;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities;

//==============================================================================

Transport::Transport(
    Thermodynamics& thermo, const std::string& viscosity,
    const std::string& lambda)
    : m_thermo(thermo), m_collisions(thermo)
{ 
    // Load the viscosity calculator
    mp_viscosity = 
        Config::Factory<ViscosityAlgorithm>::create(viscosity, m_collisions);
    
    // Load the thermal conductivity calculator
    mp_thermal_conductivity = 
        Config::Factory<ThermalConductivityAlgorithm>::create(
            lambda, m_collisions);
    
    // Load the diffusion matrix calculator
    mp_diffusion_matrix =
        new Ramshaw(thermo, m_collisions);
    
    // Allocate work array storage
    mp_wrk1 = new double [m_thermo.nSpecies()];
    mp_wrk2 = new double [m_thermo.nSpecies()];
    mp_wrk3 = new double [m_thermo.nSpecies()];
}
    
//==============================================================================
    
Transport::~Transport()
{
    delete mp_viscosity;
    delete mp_thermal_conductivity;
    delete mp_diffusion_matrix;
    
    delete [] mp_wrk1;
    delete [] mp_wrk2;
    delete [] mp_wrk3;
}

//==============================================================================

double Transport::electronThermalConductivity()
{
    if (!m_thermo.hasElectrons() || m_thermo.X()[0] < 1.0e-30) 
        return 0.0;

    // Get thermodynamic properties
    const int ns     = m_thermo.nSpecies();
    const double Th  = m_thermo.T();
    const double Te  = m_thermo.Te();
    const double nd  = m_thermo.numberDensity();
    const double me  = m_thermo.speciesMw(0)/NA;
    const double *const X = m_thermo.X();
    
    // Get collision integral information
    const RealSymMat& Q11   = m_collisions.Q11(Th, Te, nd, X);
    const RealSymMat& Q22   = m_collisions.Q22(Th, Te, nd, X);
    const RealSymMat& B     = m_collisions.Bstar(Th, Te, nd, X);
    const RealVector& Q12ei = m_collisions.Q12ei(Th, Te, nd, X);
    const RealVector& Q13ei = m_collisions.Q13ei(Th, Te, nd, X);
    const RealVector& Q14ei = m_collisions.Q14ei(Th, Te, nd, X);
    const RealVector& Q15ei = m_collisions.Q15ei(Th, Te, nd, X);
    const double      Q23ee = m_collisions.Q23ee(Th, Te, nd, X);
    const double      Q24ee = m_collisions.Q24ee(Th, Te, nd, X);
    
    // Compute the lambdas
    double fac;
    double lam11 = 0.0;
    double lam12 = 0.0;
    double lam22 = 0.0;
    
    for (int i = 1; i < ns; ++i) {
        lam11 += X[i]*Q11(i)*(25.0/4.0-3.0*B(i));
        lam12 += X[i]*(175.0/16.0*Q11(i)-315.0/8.0*Q12ei(i)+57.0*Q13ei(i)-30.0*
            Q14ei(i));
        lam22 += X[i]*(1225.0/64.0*Q11(i)-735.0/8.0*Q12ei(i)+399.0/2.0*Q13ei(i)-
            210.0*Q14ei(i)+90.0*Q15ei(i));
    }
    
    fac = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0];
    lam11 = fac*(lam11 + SQRT2*X[0]*Q22(0));
    lam12 = fac*(lam12 + SQRT2*X[0]*(7.0/4.0*Q22(0)-2.0*Q23ee));
    lam22 = fac*(lam22 + SQRT2*X[0]*(77.0/16.0*Q22(0)-7.0*Q23ee+5.0*Q24ee));
    
    // 2nd order solution
    //return (X[0]*X[0]/lam11);
    
    // 3rd order solution
    return (X[0]*X[0]*lam22/(lam11*lam22-lam12*lam12));
}

//==============================================================================

double Transport::internalThermalConductivity()
{
    const int ns     = m_thermo.nSpecies();
    const double Th  = m_thermo.T();
    const double Te  = m_thermo.Te();
    const double Tr  = m_thermo.Tr();
    const double Tv  = m_thermo.Tv();
    const double Tel = m_thermo.Tel();
    const double nd  = m_thermo.numberDensity(); 
    const double *const X = m_thermo.X();
    
    const RealSymMat& nDij = m_collisions.nDij(Th, Te, nd, X);
    
    m_thermo.speciesCpOverR(
        Th, Te, Tr, Tv, Tel, NULL, NULL, mp_wrk1, mp_wrk2, mp_wrk3);
    
    double lambda = 0.0;
    for (int i = 1; i < ns; ++i) {
        double sum = 0.0;
        for (int j = 0; j < ns; ++j)
            sum += X[j] / nDij(i,j);
        
        lambda += X[i] * (mp_wrk1[i] + mp_wrk2[i] + mp_wrk3[i]) / sum;
    }
    
    return (lambda * KB);
}

//==============================================================================

double Transport::reactiveThermalConductivity()
{
    // Compute dX_i/dT
    m_thermo.dXidT(m_thermo.T(), m_thermo.P(), m_thermo.X(), mp_wrk1);
    
    // Compute the multicomponent diffusion coefficient matrix
    const int ns = m_thermo.nSpecies();
    const RealMatrix& Dij = diffusionMatrix();
    const double* const X = m_thermo.X();
    const double* const Y = m_thermo.Y();
    
    // Store the species charges and determine the mixture charge
    double q = 0.0;
    for (int i = 0; i < ns; ++i) {
        mp_wrk3[i] = m_thermo.speciesCharge(i);
        q += mp_wrk3[i] * X[i];
    }
    
    /*cout << "\nDij" << endl;
    cout.precision(3);
    for (int i = 0; i < ns; ++i) {
        for (int j = 0; j < ns; ++j)
            cout << setw(11) << Dij(i,j);
        cout << endl;
    }
    cout << endl;*/
    
    for (int i = 0; i < ns; ++i) {
        mp_wrk2[i] = 0.0;
        for (int j = 0; j < ns; ++j) {
            mp_wrk2[i] += Dij(i,j)*(mp_wrk1[j]-(X[j]*mp_wrk3[j]-Y[j]*q)/
                (X[0]*mp_wrk3[0]-Y[0]*q)*mp_wrk1[0]);
        }
    }
    
    // Compute the species enthalpies per unit mass
    m_thermo.speciesHOverRT(mp_wrk1);
    
    double lambda = 0.0;
    for (int i = 0; i < ns; ++i)
        lambda += mp_wrk1[i] * mp_wrk2[i] * X[i];
    
    return (m_thermo.P()*lambda);
}

//==============================================================================

double Transport::sigma() 
{
    if (!m_thermo.hasElectrons() || m_thermo.X()[0] < 1.0e-30) 
        return 0.0;

    const double ns = m_thermo.nSpecies();
    const double Th = m_thermo.T();
    const double Te = m_thermo.Te();
    const double nd = m_thermo.numberDensity();
    const double* const X = m_thermo.X();
    const double me = m_thermo.speciesMw(0)/NA;

    const RealSymMat& Q11 = m_collisions.Q11(Th, Te, nd, X);
    const RealSymMat& Q22 = m_collisions.Q22(Th, Te, nd, X);
    const RealSymMat& B   = m_collisions.Bstar(Th, Te, nd, X);
    const RealVector& Cei = m_collisions.Cstei(Th, Te, nd, X);
    
    // Compute lambdas
    double lam00 = 0.0;
    double lam01 = 0.0;
    double lam11 = 0.0;
    double fac   = 0.0;
    
    for (int i = 1; i < ns; ++i) {
        fac = X[i]*Q11(i);
        lam00 += fac;
        lam01 += fac*(2.5-3.0*Cei(i));
        lam11 += fac*(25.0/4.0-3.0*B(i));
    }
    
    fac = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0];
    lam00 = fac*lam00;
    lam01 = fac*lam01;
    lam11 = fac*(lam11 + SQRT2*X[0]*Q22(0));
    
    // First order
    //return (4.0/25.0*(X[0]*X[0]*QE*QE)/(KB*KB*Te*lam00));
    
    // Second order
    return (4.0/25.0*(X[0]*X[0]*QE*QE)/(KB*KB*Te*(lam00-lam01*lam01/lam11)));
}

//==============================================================================

