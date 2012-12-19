
#include "Transport.h"
#include "Constants.h"
#include "Numerics.h"

using namespace Numerics;

//==============================================================================

Transport::Transport(
    const Thermodynamics& thermo, const std::string& viscosity,
    const std::string& lambda)
    : m_thermo(thermo), m_collisions(thermo)
{ 
    // Load the viscosity calculator
    mp_viscosity = 
        Utilities::Factory<ViscosityAlgorithm>::create(viscosity, m_collisions);
    
    // Load the thermal conductivity calculator
    mp_thermal_conductivity = 
        Utilities::Factory<ThermalConductivityAlgorithm>::create(
            lambda, std::pair<Thermodynamics&, CollisionDB&>(
                const_cast<Thermodynamics&>(m_thermo), m_collisions));
    
    // Load the diffusion matrix calculator
    mp_diffusion_matrix =
        new Ramshaw(thermo, m_collisions);
    
    // Allocate work array storage
    mp_work1 = new double [m_thermo.nSpecies()];
    mp_work2 = new double [m_thermo.nSpecies()];
}
    
//==============================================================================
    
Transport::~Transport()
{
    delete mp_viscosity;
    delete mp_thermal_conductivity;
    delete mp_diffusion_matrix;
    
    delete [] mp_work1;
    delete [] mp_work2;
}

//==============================================================================

double Transport::reactiveThermalConductivity()
{
    // Compute dX_i/dT
    m_thermo.dXidT(m_thermo.T(), m_thermo.P(), m_thermo.X(), mp_work1);
    
    // Compute the multicomponent diffusion coefficient matrix
    const int ns = m_thermo.nSpecies();
    const RealMatrix& Dij = diffusionMatrix();
    const double* const Y = m_thermo.Y();
    
    /*cout << "\nDij" << endl;
    cout.precision(3);
    for (int i = 0; i < ns; ++i) {
        for (int j = 0; j < ns; ++j)
            cout << setw(11) << Dij(i,j);
        cout << endl;
    }
    cout << endl;*/
    
    for (int i = 0; i < ns; ++i) {
        mp_work2[i] = 0.0;
        for (int j = 0; j < ns; ++j)
            mp_work2[i] += Dij(i,j) * mp_work1[j];
    }
    
    // Compute the species enthalpies per unit mass
    m_thermo.speciesHOverRT(mp_work1);
    
    double lambda = 0.0;
    for (int i = 0; i < ns; ++i)
        lambda += mp_work1[i] * mp_work2[i] * Y[i] / m_thermo.speciesMw(i);
    
    return (lambda * RU * m_thermo.T() * m_thermo.density());
}

//==============================================================================

double Transport::sigma() 
{
    const double Th = m_thermo.T();
    const double Te = m_thermo.Te();
    const double nd = m_thermo.numberDensity();
    const double* const X = m_thermo.X();
    const double Mwe = m_thermo.speciesMw(0);

    const RealSymMat& Q11 = m_collisions.Q11(Th, Te, nd, X);
    const RealSymMat& Q22 = m_collisions.Q22(Th, Te, nd, X);
    const RealSymMat& B   = m_collisions.Bstar(Th, Te, nd, X);
    const RealVector& Cei = m_collisions.Cstei(Th, Te, nd, X);
    
    double fac = 0.0;
    double lam00 = 0.0;
    double lam01 = 0.0;
    double lam11 = std::sqrt(128.0) * X[0] * Q22(0);
    
    for (int i = 1; i < m_thermo.nSpecies(); ++i) {
        fac = X[i] * Q11(i);
        lam00 += fac * 8.0;
        lam01 += fac * (20.0 - 24.0 * Cei(i));
        lam11 += fac * (50.0 - 24.0 * B(i));
    }
    
    double beta = 3.0 * std::sqrt(TWOPI * RU * Te / Mwe);
    return X[0]*QE*QE*beta/(2.0*KB*Te*(lam00-lam01*lam01/lam11));
}

//==============================================================================

