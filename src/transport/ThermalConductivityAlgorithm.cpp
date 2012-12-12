
#include "ThermalConductivityAlgorithm.h"
#include "Numerics.h"
#include "Thermodynamics.h"
#include "CollisionDB.h"
#include "Constants.h"

//==============================================================================

ThermalConductivityAlgorithm::ThermalConductivityAlgorithm(ARGS arguments)
    : m_thermo(arguments.first), m_collisions(arguments.second)
{ 
    mp_cvr  = new double [m_thermo.nSpecies()];
    mp_cvv  = new double [m_thermo.nSpecies()];
    mp_cvel = new double [m_thermo.nSpecies()];
}

//==============================================================================

ThermalConductivityAlgorithm::~ThermalConductivityAlgorithm() 
{ 
    delete [] mp_cvr;
    delete [] mp_cvv;
    delete [] mp_cvel;
}

//==============================================================================

double ThermalConductivityAlgorithm::thermalConductivity(
    double Th, double Te, double Tr, double Tv, double Tel, double nd, 
    const double* const p_x)
{
    const int ns = m_thermo.nSpecies();
    
    // Thermal conductivity of internal energy using Euken's formula
    double lambda = 0.0;
    
    // Now apply Euken's correction
    const Numerics::RealSymMat& nDij  = m_collisions.nDij(Th, Te, nd, p_x);
    
    m_thermo.speciesCpOverR(
        Th, Te, Tr, Tv, Tel, NULL, NULL, mp_cvr, mp_cvv, mp_cvel);
    
    for (int i = 1; i < ns; ++i) {
        double sum = 0.0;
        for (int j = 0; j < ns; ++j)
            sum += p_x[j] / nDij(i,j);
        
        lambda += p_x[i] * (mp_cvr[i] + mp_cvv[i] + mp_cvel[i]) / sum;
    }
    
    return (lambda * KB + elasticThermalConductivity(Th, Te, nd, p_x));
}

//==============================================================================
