#ifndef TRANSFER_TRANSFER_MODEL_H
#define TRANSFER_TRANSFER_MODEL_H

#include "Thermodynamics.h"
#include "MillikanWhite.h"

namespace Mutation {
    namespace Transfer {

class TransferModel
{
public:
    virtual ~TransferModel() { }
    virtual void source(double* const p_source) = 0;
};

/**
 * Represents a coupling between vibrational and translational energy modes.
 */
class TransferModelVT : public TransferModel
{
public:

    TransferModelVT(const Thermodynamics::Thermodynamics& thermo)
        : m_mw(thermo)
    { }
    
    virtual void source(double* const p_source)
    {
    
    }

private:

    MillikanWhite m_mw;
};

/**
 * Calculation of omegaTE
 */
class TransferModelTE : public TransferModel
{
public:

    TransferModelTE(const Thermodynamics::Thermodynamics& thermo, Transport::Transport& transport)
        : m_thermo(thermo), 
          m_collisions(transport.collisionDB()), 
          m_iar(m_thermo.speciesIndex("Ar(30)")), 
          m_iarp(m_thermo.speciesIndex("Ar+")),
          m_iem(m_thermo.speciesIndex("e-")),
          m_ns(m_thermo.nSpecies()) 
                         
    { }

virtual void source(double* const p_source)
{
    const double T = m_thermo.T();
    const double Te = m_thermo.Te();
    const double mm = m_thermo.mixtureMw();
    const double nd = m_thermo.numberDensity();    
    const double mm_Ar = m_thermo.speciesMw(m_iar);
    const double mm_Arp = m_thermo.speciesMw(m_iarp);
    const double mm_el = m_thermo.speciesMw(0);
    
 
    // Collisional integrals
    const double *const p_x=m_thermo.X();
    const Numerics::RealSymMat& Q11 =m_collisions.Q11(T, Te, nd, p_x);
    const double Q11_Arp = Q11(m_iarp);
    const double Q11_Ar = Q11(m_iar);

 
    // Mass fractions of argon and argon+
    const double* const p_y = m_thermo.Y();
    double y_Arp = p_y[m_iarp];
    double y_Ar = 1.0- y_Arp - p_y[0];


    // Densities
    double rho[m_thermo.nSpecies()];
    for (int i=0; i< m_ns; i++)
        rho[i] =  p_y[i]* m_thermo.density();   

    // Number density charged
    const double nd_plus = rho[m_iem]*NA/ mm_el;

    // Calculation of omegaTE
    double nu_Arp = sqrt(RU*8.0*Te/(PI*mm_el))*nd_plus*2.0*Q11_Arp;
    double nu_Ar = sqrt(RU*8.0*Te/(PI*mm_el))*nd*2.0*Q11_Ar;

    double sum_int = nu_Ar/mm_Ar + nu_Arp/mm_Arp;
     
    double omegaTE = 3.0*RU*rho[m_iem]*(T-Te)*sum_int;

    // Calculation of omegaEexc (electron impact excitation)
    double omegaEexc = 0.0;

    // Calculation of omegaEI (electron impact ionization)
    double omegaEI = 0.0;

    p_source[0] = omegaTE + omegaEI + omegaEexc;

 }

private:

    
    const Mutation::Thermodynamics::Thermodynamics& m_thermo;
    Transport::CollisionDB& m_collisions;
    int m_iar;
    int m_iarp;
    int m_ns;
    int m_iem;

};


    } // namespace Transfer
} // namespace Mutation

#endif // TRANSFER_TRANSFER_MODEL_H
