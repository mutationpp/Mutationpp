#ifndef TRANSPORT_TRANSPORT_H
#define TRANSPORT_TRANSPORT_H

#include "Thermodynamics.h"
#include "CollisionDB.h"
#include "ViscosityAlgorithm.h"
#include "ThermalConductivityAlgorithm.h"
#include "Ramshaw.h"
#include "AutoRegistration.h"
#include "Constants.h"

/**
 * Manages the computation of transport properties.  The class is templated to 
 * allow different transport algorithms to be used (chosen at compile time) with
 * a central instantiation of the collision database.
 */
class Transport
{
public:
    
    /**
     * Constructs a Transport object given a Thermodynamics reference.
     */
    Transport(
        const Thermodynamics& thermo, const std::string& viscosity,
        const std::string& lambda)
        : m_thermo(thermo), m_collisions(thermo)
    { 
        mp_viscosity = 
            Utilities::Factory<ViscosityAlgorithm>::create(viscosity, m_collisions);
        mp_thermal_conductivity = 
            Utilities::Factory<ThermalConductivityAlgorithm>::create(
                lambda, std::pair<Thermodynamics&, CollisionDB&>(
                    const_cast<Thermodynamics&>(m_thermo), m_collisions));
        mp_diffusion_matrix =
            new Ramshaw(thermo, m_collisions);
    }
    
    ~Transport()
    {
        delete mp_viscosity;
        delete mp_thermal_conductivity;
        delete mp_diffusion_matrix;
    }
    
    /**
     * Sets the algorithm to use when computing viscosity.
     */
    void setViscosityAlgo(const std::string& algo) {
        delete mp_viscosity;
        mp_viscosity = Utilities::Factory<ViscosityAlgorithm>::create(algo, m_collisions);
    }
    
    /**
     * Sets the algorithm to use when computing thermal conductivity.
     */
    void setThermalConductivityAlgo(const std::string& algo) {
        delete mp_thermal_conductivity;
        mp_thermal_conductivity = 
            Utilities::Factory<ThermalConductivityAlgorithm>::create(
                algo, std::pair<Thermodynamics&, CollisionDB&>(
                    const_cast<Thermodynamics&>(m_thermo), m_collisions));
    }
    
    /**
     * Returns the total number of collision pairs accounted for in the 
     * collision database.
     */
    int nCollisionPairs() const {
        return m_collisions.nCollisionPairs();
    }
    
    /**
     * Returns the mixture viscosity.
     */
    double eta(const double T, const double nd, double *const p_x) {
        return mp_viscosity->viscosity(T, nd, p_x);
    }
    
    /**
     * Returns the mixture viscosity.
     */
    double eta() {
        return mp_viscosity->viscosity(
            m_thermo.T(), m_thermo.numberDensity(), m_thermo.X());
    }
    
    /**
     * Returns the mixture translational thermal conductivity.
     */
    double lambda(const double T, const double nd, double *const p_x) {
        return mp_thermal_conductivity->thermalConductivity(
            T, T, T, T, T, nd, p_x);
    }
    
    /**
     * Returns the mixture translational thermal conductivity.
     */
    double lambda() {
        return mp_thermal_conductivity->thermalConductivity(
            m_thermo.T(), m_thermo.Te(), m_thermo.Tr(), m_thermo.Tv(),
            m_thermo.Tel(), m_thermo.numberDensity(), m_thermo.X());
    }
    
    /**
     * Returns the binary diffusion coefficient matrix.
     */
    const RealMatrix& diffusionMatrix(
        const double T, const double nd, const double *const p_x) {
        return mp_diffusion_matrix->diffusionMatrix(T, nd, p_x);
    }
    
    /**
     * Returns the binary diffusion coefficient matrix.
     */
    const RealMatrix& diffusionMatrix() {
        return mp_diffusion_matrix->diffusionMatrix(
            m_thermo.T(), m_thermo.numberDensity(), m_thermo.X());
    }
    
    /**
     * Returns the electric conductivity.
     */
    double sigma() 
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
        
        double fac = 0;
        double lam00 = 0;
        double lam01 = 0;
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
    
private:

    const Thermodynamics& m_thermo;
    CollisionDB m_collisions;
    
    ViscosityAlgorithm* mp_viscosity;
    ThermalConductivityAlgorithm* mp_thermal_conductivity;
    Ramshaw* mp_diffusion_matrix;
};

#endif // TRANSPORT_TRANSPORT_H
