#ifndef TRANSPORT_TRANSPORT_H
#define TRANSPORT_TRANSPORT_H

#include "Thermodynamics.h"
#include "CollisionDB.h"
#include "ViscosityAlgorithm.h"
#include "ThermalConductivityAlgorithm.h"
#include "Ramshaw.h"
#include "Utilities.h"

namespace Mutation {
    namespace Transport {

/**
 * Manages the computation of transport properties.
 */
class Transport
{
public:
    
    /**
     * Constructs a Transport object given a Thermodynamics reference.
     */
    Transport(
        Mutation::Thermodynamics::Thermodynamics& thermo, 
        const std::string& viscosity, const std::string& lambda);
    
    /**
     * Destructor.
     */
    ~Transport();
    
    /**
     * Sets the algorithm to use when computing viscosity.
     */
    void setViscosityAlgo(const std::string& algo) {
        delete mp_viscosity;
        mp_viscosity = Mutation::Utilities::Config::Factory<
            ViscosityAlgorithm>::create(algo, m_collisions);
    }
    
    /**
     * Sets the algorithm to use when computing thermal conductivity.
     */
    void setThermalConductivityAlgo(const std::string& algo) {
        delete mp_thermal_conductivity;
        mp_thermal_conductivity = Mutation::Utilities::Config::Factory<
            ThermalConductivityAlgorithm>::create(algo, m_collisions);
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
    double viscosity() {
        return mp_viscosity->viscosity(
            m_thermo.T(), m_thermo.numberDensity(), m_thermo.X());
    }
    
    /**
     * Returns the mixture thermal conductivity for a frozen mixture.
     */
    double frozenThermalConductivity() {
        return (
            heavyThermalConductivity() + 
            electronThermalConductivity() +
            internalThermalConductivity()
        );
    }
    
    /**
     * Returns the mixture thermal conductivity for a mixture in thermochemical
     * equilibrium.
     */
    double equilibriumThermalConductivity() {
        return (
            frozenThermalConductivity() +
            reactiveThermalConductivity()
        );
    }
    
    /**
     * Returns the heavy particle translational thermal conductivity using the 
     * set algorithm.
     */
    double heavyThermalConductivity() {
        return mp_thermal_conductivity->thermalConductivity(
            m_thermo.T(), m_thermo.Te(), m_thermo.numberDensity(), m_thermo.X());
    }
    
    /**
     * Returns the electron translational thermal conductivity.
     */
    double electronThermalConductivity();
    
    /**
     * Returns the internal energy thermal conductivity using Euken's formulas.
     */
    double internalThermalConductivity();
    
    /**
     * Returns the reactive thermal conductivity which accounts for reactions
     * for mixtures in thermochemical equilibrium.
     */
    double reactiveThermalConductivity();
    
    /**
     * Returns the binary diffusion coefficient matrix.
     */
    const Mutation::Numerics::RealMatrix& diffusionMatrix() {
        return mp_diffusion_matrix->diffusionMatrix(
            m_thermo.T(), m_thermo.numberDensity(), m_thermo.X());
    }
    
    /**
     * Returns the electric conductivity.
     */
    double sigma();
    
private:

    Mutation::Thermodynamics::Thermodynamics& m_thermo;
    CollisionDB m_collisions;
    
    ViscosityAlgorithm* mp_viscosity;
    ThermalConductivityAlgorithm* mp_thermal_conductivity;
    Ramshaw* mp_diffusion_matrix;
    
    double* mp_wrk1;
    double* mp_wrk2;
    double* mp_wrk3;
};

    } // namespace Transport
} // namespace Mutation

#endif // TRANSPORT_TRANSPORT_H
