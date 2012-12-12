#ifndef TRANSPORT_TRANSPORT_H
#define TRANSPORT_TRANSPORT_H

#include "Thermodynamics.h"
#include "CollisionDB.h"
#include "ViscosityAlgorithm.h"
#include "ThermalConductivityAlgorithm.h"
#include "Ramshaw.h"
#include "AutoRegistration.h"

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
    
private:

    const Thermodynamics& m_thermo;
    CollisionDB m_collisions;
    
    ViscosityAlgorithm* mp_viscosity;
    ThermalConductivityAlgorithm* mp_thermal_conductivity;
    Ramshaw* mp_diffusion_matrix;
};

#endif // TRANSPORT_TRANSPORT_H
