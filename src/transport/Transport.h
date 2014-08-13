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
class Transport //: public Mutation::Thermodynamics::StateModelUpdateHandler
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
    
    void stateUpdated() {
        std::cout << "stateUpdated: Transport" << std::endl;
    }
    
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
    
    void omega11ii(double* const p_omega);
    void omega22ii(double* const p_omega);
    
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
            reactiveThermalConductivity() +
            soretThermalConductivity()
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
     * Returns the Soret thermal conductivity the mixture in thermochemical 
     * equilibrium.
     */
    double soretThermalConductivity();
    
    /**
     * Returns the thermal diffusion ratios for each species.
     */
    void thermalDiffusionRatios(double* const p_k) {
        return mp_thermal_conductivity->thermalDiffusionRatios(
            m_thermo.T(), m_thermo.Te(), m_thermo.numberDensity(),
            m_thermo.X(), p_k);
    }
    
    /**
     * Returns the multicomponent diffusion coefficient matrix.
     */
    const Mutation::Numerics::RealMatrix& diffusionMatrix() {
        return mp_diffusion_matrix->diffusionMatrix(
            m_thermo.T(), m_thermo.numberDensity(), m_thermo.X());
    }
    
    /**
     * Returns mixture averaged species diffusion coefficients which are defined
     * as \f[ D_i = \frac{1 - X_i}{\sum_j X_j / \mathcal{D}_{ij}} \f] where 
     * \f$\mathcal{D}_{ij}\f$ are the binary diffusion coefficients for each
     * species pair.
     */
    void averageDiffusionCoeffs(double *const p_Di);
    
    void equilibriumFickP(double* const p_F);
    void equilibriumFickT(double* const p_F);
    void equilibriumFickXe(double* const p_F);

    /**
     * Computes the species diffusion velocities and ambipolar electric field 
     * using the Ramshaw approximation of the generalized Stefan-Maxwell 
     * equations and the supplied modified driving forces
     *
     * \f[ d^{'}_i = \frac{\nabla p_i}{n k_B T_h} - \frac{y_i p}{n k_B T_h}
     *     \nabla \ln p + k^h_{Ti} \nabla \ln T_h + k^e_{Ti} \frac{T_h}{T_e}
     *     \nabla \ln T_e \f]
     *
     * @param p_dp - the vector of modified driving forces \f$ d^{'}_i \f$
     * @param p_V  - on return, the vector of diffusion velocities
     * @param E    - on return, the ambipolar electric field
     */
    void stefanMaxwell(const double* const p_dp, double* const p_V, double& E);
        
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
