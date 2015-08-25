/**
 * @file Transport.h
 *
 * @brief Declaration of Transport class.
 */

/*
 * Copyright 2014 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Mutation++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Mutation++.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef TRANSPORT_TRANSPORT_H
#define TRANSPORT_TRANSPORT_H

#include "Thermodynamics.h"
#include "CollisionDB.h"
#include "ViscosityAlgorithm.h"
#include "ThermalConductivityAlgorithm.h"
#include "Ramshaw.h"
#include "Utilities.h"

#include <Eigen/Dense>

namespace Mutation {
    namespace Transport {

#define ERROR_IF_INTEGRALS_ARE_NOT_LOADED(__RET__)\
if (mp_collisions == NULL) {\
    std::cout << "Error! Trying to use transport without loading collision integrals!!" << std::endl;\
    return __RET__;\
}

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
        const std::string& viscosity, const std::string& lambda,
        const bool load_data = true);
    
    /**
     * Destructor.
     */
    ~Transport();

    /**
    * Returns a pointer to the collision data object.
    */
    CollisionDB* const collisionData() {
        return mp_collisions;
    }
    
    void stateUpdated() {
        std::cout << "stateUpdated: Transport" << std::endl;
    }
    
    /**
     * Sets the algorithm to use when computing viscosity.
     */
    void setViscosityAlgo(const std::string& algo) {
        delete mp_viscosity;
        mp_viscosity = Mutation::Utilities::Config::Factory<
            ViscosityAlgorithm>::create(algo, *mp_collisions);
    }
    
    /**
     * Sets the algorithm to use when computing thermal conductivity.
     */
    void setThermalConductivityAlgo(const std::string& algo) {
        delete mp_thermal_conductivity;
        mp_thermal_conductivity = Mutation::Utilities::Config::Factory<
            ThermalConductivityAlgorithm>::create(algo, *mp_collisions);
    }
    
    /**
     * Returns the total number of collision pairs accounted for in the 
     * collision database.
     */
    int nCollisionPairs() const {
        ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0)
        return mp_collisions->nCollisionPairs();
    }
    
    void omega11ii(double* const p_omega);
    void omega22ii(double* const p_omega);
    
    /**
     * Returns the mixture viscosity.
     */
    double viscosity() {
        ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)
        return mp_viscosity->viscosity(
            m_thermo.T(), m_thermo.numberDensity(), m_thermo.X());
    }
    
    /**
     * Returns the mixture thermal conductivity for a frozen mixture.
     * To be used only at thermal equilibrium.
     */
    double frozenThermalConductivity() {
        ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)
        return (
            heavyThermalConductivity() + 
            electronThermalConductivity() +
            internalThermalConductivity()
        );
    }

    /**
     * Returns the mixture thermal conductivity vector for a frozen mixture
     * according to the state model.
     */
    void frozenThermalConductivityVector(double* const p_lambda);
    
    /**
     * Returns the mixture thermal conductivity for a mixture in thermochemical
     * equilibrium.
     */
    double equilibriumThermalConductivity() {
        ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)
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
        ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)
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
     * Returns the rotational energy thermal conductivity using Euken's formulas.
     */
    double rotationalThermalConductivity();

    /**
     * Returns the vibrational energy thermal conductivity using Euken's formulas.
     */
    double vibrationalThermalConductivity();

    /**
     * Returns the electronic energy thermal conductivity using Euken's formulas.
     */
    double electronicThermalConductivity();
    
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
        ERROR_IF_INTEGRALS_ARE_NOT_LOADED()
        return mp_thermal_conductivity->thermalDiffusionRatios(
            m_thermo.T(), m_thermo.Te(), m_thermo.numberDensity(),
            m_thermo.X(), p_k);
    }
    
    /**
     * Returns the multicomponent diffusion coefficient matrix.
     */
    const Eigen::MatrixXd& diffusionMatrix() {
        static Eigen::MatrixXd empty;
        ERROR_IF_INTEGRALS_ARE_NOT_LOADED(empty)
        return mp_diffusion_matrix->diffusionMatrix(
            m_thermo.T(), m_thermo.Te(), m_thermo.numberDensity(), m_thermo.X());
    }
    
    /**
     * Returns mixture averaged species diffusion coefficients which are defined
     * as \f[ D_i = \frac{1 - X_i}{\sum_j X_j / \mathcal{D}_{ij}} \f] where 
     * \f$\mathcal{D}_{ij}\f$ are the binary diffusion coefficients for each
     * species pair.
     */
    void averageDiffusionCoeffs(double *const p_Di);
    
    /**
     * Computes the vector of elemental mass and energy diffusion fluxes per
     * pressure gradient for a mixture in equilibrium, \f$\mathcal{F}_k^p\f$ and
     * \f$\mathcal{F}_h^p\f$,
     * where
     * \f[
     * \mathcal{F}_k^p = -\sum_{i\in\mathcal{S}} \rho_i \frac{M_{w,k}}{M_{w,i}}
     * 	   \nu_i^k \sum_{j\in\mathcal{S}} D_{ij} \left(
     * 	   \frac{\partial x_j}{\partial p} +
     * 	   \frac{x_j - y_j}{p} \right), \quad \forall\;k\in\mathcal{E},
     * \f]
     * and
     * \f[
     * \mathcal{F}_h^p = -\sum_{i\in\mathcal{S}} \rho_i h_i
     * 	   \sum_{j\in\mathcal{S}} D_{ij} \left(
     * 	   \frac{\partial x_j}{\partial p} +
     * 	   \frac{x_j - y_j}{p} \right).
     * \f]
     *
     * The total mass diffusion flux of element \f$k\f$ for a mixture at
     * equilibrium can then be computed with
     * \f[
     * \rho_k V_k = \mathcal{F}_k^p \nabla p + \mathcal{F}_k^T \nabla T +
     *     \sum_{l\in\mathcal{E}} \mathcal{F}_k^{z_l} \nabla z_l,
     *     \quad \forall\;k\in\mathcal{E}.
     * \f]
     * The total energy diffusion flux is computed by
     * \f[
     * \sum_{i\in\mathcal{S}} \rho_i h_i V_i = \mathcal{F}_h^p \nabla p +
     * 	   \mathcal{F}_h^T \nabla T +
     *     \sum_{l\in\mathcal{E}} \mathcal{F}_h^{z_l} \nabla z_l.
     * \f]
     *
     * @param p_F - A vector with length greater than the number of elements
     * plus one.  On return, the vector will correspond to
     * \f$\{\mathcal{F}^p_1, \cdots, \mathcal{F}^p_{N_e}, \mathcal{F}^p_h\}\f$.
     *
     * @see equilDiffFluxFacsT
     *     for \f$\mathcal{F}_k^T\f$ and \f$\mathcal{F}_h^T\f$
     * @see equilDiffFluxFacsZ
     *     for \f$\mathcal{F}_k^{z_l}\f$ and \f$\mathcal{F}_h^{z_l}\f$
     */
    void equilDiffFluxFacsP(double* const p_F);

    /**
	 * Computes the vector of elemental mass and energy diffusion fluxes per
	 * temperature gradient for a mixture in equilibrium, \f$\mathcal{F}_k^T\f$
	 * and \f$\mathcal{F}_h^T\f$,
	 * where
	 * \f[
	 * \mathcal{F}_k^T = -\sum_{i\in\mathcal{S}} \rho_i \frac{M_{w,k}}{M_{w,i}}
	 * 	   \nu_i^k \sum_{j\in\mathcal{S}} D_{ij} \left(
	 * 	   \frac{\partial x_j}{\partial T} +
	 * 	   \frac{k_{Tj}}{T} \right), \quad \forall\;k\in\mathcal{E},
	 * \f]
	 * and
	 * \f[
	 * \mathcal{F}_h^T = -\sum_{i\in\mathcal{S}} \rho_i h_i
	 * 	   \sum_{j\in\mathcal{S}} D_{ij} \left(
	 * 	   \frac{\partial x_j}{\partial T} +
	 * 	   \frac{k_{Tj}}{T} \right).
	 * \f]
	 *
	 * The total mass diffusion flux of element \f$k\f$ for a mixture at
	 * equilibrium can then be computed with
	 * \f[
	 * \rho_k V_k = \mathcal{F}_k^p \nabla p + \mathcal{F}_k^T \nabla T +
	 *     \sum_{l\in\mathcal{E}} \mathcal{F}_k^{z_l} \nabla z_l,
	 *     \quad \forall\;k\in\mathcal{E}.
	 * \f]
	 * The total energy diffusion flux is computed by
	 * \f[
	 * \sum_{i\in\mathcal{S}} \rho_i h_i V_i = \mathcal{F}_h^p \nabla p +
	 * 	   \mathcal{F}_h^T \nabla T +
	 *     \sum_{l\in\mathcal{E}} \mathcal{F}_h^{z_l} \nabla z_l.
	 * \f]
	 *
	 * @param p_F - A vector with length greater than the number of elements
	 * plus one.  On return, the vector will correspond to
	 * \f$\{\mathcal{F}^T_1, \cdots, \mathcal{F}^T_{N_e}, \mathcal{F}^T_h\}\f$.
	 *
	 * @see equilDiffFluxFacsP
	 *     for \f$\mathcal{F}_k^p\f$ and \f$\mathcal{F}_h^p\f$
	 * @see equilDiffFluxFacsZ
	 *     for \f$\mathcal{F}_k^{z_l}\f$ and \f$\mathcal{F}_h^{z_l}\f$
	 */
    void equilDiffFluxFacsT(double* const p_F);

    /**
	 * Computes the matrix of elemental mass and energy diffusion fluxes per
	 * elemental mole fraction gradients for a mixture in equilibrium,
	 * \f$\mathcal{F}_k^{z_l}\f$ and \f$\mathcal{F}_h^{z_l}\f$,
	 * where
	 * \f[
	 * \mathcal{F}_k^{z_l} = -\sum_{i\in\mathcal{S}} \rho_i
	 * 	   \frac{M_{w,k}}{M_{w,i}}
	 * 	   \nu_i^k \sum_{j\in\mathcal{S}} D_{ij}
	 * 	   \frac{\partial x_j}{\partial z_l}, \quad \forall\;k,l\in\mathcal{E},
	 * \f]
	 * and
	 * \f[
	 * \mathcal{F}_h^{z_l} = -\sum_{i\in\mathcal{S}} \rho_i h_i
	 * 	   \sum_{j\in\mathcal{S}} D_{ij} \frac{\partial x_j}{\partial z_l},
	 * 	   \quad \forall\;l\in\mathcal{E}.
	 * \f]
	 *
	 * The total mass diffusion flux of element \f$k\f$ for a mixture at
	 * equilibrium can then be computed with
	 * \f[
	 * \rho_k V_k = \mathcal{F}_k^p \nabla p + \mathcal{F}_k^T \nabla T +
	 *     \sum_{l\in\mathcal{E}} \mathcal{F}_k^{z_l} \nabla z_l,
	 *     \quad \forall\;k\in\mathcal{E}.
	 * \f]
	 * The total energy diffusion flux is computed by
	 * \f[
	 * \sum_{i\in\mathcal{S}} \rho_i h_i V_i = \mathcal{F}_h^p \nabla p +
	 * 	   \mathcal{F}_h^T \nabla T +
	 *     \sum_{l\in\mathcal{E}} \mathcal{F}_h^{z_l} \nabla z_l.
	 * \f]
	 *
	 * @param p_F - A pointer to a (\f$N_e+1\times N_e\f$)
	 * matrix stored in row-major order.  On return, the matrix will correspond
	 * to
	 * \f[
	 * \begin{bmatrix}
	 * \mathcal{F}_1^{z_1} & \cdots & \mathcal{F}_1^{z_{N_e}} \\
	 * \vdots & \ddots & \vdots \\
	 * \mathcal{F}_{N_e}^{z_1} & \cdots & \mathcal{F}_{N_e}^{z_{N_e}} \\
	 * \mathcal{F}_h^{z_1} & \cdots & \mathcal{F}_h^{z_{N_e}}
	 * \end{bmatrix}.
	 * \f]
	 *
	 * @see equilDiffFluxFacsP
	 *     for \f$\mathcal{F}_k^p\f$ and \f$\mathcal{F}_h^p\f$
	 * @see equilDiffFluxFacsT
	 *     for \f$\mathcal{F}_k^T\f$ and \f$\mathcal{F}_h^T\f$
	 */
    void equilDiffFluxFacsZ(double* const p_F);

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
        
    /// Electric conductivity in S/m (no magnetic field).
    double sigma();
    /// Electric conductivity parallel to the magnetic field in S/m.
    double sigmaParallel();
    /// Electric conductivity perpendicular to the magnetic field in S/m.
    double sigmaPerpendicular();
    /// Electriic conductivity transverse to the magnetic field in S/m.
    double sigmaTransverse();

    /// Mean free path of the mixture in m.
    double meanFreePath();
    /// Mean free path of electrons in m.
    double electronMeanFreePath();

    /// Average heavy particle thermal speed of mixture in m/s.
    double averageHeavyThermalSpeed();
    /// Electron thermal speed of mixture in m/s.
    double electronThermalSpeed();

    /// Electron-heavy collision frequency in 1/s.
    double electronHeavyCollisionFreq();
    /// Average collision frequency of heavy particles in mixture in 1/s.
    double averageHeavyCollisionFreq();

    /// Coulomb mean collision time of the mixture in s.
    double coulombMeanCollisionTime();
    /// Hall parameter.
    double hallParameter();

    // Anisotropic Diffusion Coefficient
    double parallelDiffusionCoefficient();
    double perpDiffusionCoefficient();
    double transverseDiffusionCoefficient();

    // Anisotropic Thermal Diffusion Coefficient
    double parallelThermalDiffusionCoefficient();
    double perpThermalDiffusionCoefficient();
    double transverseThermalDiffusionCoefficient();


    //Anisotropic Electron Thermal Conductivity
    double parallelElectronThermalConductivity();
    double perpElectronThermalConductivity();
    double transverseElectronThermalConductivity();

    // Thermal diffusion ratios
    std::vector<double> parallelThermalDiffusionRatio();
    std::vector<double> perpThermalDiffusionRatio();
    std::vector<double> transverseThermalDiffusionRatio();

    // Return ratios of anisotropic properties
    double ratioSigmaPerpPar();
    double ratioSigmaTransPar();
    double ratioLambdaPerpPar();
    double ratioLambdaTransPar();
    std::vector<double> ratiokTPerpPar();
    std::vector<double> ratiokTTransPar();

private:

    /**
	 * Helper function for computing
	 * \f[
	 * F_k^T = \sum_{i=1}^{N_S}\rho_i\frac{M_{wk}}{M_{wi}}\nu_{ik}
	 *     \sum_{j=1}^{N_S} D_{ij}\alpha_j,
	 * \f]
	 * and
	 * \f[
	 * F^h = \sum_{i=1}^{N_S} h_i \rho_i \sum_{j=1}^{N_S} D_{ij}\alpha_j,
	 * \f]
	 * where \f$\alpha_j\f$ is some species-valued vector.
	 *
	 * @see equilibriumFickP
	 * @see equilibriumFickT
	 * @see equilibriumFickXe
	 */
	void equilDiffFluxFacs(double* const p_F);

private:

    Mutation::Thermodynamics::Thermodynamics& m_thermo;
    CollisionDB* mp_collisions;
    
    ViscosityAlgorithm* mp_viscosity;
    ThermalConductivityAlgorithm* mp_thermal_conductivity;
    Ramshaw* mp_diffusion_matrix;
    
    double* mp_wrk1;
    double* mp_wrk2;
    double* mp_wrk3;
    int* mp_tag;
};

    } // namespace Transport
} // namespace Mutation

#endif // TRANSPORT_TRANSPORT_H
