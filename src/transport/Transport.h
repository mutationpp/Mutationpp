/**
 * @file Transport.h
 *
 * @brief Declaration of Transport class.
 */

/*
 * Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
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
#include "ElectronSubSystem.h"
#include "Utilities.h"

#include <Eigen/Dense>

namespace Mutation {
    namespace Transport {

class DiffusionMatrix;
class ThermalConductivityAlgorithm;
class ViscosityAlgorithm;

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
    
    /// Provides reference to the underlying collision integral database.
    CollisionDB& collisionDB() { return m_collisions; }
    
    /// Sets the viscosity algorithm.
    void setViscosityAlgo(const std::string& algo);
    
    /// Sets the heavy particle, thermal conductivity algorithm.
    void setThermalConductivityAlgo(const std::string& algo);
    
    /// Sets the diffusion matrix algorithm.
    void setDiffusionMatrixAlgo(const std::string& algo);

    /// Returns the number of collision pairs accounted for in this mixture.
    int nCollisionPairs() const { return m_collisions.size(); }
    
    //void omega11ii(double* const p_omega);
    //void omega22ii(double* const p_omega);
    
    /// Returns the mixture viscosity.
    double viscosity();
    
    /**
     * Returns the mixture thermal conductivity for a frozen mixture.
     * To be used only at thermal equilibrium.
     */
    double frozenThermalConductivity() {
        return
            heavyThermalConductivity() + 
            electronThermalConductivity() +
            internalThermalConductivity(m_thermo.T());
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
        return
            frozenThermalConductivity() +
            reactiveThermalConductivity() +
            soretThermalConductivity();
    }
    
    /**
     * Returns the heavy particle translational thermal conductivity using the 
     * set algorithm.
     */
    double heavyThermalConductivity();
    
    /**
     * Returns the thermal conductivity of an internal energy mode using
     * Euken's formula.
     */
    template <typename E>
    double euken(const Eigen::ArrayBase<E>& cp)
    {
        const int ns = m_thermo.nGas();  assert(cp.size() == ns);
        const int nh = m_thermo.nHeavy();
        const int k  = ns-nh;

        Eigen::Map<const Eigen::ArrayXd> X(m_thermo.X()+k, nh);
        static Eigen::ArrayXd avDij; avDij.resize(nh);
        const Eigen::ArrayXd& nDij = m_collisions.nDij();

        avDij.setZero();
        for (int i = 0, index = 0; i < nh; ++i) {
            // Account for diagonal term added twice
            avDij(i) -= X(i)/nDij(index);
            for (int j = i; j < nh; ++j, ++index) {
                avDij(i) += X(j)/nDij(index);
                avDij(j) += X(i)/nDij(index);
            }
        }

        return KB * (X * cp.tail(nh) / avDij).sum();
    }

    /**
     * Returns the internal energy thermal conductivity using Euken's formulas.
     */
    double internalThermalConductivity(double T);

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
     * Returns the reactive thermal conductivity which accounts for reactions
     * for mixtures in thermochemical equilibrium using the Butler-Brokaw
     * formula.
     */
    double butlerBrokawThermalConductivity();

    /**
     * Returns the Soret thermal conductivity the mixture in thermochemical 
     * equilibrium.
     */
    double soretThermalConductivity();
    
    /// Returns the heavy thermal diffusion ratios for each species.
    void heavyThermalDiffusionRatios(double* const p_k);
    
    /// Returns the multicomponent diffusion coefficient matrix.
    const Eigen::MatrixXd& diffusionMatrix();

    /**
     * Returns the average diffusion coefficients.
     * \f[ D_{im} = \frac{(1-x_i)}{\sum_{j\ne i}x_j/\mathscr{D}_{ij}} \f]
     */
    void averageDiffusionCoeffs(double *const p_Di) {
        Eigen::Map<Eigen::ArrayXd>(p_Di, m_thermo.nGas()) =
            m_collisions.Dim();
    }
    
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
     * @todo Check that this is correct for thermal nonequilibrium.
     *
     * @param p_dp - the vector of modified driving forces \f$ d^{'}_i \f$
     * @param p_V  - on return, the vector of diffusion velocities
     * @param E    - on return, the ambipolar electric field
     */
    void stefanMaxwell(
        const double* const p_dp, double* const p_V, double& E, int order = 1);

    /**
     * Computes the species diffusion velocities and ambipolar electric field 
     * using the Ramshaw approximation of the generalized Stefan-Maxwell 
     * equations and the supplied modified driving forces
     *
     * \f[ d^{'}_i = \frac{\nabla p_i}{n k_B T_h} - \frac{y_i p}{n k_B T_h}
     *     \nabla \ln p + k^h_{Ti} \nabla \ln T_h + k^e_{Ti} \frac{T_h}{T_e}
     *     \nabla \ln T_e \f]
     *
     * @param Th   - heavy particle translational temperature
     * @param Te   - free electron translational temperature
     * @param p_dp - the vector of modified driving forces \f$ d^{'}_i \f$
     * @param p_V  - on return, the vector of diffusion velocities
     * @param E    - on return, the ambipolar electric field
     */
    void stefanMaxwell(double Th, double Te,
        const double* const p_dp, double* const p_V, double& E, int order = 1);

    /// Computes the electron corrections for the Stefan-Maxwell equations.
    void smCorrectionsElectron(int order, Eigen::ArrayXd& phi);
    
    /// Computes the heavy corrections for the Stefan-Maxwell equations.
    void smCorrectionsHeavy(int order, Eigen::ArrayXd& phi);
        
    //==========================================================================
    // Electron subsystem functions
    //==========================================================================

    /// Isotropic electric conductivity in S/m (no magnetic field).
    double electricConductivity(int order = 3) {
        return mp_esubsyst->electricConductivity(order);
    }

    /// Anisotropic electric conductivity in S/m (with magnetic field).
    Eigen::Vector3d electricConductivityB(int order = 3) {
        return mp_esubsyst->electricConductivityB(order);
    }

    /// Returns the electron thermal conductivity in W/m-K
    double electronThermalConductivity(int order = 3) {
        return mp_esubsyst->electronThermalConductivity(order);
    }

    /// Anisotropic electron thermal conductivity in W/m-K.
    Eigen::Vector3d electronThermalConductivityB(int order = 3) {
        return mp_esubsyst->electronThermalConductivityB(order);
    }

    /// Isotropic electron diffusion coefficient.
    double electronDiffusionCoefficient(int order = 3) {
        return mp_esubsyst->electronDiffusionCoefficient(order);
    }

    /// Anisotropic electron diffusion coefficient.
    Eigen::Vector3d electronDiffusionCoefficientB(int order = 3) {
        return mp_esubsyst->electronDiffusionCoefficientB(order);
    }

    /// Isotropic second-order electron diffusion coefficient.
    double electronDiffusionCoefficient2(int order = 3)
    {
        const int nh = m_thermo.nHeavy();
        return mp_esubsyst->electronDiffusionCoefficient2(
            diffusionMatrix().bottomRightCorner(nh,nh), order);
    }

    /// Anisotropic second-order electron diffusion coefficient.
    Eigen::Vector3d electronDiffusionCoefficient2B(int order = 3)
    {
        const int nh = m_thermo.nHeavy();
        return mp_esubsyst->electronDiffusionCoefficient2B(
            diffusionMatrix().bottomRightCorner(nh,nh), order);
    }

    /// Isotropic alpha coefficients.
    const Eigen::VectorXd& alpha(int order = 3) {
        return mp_esubsyst->alpha(order);
    }

    /// Anisotropic alpha coefficients.
    const Eigen::Matrix<double,-1,3>& alphaB(int order = 3) {
        return mp_esubsyst->alphaB(order);
    }

    /// Isotropic electron thermal diffusion ratio.
    double electronThermalDiffusionRatio(int order = 3) {
        return mp_esubsyst->electronThermalDiffusionRatio(order);
    }

    /// Anisotropic electron thermal diffusion ratio.
    Eigen::Vector3d electronThermalDiffusionRatioB(int order = 3) {
        return mp_esubsyst->electronThermalDiffusionRatioB(order);
    }

    /// Isotropic second-order electron thermal diffusion ratios.
    const Eigen::VectorXd& electronThermalDiffusionRatios2(int order = 3) {
        return mp_esubsyst->electronThermalDiffusionRatios2(order);
    }

    /// Anisotropic second-order electron thermal diffusion ratios.
    const Eigen::Matrix<double,-1,3>& electronThermalDiffusionRatios2B(int order = 3) {
        return mp_esubsyst->electronThermalDiffusionRatios2B(order);
    }


    /// Mean free path of the mixture in m.
    double meanFreePath();
    /// Mean free path of electrons in m.
    double electronMeanFreePath();
    /// Thermal speed of species i in m/s
    double speciesThermalSpeed(const int& i) const;
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

    // @todo Determine if these functions are necessary
    
    // // Anisotropic Diffusion Coefficient
    // double parallelDiffusionCoefficient();
    // double perpDiffusionCoefficient();
    // double transverseDiffusionCoefficient();

    // // Anisotropic Thermal Diffusion Coefficient
    // double parallelThermalDiffusionCoefficient();
    // double perpThermalDiffusionCoefficient();
    // double transverseThermalDiffusionCoefficient();

    // //Anisotropic Electron Thermal Conductivity
    // double parallelElectronThermalConductivity();
    // double perpElectronThermalConductivity();
    // double transverseElectronThermalConductivity();

    // // Thermal diffusion ratios
    // std::vector<double> parallelThermalDiffusionRatio2();
    // double parallelThermalDiffusionRatio();
    // std::vector<double> perpThermalDiffusionRatio2();
    // double perpThermalDiffusionRatio();
    // std::vector<double> transverseThermalDiffusionRatio2();
    // double transverseThermalDiffusionRatio();

    // double taueLambda();
    // double taueLambdaBr();
    // double perpLambdaeBr();
    // double transLambdaeBr();

    // double coefficientFriction();
    // double perpfriccoeffBr();
    // double transfriccoeffBr();

    // double tauViscosity();
    // double tauViscosityBr();
    // double tauLambdaHeavy();
    // double tauLambdaHeavyBr();
    // double tauEnergy();
    // double tauEnergyBr();
    // double etaohmbraginskii();
    // double transetaohmBr();
    // double perpetaohmBr();


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
     * @todo Document how this modifies the mp_wrk1 and mp_wrk2 arrays.
	 *
	 * @see equilibriumFickP
	 * @see equilibriumFickT
	 * @see equilibriumFickXe
	 */
	void equilDiffFluxFacs(double* const p_F);

private:

    Mutation::Thermodynamics::Thermodynamics& m_thermo;
    CollisionDB m_collisions;
    ElectronSubSystem* mp_esubsyst;
    
    ViscosityAlgorithm* mp_viscosity;
    ThermalConductivityAlgorithm* mp_thermal_conductivity;
    DiffusionMatrix* mp_diffusion_matrix;
    
    double* mp_wrk1;
    double* mp_wrk2;
    double* mp_wrk3;
    int* mp_tag;
};

    } // namespace Transport
} // namespace Mutation

#endif // TRANSPORT_TRANSPORT_H
