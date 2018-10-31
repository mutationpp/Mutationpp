/**
 * @file SurfaceBalanceSolverMassEnergy.cpp
 *
 * @brief Class which solves the mass balance and total energy
 * for an interface.
 */

/*
 * Copyright 2014-2018 von Karman Institute for Fluid Dynamics (VKI)
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


#include "Errors.h"
#include "NewtonSolver.h"
#include "Thermodynamics.h"
#include "Transport.h"
#include "Utilities.h"

#include "DiffusionVelocityCalculator.h"
#include "GasFourierHeatFluxCalculator.h"
#include "MassBlowingRate.h"
#include "SolidProperties.h"
#include "Surface.h"
#include "SurfaceChemistry.h"
#include "SurfaceRadiation.h"
#include "SurfaceState.h"

using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace GasSurfaceInteraction {

class SurfaceBalanceSolverMassEnergy :
    public Surface,
    public Mutation::Numerics::NewtonSolver<
        Eigen::VectorXd, SurfaceBalanceSolverMassEnergy>
{
public:
    SurfaceBalanceSolverMassEnergy(ARGS args)
        : m_thermo(args.s_thermo),
          m_surf_state(args.s_surf_state),
          mp_surf_chem(NULL),
          mp_surf_rad(NULL),
          mp_diff_vel_calc(NULL),
          mp_mass_blowing_rate(NULL),
          m_ns(m_thermo.nSpecies()),
          m_nE(m_thermo.nEnergyEqns()),
          m_neqns(m_thermo.nSpecies()+1),
          mv_wdot(m_ns),
          mv_rhoi(m_ns),
          mv_hi(m_ns),
          mv_Vdiff(m_ns),
          mv_X(m_neqns),
          mv_dX(m_neqns),
          mv_f(m_neqns),
          mv_f_unpert(m_neqns),
          m_jac(m_neqns, m_neqns),
          m_pert(1.e-2),
          pos_E(m_ns),
          pos_T_trans(0),
          m_phi(m_surf_state.solidProps().getPhiRatio()),
          m_h_v(m_surf_state.solidProps().getEnthalpyVirginMaterial()),
          set_state_with_rhoi_T(1),
          mv_surf_reac_rates(m_ns),
          is_surf_in_thermal_eq(false),
          is_gas_rad_on(false)
    {
        // Initializing surface chemistry
        mp_surf_chem = new SurfaceChemistry(
            m_thermo,
            args.s_transport,
            args.s_gsi_mechanism,
            args.xml_surf_chem,
            m_surf_state);

        // DiffusionVelocityCalculator
        mp_diff_vel_calc = new DiffusionVelocityCalculator(
            m_thermo, args.s_transport);
        // GasFourierHeatFluxCalculator
        mp_gas_heat_flux_calc = new GasFourierHeatFluxCalculator(
            m_thermo, args.s_transport);

        // Impose thermal equilibrium on the surface
        args.xml_feats.getAttribute(
            "surface_in_thermal_equil", is_surf_in_thermal_eq, false);

        // MassBlowingRate
        DataMassBlowingRate data_mass_blowing_rate = {m_thermo, *mp_surf_chem};
        const std::string s_mass_blowing = "isOn";
        mp_mass_blowing_rate = Factory<MassBlowingRate>::create(
            s_mass_blowing, data_mass_blowing_rate);

        // Surface Radiation
        if (args.xml_surf_rad.tag() == "surface_radiation"){
            args.xml_feats.getAttribute(
                "gas_radiation", is_gas_rad_on, false);
            mp_surf_rad = new SurfaceRadiation(
                m_thermo, args.xml_surf_rad, m_surf_state, is_gas_rad_on);
        }

        // Setup NewtonSolver
        setMaxIterations(5);
        setWriteConvergenceHistory(false);
        setEpsilon(1.e-18);
    }

//=============================================================================

    ~SurfaceBalanceSolverMassEnergy()
    {
        if (mp_surf_chem != NULL) { delete mp_surf_chem; }
        if (mp_surf_rad != NULL) { delete mp_surf_rad; }
        if (mp_diff_vel_calc != NULL) { delete mp_diff_vel_calc; }
        if (mp_gas_heat_flux_calc != NULL) { delete mp_gas_heat_flux_calc; }
        if (mp_mass_blowing_rate != NULL) { delete mp_mass_blowing_rate; }
    }

//=============================================================================

    void computeSurfaceReactionRates(Eigen::VectorXd& v_surf_reac_rates)
    {
        errorSurfaceStateNotSet();

        v_surf_reac_rates.setZero();
        if (mp_surf_chem != NULL)
            mp_surf_chem->surfaceReactionRates(v_surf_reac_rates);
    }

//=============================================================================

    Eigen::VectorXd computeSurfaceReactionRatesPerReaction()
    {
        const int nr = nSurfaceReactions();
        Eigen::VectorXd v_wrk(nr);

        if (mp_surf_chem != NULL && nr > 0){
            mp_surf_chem->surfaceReactionRatesPerReaction(v_wrk);
        }
        return v_wrk;
    }

//=============================================================================

    int nSurfaceReactions()
    {
        if (mp_surf_chem != NULL)
            return mp_surf_chem->nSurfaceReactions();

        return 0;
    }

//=============================================================================

    void setDiffusionModel(
        const Eigen::VectorXd& v_mole_frac_edge, const double& dx)
    {
        mp_diff_vel_calc->setDiffusionModel(v_mole_frac_edge, dx);
    }

//=============================================================================

    void setGasFourierHeatFluxModel(
        const Eigen::VectorXd& v_T_edge, const double& dx){
        mp_gas_heat_flux_calc->setGasFourierHeatFluxModel(v_T_edge, dx);
    }

//=============================================================================

    virtual void setGasRadHeatFlux(const double& gas_rad_heat_flux)
    {
        if (mp_surf_rad != NULL)
            mp_surf_rad->gasRadiativeHeatFlux(gas_rad_heat_flux);
    }

//=============================================================================

    void solveSurfaceBalance()
    {
        // errorUninitializedDiffusionModel
        errorSurfaceStateNotSet();

    	// Getting the state
        mv_rhoi = m_surf_state.getSurfaceRhoi();
        double T_trans = m_surf_state.getSurfaceT()(pos_T_trans);

        // Impose equilibrium
        if (is_surf_in_thermal_eq){
             mv_X.tail(m_nE).setConstant(T_trans);
        } else {
            mv_X.tail(m_nE) = m_surf_state.getSurfaceT();
        }

        saveUnperturbedPressure(mv_rhoi, mv_X.tail(m_nE));

        // Changing to the solution variables and solving
        computeMoleFracfromPartialDens(mv_rhoi, mv_X.tail(m_nE), mv_X);
        mv_X = solve(mv_X);

        computePartialDensfromMoleFrac(
            mv_X.head(m_ns), mv_X.tail(m_nE), mv_rhoi);

        // Setting the state again
        m_surf_state.setSurfaceState(
            mv_rhoi.data(), mv_X.tail(m_nE).data(), set_state_with_rhoi_T);
    }

//==============================================================================

    void setIterationsSurfaceBalance(const int& iter){ setMaxIterations(iter); }

//==============================================================================

    double massBlowingRate()
    {
        if (mp_surf_chem != NULL)
            return mp_mass_blowing_rate->computeBlowingFlux();
        return 0.;
    }

//==============================================================================

    void updateFunction(Eigen::VectorXd& v_X)
    {
        mv_f.setZero();

    	// Setting Initial Gas and Surface State;
        computePartialDensfromMoleFrac(
            v_X.head(m_ns), v_X.tail(m_nE), mv_rhoi);

        m_thermo.setState(
            mv_rhoi.data(), v_X.tail(m_nE).data(), set_state_with_rhoi_T);
        m_surf_state.setSurfaceState(
            mv_rhoi.data(), v_X.tail(m_nE).data(), set_state_with_rhoi_T);

        // Diffusion Fluxes
        mp_diff_vel_calc->computeDiffusionVelocities(
            v_X.head(m_ns), mv_Vdiff);
        mv_f = mv_rhoi.cwiseProduct(mv_Vdiff);

        // Chemical Production Rates
        computeSurfaceReactionRates(mv_surf_reac_rates);
        mv_f -= mv_surf_reac_rates;

        // Blowing flux
        double mass_blow = mp_mass_blowing_rate->computeBlowingFlux(
            mv_surf_reac_rates);
        mv_f.head(m_ns) += mv_rhoi * mass_blow / mv_rhoi.sum();

        // Energy
        m_thermo.getEnthalpiesMass(mv_hi.data());
        double hmix = m_thermo.mixtureHMass();

        mv_f(pos_E) += mv_hi.dot(mv_Vdiff.cwiseProduct(mv_rhoi));
        mv_f(pos_E) += mp_gas_heat_flux_calc->
                           computeGasFourierHeatFlux(v_X.tail(m_nE));
        mv_f(pos_E) += hmix * mass_blow;

        // Radiation
        if (mp_surf_rad != NULL)
            mv_f(pos_E) += mp_surf_rad->surfaceNetRadiativeHeatFlux();

        // Steady state assumption virgin material enthalpy
        mv_f(pos_E) -= mass_blow * (1+m_phi) * m_h_v;
    }

//==============================================================================

    void updateJacobian(Eigen::VectorXd& v_X)
    {
        // Perturbing Mass
        mv_f_unpert = mv_f;
        for (int i_ns = 0; i_ns < m_ns; i_ns++){
            double X_unpert = v_X(i_ns);
            double pert = m_pert;
            v_X(i_ns) += pert;

            updateFunction(v_X);

            // Update Jacobian column
            m_jac.col(i_ns) = (mv_f-mv_f_unpert) / pert;

            // Unperturb mole fractions
            v_X(i_ns) = X_unpert;
        }

        // Perturbing Energy
        double T_pert = m_pert;
        double X_unpert = v_X(pos_E);
        v_X(pos_E) += T_pert;

        updateFunction(v_X);
        m_jac.col(pos_E) = (mv_f-mv_f_unpert) / T_pert;

        v_X(pos_E) = X_unpert;
    }

//==============================================================================

    Eigen::VectorXd& systemSolution()
    {
        double a = m_jac.diagonal().maxCoeff();

        mv_dX = (m_jac + a*Eigen::MatrixXd::Ones(m_ns,m_ns)).
                    fullPivLu().solve((1 + a)*mv_f_unpert);
        return mv_dX;
    }
//==============================================================================

    double norm() {
        // return mv_dX.lpNorm<Eigen::Infinity>();
        return mv_f.lpNorm<Eigen::Infinity>();
    }

//==============================================================================
private:
    void saveUnperturbedPressure(
        const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_T)
    {
        m_thermo.setState(
            v_rhoi.data(), v_T.data(), set_state_with_rhoi_T);
        m_Psurf = m_thermo.P();
    }
//==============================================================================

    void computeMoleFracfromPartialDens(
        const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_T,
        Eigen::VectorXd& v_xi)
    {
        m_thermo.setState(
            v_rhoi.data(), v_T.data(), set_state_with_rhoi_T);
        v_xi.head(m_ns) = Eigen::Map<const Eigen::VectorXd>(
            m_thermo.X(), m_ns);
    }
//==============================================================================

    void computePartialDensfromMoleFrac(
        const Eigen::VectorXd& v_xi, const Eigen::VectorXd& v_T,
        Eigen::VectorXd& v_rhoi)
    {
    	v_rhoi = v_xi.cwiseProduct(m_thermo.speciesMw().matrix()) *
    			  m_Psurf / (v_T(pos_T_trans) * RU);
    }

//==============================================================================

    void errorSurfaceStateNotSet() const
    {
        if (!m_surf_state.isSurfaceStateSet()) {
            throw LogicError()
            << "The surface state must have been set!";
        }
    }

//==============================================================================
private:
    Mutation::Thermodynamics::Thermodynamics& m_thermo;

    SurfaceState& m_surf_state;

    SurfaceChemistry* mp_surf_chem;
    SurfaceRadiation* mp_surf_rad;
    DiffusionVelocityCalculator* mp_diff_vel_calc;
    GasFourierHeatFluxCalculator* mp_gas_heat_flux_calc;

    MassBlowingRate* mp_mass_blowing_rate;

    bool is_surf_in_thermal_eq;
    bool is_gas_rad_on;

    const size_t m_ns;
    const size_t m_nE;
    const size_t m_neqns;

    Eigen::VectorXd mv_Tsurf;
    double m_Psurf;

    Eigen::VectorXd mv_wdot;
    Eigen::VectorXd mv_Vdiff;
    Eigen::VectorXd mv_hi;

    Eigen::VectorXd mv_rhoi;
    Eigen::VectorXd mv_X;
    Eigen::VectorXd mv_dX;
    Eigen::VectorXd mv_f;
    Eigen::MatrixXd m_jac;
    double m_pert;
    Eigen::VectorXd mv_f_unpert;
    Eigen::VectorXd mv_surf_reac_rates;

    const double m_phi;
    const double m_h_v;

    const int pos_E;
    const int pos_T_trans;
    const int set_state_with_rhoi_T;
};

ObjectProvider<
    SurfaceBalanceSolverMassEnergy, Surface>
    surface_balance_solver_mass_energy_gamma("gamma_energy");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
