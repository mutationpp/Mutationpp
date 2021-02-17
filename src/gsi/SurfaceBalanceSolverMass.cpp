/**
 * @file SurfaceBalanceSolverMass.cpp
 *
 * @brief Class which solves the mass balance for an interface.
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


#include "Errors.h"
#include "NewtonSolver.h"
#include "Thermodynamics.h"
#include "Transport.h"
#include "Utilities.h"

#include "DiffusionVelocityCalculator.h"
#include "MassBlowingRate.h"
#include "Surface.h"
#include "SurfaceChemistry.h"
#include "SurfaceState.h"

using namespace Mutation::Utilities::Config;
using namespace Mutation::Utilities::IO;

namespace Mutation {
    namespace GasSurfaceInteraction {

class SurfaceBalanceSolverMass :
    public Surface,
    public Mutation::Numerics::NewtonSolver<
        Eigen::VectorXd, SurfaceBalanceSolverMass>
{
public:
    SurfaceBalanceSolverMass(ARGS args)
        : m_thermo(args.s_thermo),
          m_surf_state(args.s_surf_state),
          mp_surf_chem(NULL),
          mp_diff_vel_calc(NULL),
          mp_mass_blowing_rate(NULL),
          m_ns(m_thermo.nSpecies()),
          m_nE(m_thermo.nEnergyEqns()),
          mv_wdot(m_ns),
          mv_rhoi(m_ns),
          mv_Tsurf(m_nE),
          mv_X(m_ns),
          mv_dX(m_ns),
          mv_f(m_ns),
          mv_f_unpert(m_ns),
          m_jac(m_ns, m_ns),
          m_tol(1.e-13),
          m_pert(1.e-2),
          mv_X_unpert(m_ns),
          pos_T_trans(0),
          set_state_with_rhoi_T(1),
          mv_surf_reac_rates(m_ns)
    {
        // Initializing surface chemistry
        if (args.xml_surf_chem.tag() == "surface_chemistry"){
            mp_surf_chem = new SurfaceChemistry(
                m_thermo,
                args.s_transport,
                args.s_gsi_mechanism,
                args.xml_surf_chem,
                m_surf_state);
        }

        // DiffusionVelocityCalculator
        mp_diff_vel_calc = new DiffusionVelocityCalculator(
            m_thermo, args.s_transport);

        // MassBlowingRate
        DataMassBlowingRate data_mass_blowing_rate = {m_thermo, *mp_surf_chem};
        const std::string s_mass_blowing = "isOn";
        mp_mass_blowing_rate = Factory<MassBlowingRate>::create(
            s_mass_blowing, data_mass_blowing_rate);

        // Setup NewtonSolver
        setMaxIterations(5);
        setWriteConvergenceHistory(false);
        setEpsilon(m_tol);
    }

//=============================================================================

    ~SurfaceBalanceSolverMass()
    {
        if (mp_surf_chem != NULL) { delete mp_surf_chem; }
        if (mp_diff_vel_calc != NULL) { delete mp_diff_vel_calc; }
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

        if (mp_surf_chem != NULL){
            mp_surf_chem->surfaceReactionRatesPerReaction(v_wrk);
            return v_wrk;
        }
        throw LogicError()
            << "computeGSIReactionRatePerReaction cannot be invoked "
            << "without defining a surface_chemistry option in "
            << "Gas-Surface Interaction input file.";

        return v_wrk.setZero();
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

    void solveSurfaceBalance()
    {
        // errorUninitializedDiffusionModel
        errorSurfaceStateNotSet();

    	// Getting the state
        mv_rhoi = m_surf_state.getSurfaceRhoi();
        mv_Tsurf = m_surf_state.getSurfaceT();

        saveUnperturbedPressure(mv_rhoi);

        // Changing to the solution variables
        computeMoleFracfromPartialDens(mv_rhoi, mv_X);

        // Solving
        mv_X = solve(mv_X);

        applyTolerance(mv_X);
        computePartialDensfromMoleFrac(mv_X, mv_rhoi);

        // Setting the state again
        m_surf_state.setSurfaceState(
            mv_rhoi.data(), mv_Tsurf.data(), set_state_with_rhoi_T);
    }

//==============================================================================

    void setIterationsSurfaceBalance(const int& iter){ setMaxIterations(iter); }

//==============================================================================

    double massBlowingRate() {
        return mp_mass_blowing_rate->computeBlowingFlux();
    }

//==============================================================================

    void updateFunction(Eigen::VectorXd& v_mole_frac)
    {
        applyTolerance(v_mole_frac);

    	// Setting Initial Gas and Surface State;
        computePartialDensfromMoleFrac(v_mole_frac, mv_rhoi);
        m_thermo.setState(
            mv_rhoi.data(), mv_Tsurf.data(), set_state_with_rhoi_T);

        m_surf_state.setSurfaceState(
            mv_rhoi.data(), mv_Tsurf.data(), set_state_with_rhoi_T);

        // Diffusion Fluxes
        mp_diff_vel_calc->computeDiffusionVelocities(v_mole_frac, mv_f);
        applyTolerance(mv_f);
        mv_f = mv_rhoi.cwiseProduct(mv_f);

        // Chemical Production Rates
        computeSurfaceReactionRates(mv_surf_reac_rates);
        mv_f -= mv_surf_reac_rates;

        // Blowing Fluxes
        double mass_blow = mp_mass_blowing_rate->computeBlowingFlux(
            mv_surf_reac_rates);
        mv_f += mv_rhoi * mass_blow / mv_rhoi.sum();
    }

//==============================================================================

    void updateJacobian(Eigen::VectorXd& v_mole_frac)
    {
        mv_f_unpert = mv_f;
        for (int i_ns = 0 ; i_ns < m_ns ; i_ns++){
            mv_X_unpert = v_mole_frac;
            double pert = m_pert;
            v_mole_frac(i_ns) += pert;

            updateFunction(v_mole_frac);

            // Update Jacobian column
            m_jac.col(i_ns) = (mv_f - mv_f_unpert) / pert;

            // Unperturbed mole fractions
            v_mole_frac = mv_X_unpert;
        }
    }

//==============================================================================

    Eigen::VectorXd& systemSolution()
    {
        double a = (m_jac.diagonal()).cwiseAbs().maxCoeff();
        mv_dX = (m_jac + a*Eigen::MatrixXd::Ones(m_ns,m_ns)).
            fullPivLu().solve(mv_f_unpert);

        applyTolerance(mv_dX);
        return mv_dX;
    }
//==============================================================================

    double norm()
    {
        // return mv_f_unpert.lpNorm<Eigen::Infinity>();
        return mv_dX.lpNorm<Eigen::Infinity>();
    }

//==============================================================================
private:
    void saveUnperturbedPressure(Eigen::VectorXd& v_rhoi) {
        m_thermo.setState(
            v_rhoi.data(), mv_Tsurf.data(), set_state_with_rhoi_T);
        m_Psurf = m_thermo.P();
    }
//==============================================================================

    void computeMoleFracfromPartialDens(
        Eigen::VectorXd& v_rhoi, Eigen::VectorXd& v_xi) {
        m_thermo.setState(
            v_rhoi.data(), mv_Tsurf.data(), set_state_with_rhoi_T);
        v_xi = Eigen::Map<const Eigen::VectorXd>(m_thermo.X(), m_ns);
    }
//==============================================================================

    void computePartialDensfromMoleFrac(
        Eigen::VectorXd& v_xi, Eigen::VectorXd& v_rhoi) {
    	v_rhoi = v_xi.cwiseProduct(m_thermo.speciesMw().matrix()) *
    			  m_Psurf / (mv_Tsurf(pos_T_trans) * RU);
    }
//==============================================================================

    void errorSurfaceStateNotSet() const {
        if (!m_surf_state.isSurfaceStateSet()) {
            throw LogicError()
                << "The surface state must have been set!";
        }
    }
//==============================================================================

    inline void applyTolerance(Eigen::VectorXd& v_x) const {
        for (int i = 0; i < m_ns; i++)
            if (std::abs(v_x(i)) < m_tol) v_x(i) = 0.;
    }
//==============================================================================
private:
    Mutation::Thermodynamics::Thermodynamics& m_thermo;

    SurfaceChemistry* mp_surf_chem;
    DiffusionVelocityCalculator* mp_diff_vel_calc;
    MassBlowingRate* mp_mass_blowing_rate;
    SurfaceState& m_surf_state;

    const size_t m_ns;
    const size_t m_nE;

    Eigen::VectorXd mv_Tsurf;
    double m_Psurf;

    Eigen::VectorXd mv_wdot;

    Eigen::VectorXd mv_rhoi;
    Eigen::VectorXd mv_X;
    Eigen::VectorXd mv_dX;
    Eigen::VectorXd mv_f;
    Eigen::MatrixXd m_jac;
    const double m_tol;
    double m_pert;
    Eigen::VectorXd mv_X_unpert;
    Eigen::VectorXd mv_f_unpert;
    Eigen::VectorXd mv_surf_reac_rates;


    const int pos_T_trans;
    const int set_state_with_rhoi_T;
};

ObjectProvider<
    SurfaceBalanceSolverMass, Surface>
    surface_balance_solver_phenomenological_mass("phenomenological_mass");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
