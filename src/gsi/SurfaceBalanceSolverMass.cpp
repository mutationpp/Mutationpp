#include "Errors.h"
#include "NewtonSolver.h"
#include "Thermodynamics.h"
#include "Transport.h"
#include "Utilities.h"

#include "DiffusionVelocityCalculator.h"
#include "MassBlowingRate.h"
#include "WallProductionTerms.h"
#include "SurfaceBalanceSolver.h"
#include "SurfaceProperties.h"
#include "WallState.h"

using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace GasSurfaceInteraction {

class SurfaceBalanceSolverMass :
    public SurfaceBalanceSolver,
    public Mutation::Numerics::NewtonSolver<
        Eigen::VectorXd, SurfaceBalanceSolverMass>
{
public:
SurfaceBalanceSolverMass(ARGS args)
    : m_thermo(args.s_thermo),
      m_surf_props(args.s_surf_props),
      m_wall_state(args.s_wall_state),
      mp_diff_vel_calc(NULL),
      mp_mass_blowing_rate(NULL),
      m_ns(m_thermo.nSpecies()),
      m_nE(m_thermo.nEnergyEqns()),
      mv_rhoi(m_ns),
      mv_Twall(m_nE),
      mv_X(m_ns),
      mv_dX(m_ns),
      mv_f(m_ns),
      mv_f_unpert(m_ns),
      mv_jac(m_ns, m_ns),
      m_pert(1.e-7),
      pos_T_trans(0),
      set_state_with_rhoi_T(1),
      mv_sep_mass_prod_rate(m_ns)
{
	Mutation::Utilities::IO::XmlElement::const_iterator iter_prod_terms =
                                   args.s_node_prod_terms.begin();

    std::string s_tag;
    for(; iter_prod_terms != args.s_node_prod_terms.end(); ++iter_prod_terms)
    {
        DataWallProductionTerms data_wall_prod_terms = { m_thermo,
                                                         args.s_transport,
                                                         args.s_gsi_mechanism,
                                                         *iter_prod_terms,
                                                         m_surf_props,
                                                         m_wall_state,
                                                         &mv_surf_prod,
                                                         &m_Pwall};

        s_tag = iter_prod_terms->tag();
        addSurfaceProductionTerm(Factory<WallProductionTerms>::create(
            s_tag, data_wall_prod_terms));
    }
    errorEmptyWallProductionTerms();
 
    // DiffusionVelocityCalculator
    mp_diff_vel_calc = new DiffusionVelocityCalculator(
        m_thermo, args.s_transport);

    // MassBlowingRate
    DataMassBlowingRate data_mass_blowing_rate = {m_thermo, mv_surf_prod};
    const std::string s_mass_blowing = "isOn";
    mp_mass_blowing_rate = Factory<MassBlowingRate>::create(
    s_mass_blowing, data_mass_blowing_rate);

    // Setup NewtonSolver
    setMaxIterations(5);
    setWriteConvergenceHistory(false);
    setEpsilon(1.e-18);
}

//=============================================================================

    ~SurfaceBalanceSolverMass()
    {
        for (std::vector<WallProductionTerms*>::iterator iter =
                 mv_surf_prod.begin();
             iter != mv_surf_prod.end();
             ++iter)
        {
            delete (*iter);
        }
        mv_surf_prod.clear();

        if (mp_diff_vel_calc != NULL) { delete mp_diff_vel_calc; }
        if (mp_mass_blowing_rate != NULL) { delete mp_mass_blowing_rate; }
    }

//=============================================================================

    Eigen::VectorXd computeGSIProductionRates()
    {
        errorWallStateNotSet();
        static Eigen::VectorXd v_wrk(m_ns);

        mv_sep_mass_prod_rate.setZero();
        for (int i_prod_terms = 0;
             i_prod_terms < mv_surf_prod.size();
             ++i_prod_terms) {
            v_wrk.setZero();
            mv_surf_prod[i_prod_terms]->productionRate(v_wrk);
            mv_sep_mass_prod_rate += v_wrk;
        }
        return mv_sep_mass_prod_rate;
    }

//=============================================================================

    void setDiffusionModel(
        const Eigen::VectorXd& v_mole_frac_edge,
        const double& dx)
    {
        mp_diff_vel_calc->setDiffusionModel(v_mole_frac_edge, dx);
    }

//=============================================================================

    void solveSurfaceBalance()
    {
        // errorUninitializedDiffusionModel
        // errorWallStateNotSetYet

    	// Getting the state
        mv_rhoi = m_wall_state.getWallRhoi();

        for (int i_E = 0; i_E < m_nE; i_E++) {
            mv_Twall(i_E) = m_wall_state.getWallT()(i_E);
        }
        saveUnperturbedPressure(mv_rhoi);

        // Changing to the solution variables and solving
        computeMoleFracfromPartialDens(mv_rhoi, mv_X);

        mv_X = solve(mv_X);

        computePartialDensfromMoleFrac(mv_X, mv_rhoi);

        // Setting the state again
        m_wall_state.setWallState(
            mv_rhoi.data(),
            mv_Twall.data(),
            set_state_with_rhoi_T);
    }

//==============================================================================

    double massBlowingRate()
    {
        double mdot;
        mdot = mp_mass_blowing_rate->computeBlowingFlux();

        return mdot;
    }

//==============================================================================

    void updateFunction(Eigen::VectorXd& v_mole_frac)
    {
    	// Setting Initial Gas and Wall State;
        computePartialDensfromMoleFrac(v_mole_frac, mv_rhoi);
        m_thermo.setState(mv_rhoi.data(),
                          mv_Twall.data(),
                          set_state_with_rhoi_T);
        m_wall_state.setWallState(mv_rhoi.data(),
                                  mv_Twall.data(),
                                  set_state_with_rhoi_T);

        // Diffusion Fluxes
        mp_diff_vel_calc->computeDiffusionVelocities(v_mole_frac, mv_f);
        mv_f = mv_rhoi.cwiseProduct(mv_f);

        // Chemical Production Rates
        mv_f -= computeGSIProductionRates();

        // Blowing Fluxes
        mv_f += mv_rhoi*mp_mass_blowing_rate->computeBlowingFlux()
        		/mv_rhoi.sum();
    }
    
//==============================================================================

    void updateJacobian(Eigen::VectorXd& v_mole_frac)
    {
        mv_f_unpert = mv_f;
        for ( int i_ns = 0 ; i_ns < m_ns ; i_ns++){
            m_X_unpert = v_mole_frac(i_ns);
            v_mole_frac(i_ns) += m_pert;

            updateFunction(v_mole_frac);

            // Update Jacobian column
            //@todo Can be improved using .col() function Eigen library
            for( int i_eq = 0 ; i_eq < m_ns ; i_eq++ ){
                mv_jac(i_eq, i_ns) = (mv_f(i_eq) - mv_f_unpert(i_eq))/m_pert;
            }

            // Unperturbed mole fractions
            v_mole_frac(i_ns) = m_X_unpert;
        }
    }

//==============================================================================

    Eigen::VectorXd& systemSolution()
    {
        mv_dX = mv_jac.fullPivLu().solve(mv_f_unpert);
        return mv_dX;
    }

//==============================================================================

    double norm()
    {
        // return v_f.lpNorm<Eigen::Infinity>();
        return mv_dX.lpNorm<Eigen::Infinity>();
    }

private:
    void addSurfaceProductionTerm(WallProductionTerms* p_wall_prod_term){
        mv_surf_prod.push_back(p_wall_prod_term);
    }

//==============================================================================

    void saveUnperturbedPressure(Eigen::VectorXd& v_rhoi)
    {
        m_thermo.setState(
            v_rhoi.data(),
            mv_Twall.data(),
            set_state_with_rhoi_T);
        m_Pwall = m_thermo.P();
    }

//==============================================================================

    void computeMoleFracfromPartialDens(
        Eigen::VectorXd& v_rhoi, Eigen::VectorXd& v_xi)
    {
        m_thermo.setState(v_rhoi.data(),
                          mv_Twall.data(),
                          set_state_with_rhoi_T);
        v_xi = Eigen::Map<const Eigen::VectorXd>(m_thermo.X(), m_ns);
    }

//==============================================================================

    void computePartialDensfromMoleFrac(
        Eigen::VectorXd& v_xi, Eigen::VectorXd& v_rhoi)
    {
    	v_rhoi = v_xi.cwiseProduct(m_thermo.speciesMw().matrix()) *
    			  m_Pwall / (mv_Twall(pos_T_trans) * Mutation::RU);
    }
//==============================================================================

    void errorEmptyWallProductionTerms() const {
        if (mv_surf_prod.size() == 0) {
            throw LogicError()
            << "At least one wall production term should be provided "
            << "in the input file!";
        }
    }

//==============================================================================

    void errorWallStateNotSet() const {
        if (m_wall_state.isWallStateSet() == 0) {
            throw LogicError()
            << "The wall state must have been set!";
        }
    }

//==============================================================================
private:
    Mutation::Thermodynamics::Thermodynamics& m_thermo;

    SurfaceProperties& m_surf_props;
    WallState& m_wall_state;
    
    std::vector<WallProductionTerms*> mv_surf_prod;

    DiffusionVelocityCalculator* mp_diff_vel_calc;
    MassBlowingRate* mp_mass_blowing_rate;

    // VARIABLES FOR SOLVER
    const size_t m_ns;
    const size_t m_nE;
    Eigen::VectorXd mv_Twall;
    double m_Pwall;
    Eigen::VectorXd mv_rhoi;
    Eigen::VectorXd mv_X;
    Eigen::VectorXd mv_dX;
    Eigen::VectorXd mv_f;
    Eigen::MatrixXd mv_jac;
    double m_pert;
    double m_X_unpert;
    Eigen::VectorXd mv_f_unpert;
    Eigen::VectorXd mv_sep_mass_prod_rate;

    const int pos_T_trans;
    const int set_state_with_rhoi_T;
}; // class SurfaceBalanceSolverMass

ObjectProvider<
    SurfaceBalanceSolverMass, SurfaceBalanceSolver>
    surface_balance_solver_mass_gamma("gamma");
 
    } // namespace GasSurfaceInteraction
} // namespace Mutation
