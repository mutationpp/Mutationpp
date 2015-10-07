#ifndef SURFACEBALANCESOLVER_H
#define SURFACEBALANCESOLVER_H

#include <Eigen/Dense>

#include "NewtonSolver.h"
#include "Thermodynamics.h"
#include "Transport.h"
#include "Utilities.h"

#include "DiffusionVelocityCalculator.h"
#include "MassBlowingRate.h"
#include "SurfaceDescription.h" 
#include "WallProductionTerms.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class SurfaceBalanceSolver : public Mutation::Numerics::NewtonSolver<Eigen::VectorXd, SurfaceBalanceSolver> {

public:
    SurfaceBalanceSolver( Mutation::Thermodynamics::Thermodynamics& l_thermo, Mutation::Transport::Transport& l_transport, const std::string& l_gsi_mechanism, const Mutation::Utilities::IO::XmlElement& l_node_diff_model, const Mutation::Utilities::IO::XmlElement& l_node_prod_terms, SurfaceDescription& l_surf_descr );
    ~SurfaceBalanceSolver();

    void setWallState( const double* const l_mass, const double* const l_energy, const int state_variable );
    void computeGSIProductionRate( Eigen::VectorXd& lv_mass_prod_rate );

    void setDiffusionModel( const Eigen::VectorXd& lv_mole_frac_edge, const double& l_dx );
    void solveSurfaceBalance( const Eigen::VectorXd& lv_rhoi, const Eigen::VectorXd& lv_T );

    void updateFunction( Eigen::VectorXd& lv_mole_frac_wall );
    void updateJacobian( Eigen::VectorXd& lv_mole_frac_wall );
    Eigen::VectorXd& systemSolution();
    double norm();

private:
    Mutation::Thermodynamics::Thermodynamics& m_thermo;

    SurfaceDescription& m_surf_descr;
    
    std::vector<WallProductionTerms*> v_surf_prod;

    DiffusionVelocityCalculator* mp_diff_vel_calc;
    MassBlowingRate* mp_mass_blowing_rate;

    void addSurfaceProductionTerm( WallProductionTerms* lp_wall_prod_term );

    // Solver Functions
    void saveUnperturbedPressure( Eigen::VectorXd& lv_rhoi );
    void computeMoleFracfromPartialDens( Eigen::VectorXd& lv_rhoi, Eigen::VectorXd& lv_xi );
    void computePartialDensfromMoleFrac( Eigen::VectorXd& lv_xi, Eigen::VectorXd& lv_rhoi );

    // Error Functions
    void errorEmptyWallProductionTerms() const;
    void errorWallStateNotSet() const;

    // VARIABLES FOR SOLVER
    const size_t m_ns;
    double m_Twall;
    double m_Pwall;
    Eigen::VectorXd v_rhoi;
    Eigen::VectorXd v_work;
    Eigen::VectorXd v_X;
    Eigen::VectorXd v_dX;
    Eigen::VectorXd v_f;
    Eigen::MatrixXd v_jac;
    double m_pert;
    double m_X_unpert;
    Eigen::VectorXd v_f_unpert;

    const int set_state_with_rhoi_T;

    // TEMPORARY
    Eigen::VectorXd v_sep_mass_prod_rate;

};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // SURFACEBALANCESOLVER_H
