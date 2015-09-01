#ifndef SURFACEBALANCESOLVER_H
#define SURFACEBALANCESOLVER_H

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

class SurfaceBalanceSolver : public Mutation::Numerics::NewtonSolver<Mutation::Numerics::RealVector, SurfaceBalanceSolver> {

public:
    SurfaceBalanceSolver( Mutation::Thermodynamics::Thermodynamics& l_thermo, Mutation::Transport::Transport& l_transport, const std::string& l_gsi_mechanism, const Mutation::Utilities::IO::XmlElement& l_node_diff_model, const Mutation::Utilities::IO::XmlElement& l_node_prod_terms, SurfaceDescription& l_surf_descr );
    ~SurfaceBalanceSolver();

    void setWallState( const double* const l_mass, const double* const l_energy, const int state_variable );
    void computeGSIProductionRate( Mutation::Numerics::RealVector& lv_mass_prod_rate );

    void setDiffusionModel( const Mutation::Numerics::RealVector& lv_mole_frac_edge, const double& l_dx );
    void solveSurfaceBalance( const Mutation::Numerics::RealVector& lv_rhoi, const Mutation::Numerics::RealVector& lv_T );

    void updateFunction( Mutation::Numerics::RealVector& lv_mole_frac_wall );
    void updateJacobian( Mutation::Numerics::RealVector& lv_mole_frac_wall );
    Mutation::Numerics::RealVector& systemSolution();
    double norm();

private:
    Mutation::Thermodynamics::Thermodynamics& m_thermo;

    SurfaceDescription& m_surf_descr;
    
    std::vector<WallProductionTerms*> v_surf_prod;

    DiffusionVelocityCalculator* mp_diff_vel_calc;
    MassBlowingRate* mp_mass_blowing_rate;

    void addSurfaceProductionTerm( WallProductionTerms* lp_wall_prod_term );

    // Solver Functions
    void saveUnperturbedPressure( Mutation::Numerics::RealVector& lv_rhoi );
    void computeMoleFracfromPartialDens( Mutation::Numerics::RealVector& lv_rhoi, Mutation::Numerics::RealVector& lv_xi );
    void computePartialDensfromMoleFrac( Mutation::Numerics::RealVector& lv_xi, Mutation::Numerics::RealVector& lv_rhoi );

    // Error Functions
    void errorEmptyWallProductionTerms() const;
    void errorWallStateNotSet() const;

    // VARIABLES FOR SOLVER
    const size_t m_ns;
    double m_Twall;
    double m_Pwall;
    Mutation::Numerics::RealVector v_rhoi;
    Mutation::Numerics::RealVector v_work;
    Mutation::Numerics::RealVector v_X;
    Mutation::Numerics::RealVector v_dX;
    Mutation::Numerics::RealVector v_f;
    Mutation::Numerics::RealMatrix v_jac;
    double m_pert;
    double m_X_unpert;
    Mutation::Numerics::RealVector v_f_unpert;

    const int set_state_with_rhoi_T;

    // TEMPORARY
    Mutation::Numerics::RealVector v_sep_mass_prod_rate;

};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // SURFACEBALANCESOLVER_H
