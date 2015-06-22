#ifndef WALLSOLVER_H
#define WALLSOLVER_H

#include "NewtonSolver.h"
#include "QR.h"
#include "Vector.h"

#include "Thermodynamics.h"
#include "Transport.h"

#include "CatalysisRateManager.h"
#include "SurfaceProperties.h"

namespace Mutation{
    namespace gsi{

class WallSolver : public Mutation::Numerics::NewtonSolver<Mutation::Numerics::RealVector, WallSolver>{
public:

    /**
     * Default Constructor
     */
    WallSolver( Mutation::Thermodynamics::Thermodynamics& thermo,
                Mutation::Transport::Transport& transport,
                WallState* wallstate,
                CatalysisRateManager* ratemanager);

    /**
     * Default Destructor
     */
    ~WallSolver(){}

public: // Private?

    /**
     * @brief updateFunction
     * @param lv_rhoi
     */
    void updateFunction( Mutation::Numerics::RealVector& lv_mole_fractions );

    /**
     * @brief updateJacobian
     * @param lv_rhoi
     */
    void updateJacobian( Mutation::Numerics::RealVector& lv_mole_fractions );

    /**
     * @brief systemSolution
     * @return
     */
    Mutation::Numerics::Vector<double>& systemSolution();
// Vector<double>& ??? -> RealVector

    /**
     * @brief norm
     * @return
     */
    double norm();

    /**
     * 
     */
    void setDataSolver(const double * const l_rhoi_edge, const double * const l_Twall, const double& l_dx_gradient_distance_wall_edge );

    void initializeSolutionwithPartialDensities( Mutation::Numerics::RealVector& lv_rhoi_wall );

    void retrieveInitialMolarFractions( Mutation::Numerics::RealVector& l_mole_fraction_wall );
    void returnSolutioninPartialDensities( Mutation::Numerics::RealVector& lv_solution_rhoi_wall );

private:
    /**
     * @brief computeDerivative
     * @param l_rhoi
     */
     void computeMoleFractionGradientatWall( const Mutation::Numerics::RealVector& l_mole_fractions_wall );
 
     void saveUnperturbedPressure( Mutation::Numerics::RealVector& lv_rhoi_wall);
     void computeMoleFractionsfromPartialDensities( Mutation::Numerics::RealVector& lv_rhoi_wall, Mutation::Numerics::RealVector& lv_mole_fractions );
     void computeMoleFractionsfromPartialDensities(const double * const p_rhoi, Mutation::Numerics::RealVector& lv_mole_fractions );
     void computePartialDensitiesfromMoleFractions( Mutation::Numerics::RealVector& lv_mole_fractions, Mutation::Numerics::RealVector& lv_rhoi_wall );

private:
    Mutation::Transport::Transport& m_transport;
    Mutation::Thermodynamics::Thermodynamics& m_thermo;

    CatalysisRateManager* mp_catalysis_rates;
    WallState* mp_wall_state;

    Mutation::Numerics::RealVector v_rhoi_wall;

    Mutation::Numerics::RealVector v_mole_fractions_wall;
    Mutation::Numerics::RealVector v_mole_fractions_edge;
    Mutation::Numerics::RealVector v_mole_fraction_gradients;
    Mutation::Numerics::RealVector v_delta_mole_fractions;

    double m_Twall;
    double m_Pwall;

    Mutation::Numerics::RealVector v_residual_function;
    Mutation::Numerics::RealVector v_unperturbed_residual_function;
    Mutation::Numerics::RealVector v_diffusion_velocities;
    Mutation::Numerics::RealVector v_wall_production_rates;
    Mutation::Numerics::RealMatrix v_jacobian;

    double m_perturbation;

    double m_dx_gradient_distance_wall_edge;
    double zero_electric_conductivity;

    const int m_ns;

    static const size_t set_state_with_rhoi_Twall = 1;

};

    } // namespace gsi
} // namespace Mutation

#endif // WALLSOLVER_H
