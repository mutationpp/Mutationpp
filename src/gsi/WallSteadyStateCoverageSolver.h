#ifndef WALLSTEADYSTATECOVERAGESOLVER_H
#define WALLSTEADYSTATECOVERAGESOLVER_H

#include "NewtonSolver.h"
#include "QR.h"

#include "Thermodynamics.h"

#include "SurfaceProperties.h"

namespace Mutation{
    namespace gsi{

class WallSteadyStateCoverageSolver : public Mutation::Numerics::NewtonSolver<Mutation::Numerics::RealVector, WallSteadyStateCoverageSolver>{
public:

WallSteadyStateCoverageSolver( const Mutation::Thermodynamics::Thermodynamics& l_thermo, WallState& l_wall_state );
~WallSteadyStateCoverageSolver();

void updateFunction();
void updateJacobian( Mutation::Numerics::RealVector& );
Mutation::Numerics::RealVector& systemSolution();
double norm();

private:
//Mutation::Thermodynamics::Thermodynamics& m_thermo;
WallState m_wall_state;

Mutation::Numerics::RealMatrix v_jacobian;
Mutation::Numerics::RealMatrix v_residual_function;

void computeProductionRates();
void convertWallNumberDensitiestoNumberDensityFractions();
void convertWallNumberDensityFractionstoNumberDensities();

};

    } // namespace gsi
} // namespace Mutation

#endif // WALLSTEADYSTATECOVERAGESOLVER_H
