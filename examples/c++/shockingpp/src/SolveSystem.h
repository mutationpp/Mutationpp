#ifndef SOLVESYSTEM_H
#define SOLVESYSTEM_H

#include<vector>

#include "MeshOptions.h"
#include "PostShockConditions.h"

#include "mutation++.h"

#include <boost/numeric/odeint.hpp>

typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;

class RHS{
public:
    RHS( Mutation::Mixture* lp_mixture, PostShockConditionsColdGas& l_post_shock_conditions );

    void operator()( const vector_type& x, vector_type& dxdt, double t );

private:
    Mutation::Mixture* mp_mixture;

};

class Jacobian{
public:
    void operator()( const vector_type& x, matrix_type& J, const double& t, vector_type& dfdt ){ }

};

class SolveSystem{
public:
    SolveSystem( Mutation::Mixture* lp_mixture, PostShockConditionsColdGas& l_post_shock_conditions, MeshOptions* lp_mesh_options );
    ~SolveSystem(){ }

    void solve();

private:
    vector_type x;

    const int n_momentum_eq;
    const int m_ns;
    const int m_nEnergyEqns;

    RHS m_rhs;
    Jacobian m_jacobian;
      
    MeshOptions* mp_mesh_options;

};

#endif // SOLVESYSTEM_H
