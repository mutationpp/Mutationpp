#ifndef SOLVESYSTEM_H
#define SOLVESYSTEM_H

#include<vector>

#include "MeshOptions.h"
#include "PostShockConditions.h"

#include "mutation++.h"

#include <boost/numeric/odeint.hpp>

typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;

// THIS WILL BE PLACED IN A SEPARATE FILE. This is going to be the base class... Derive for ST, MT, and MT+electrons
class Problem {
public:
    Problem( const Mutation::Mixture* const lp_mixture,  const PostShockConditionsColdGas& l_post_shock_conditions ){}
    
    void problem( const vector_type& x, vector_type& dx_dt, const double& t ){ }

private:

};
// TILL HERE

class RHS{
public:
    RHS( Problem& l_problem );

    void operator()( const vector_type& x, vector_type& dxdt, double t );

private:
    Problem& m_problem;

};

// class Jacobian{
// public:
//     Jacobian( Problem& l_problem );
//
//     void operator()( const vector_type& x, matrix_type& J, const double& t, vector_type& dfdt );
// private:
//
//     Problem& m_problem;
//
// };

class SolveSystem{
public:
    SolveSystem( const Mutation::Mixture* const lp_mixture, PostShockConditionsColdGas& l_post_shock_conditions, MeshOptions* lp_mesh_options );
    ~SolveSystem(){ }

    void solve();

private:
    vector_type x; // This does not necessary need to be here! Make it private member of RHS and Jacobian?

    const int m_nMomentumEqns;
    const int m_ns;
    const int m_nEnergyEqns;

    Problem m_problem; // RENAME!!! STRATEGY. This will be changing depending on the physicochemical model. This is going to be an abstract interface. This is of course going to be replaced by a factory.
//    RHS m_rhs;
    Jacobian m_jacobian;
      
    MeshOptions* mp_mesh_options;

};

#endif // SOLVESYSTEM_H
