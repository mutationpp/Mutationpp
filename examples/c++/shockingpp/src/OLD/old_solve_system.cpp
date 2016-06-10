#include "SolveSystem.h"

SolveSystem::SolveSystem( const Mutation::Mixture* const lp_mixture, PostShockConditionsColdGas& l_post_shock_conditions, MeshOptions* lp_mesh_options  )
                        : m_nMomentumEqns( 1 ),
                          m_ns( lp_mixture->nSpecies() ),
                          m_nEnergyEqns( lp_mixture->nEnergyEqns() ),
                          x(m_ns + m_nEnergyEqns + m_nMomentumEqns, 0.0),
                          mp_mesh_options( lp_mesh_options ),
                          m_problem( lp_mixture, l_post_shock_conditions ),
                          m_jacobian( m_problem ) {

    // Initial Solution
    lp_mixture->convert<Mutation::Thermodynamics::X_TO_Y>( &l_post_shock_conditions.getPostShockMoleFrac()[0] , &x(0) );

    // ADD A TOLERENCE;

    x( m_ns ) = l_post_shock_conditions.getVsminusV2();
    for (int i_nEn = 0; i_nEn < m_nEnergyEqns; ++i_nEn ){
        x( m_ns + m_nMomentumEqns + i_nEn ) = l_post_shock_conditions.getPostShockTemperature()[ i_nEn ];
    }

}

//=============================================================================================================

void SolveSystem::solve(){
/*
    size_t num_of_steps = integrate_times( make_dense_output<rosenbrock4<double>>(1.0e-8,1.0e-8), 
                                           make_pair( m_computeRHS(), m_computeJacobian() ) ,
                                           &lp_mesh_options->getMesh()[0], &lp_mesh_options->getMesh()[ &lp_mesh_options->getNumberMeshPoints - 1 ],
                                           &lp_mesh_options->getMeshdx);
*/
}

//=============================================================================================================

RHS::RHS( Problem& l_problem ) 
        : m_problem( l_problem ){}

//=============================================================================================================

void RHS::operator()( const vector_type& x , vector_type& dxdt , double t ){
    m_problem.problem( x, dxdt, t );
}

//=============================================================================================================

//   Jacobian::Jacobian( Problem& l_problem )
//                     : m_problem( l_problem ){}
//
//   //=============================================================================================================
//
//   void Jacobian::operator()( const vector_type &x, matrix_type &J, const double &t, vector_type &dfdt ){
//
//   /*
//       vector_type dxdt(numberOfVariables,0.0); // Define them as class members...
//       vector_type x_pert(numberOfVariables,0.0);
//       vector_type dxdt_pert(numberOfVariables,0.0);
//       double pert = 1e-8;
//       matrix_type J_fd(numberOfVariables,numberOfVariables);
//
//       for (int i=0; i<numberOfVariables; i++){
//           x_pert[i] = x[i];
//       }
//       problem(x,dxdt,t,mix,mdot,overmdot);
//
//       for (int i=0; i<numberOfVariables; i++){
//       vector_type dxdt(numberOfVariables,0.0);
//       vector_type x_pert(numberOfVariables,0.0);
//       vector_type dxdt_pert(numberOfVariables,0.0);
//       double pert = 1e-8;
//       matrix_type J_fd(numberOfVariables,numberOfVariables);
//
//       for (int i=0; i<numberOfVariables; i++){
//           x_pert[i] = x[i];
//       }
//       problem(x,dxdt,t,mix,mdot,overmdot);
//
//       for (int i=0; i<numberOfVariables; i++){
//           x_pert[i]*=(1.0+pert);//*=(1.0+pert);//pert;
//           problem(x_pert,dxdt_pert,t,mix,mdot,overmdot);
//           x_pert[i]/=(1.0+pert);//pert;//(1.0+pert);
//           for (int j=0; j<numberOfVariables; j++){
//               J_fd(j,i) = (dxdt_pert[j]-dxdt[j])/((x_pert[i]*pert+1e-16));//pert(x_pert[i]*pert+1e-16);
//               J(j,i) = J_fd(j,i);
//           }
//   //          x_pert[i]-=pert;
//   //          dfdt[i] = 0.0;
//       }
//           x_pert[i]*=(1.0+pert);//*=(1.0+pert);//pert;
//           problem(x_pert,dxdt_pert,t,mix,mdot,overmdot);
//           x_pert[i]/=(1.0+pert);//pert;//(1.0+pert);
//           for (int j=0; j<numberOfVariables; j++){
//               J_fd(j,i) = (dxdt_pert[j]-dxdt[j])/((x_pert[i]*pert+1e-16));//pert(x_pert[i]*pert+1e-16);
//               J(j,i) = J_fd(j,i);
//           }
//   //          x_pert[i]-=pert;
//   //          dfdt[i] = 0.0;
//       }
//   */
//
//   }
