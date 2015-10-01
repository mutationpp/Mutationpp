#include "SolveSystem.h"

SolveSystem::SolveSystem( Mutation::Mixture* lp_mixture, PostShockConditionsColdGas& l_post_shock_conditions, MeshOptions* lp_mesh_options  )
                        : n_momentum_eq( 1 ),
                          m_ns( lp_mixture->nSpecies() ),
                          m_nEnergyEqns( lp_mixture->nEnergyEqns() ),
                          x(m_ns + m_nEnergyEqns + 1, 0.0),
                          mp_mesh_options( lp_mesh_options ) {

    // Initial Solution
    lp_mixture->convert<Mutation::Thermodynamics::X_TO_Y>( &l_post_shock_conditions.getPostShockMoleFrac()[0] , &x(0) );

    // ADD A TOLERENCE;

    x( m_ns ) = l_post_shock_conditions.getVsminusV2();
    for (int i_nEn = 0; i_nEn < m_nEnergyEqns; ++i_nEn ){
        x( m_ns + 1 + i_nEn ) = l_post_shock_conditions.getPostShockTemperature()[ i_nEn ];
    }

}

void SolveSystem::solve(){

//        size_t num_of_steps = integrate_times( make_dense_output<rosenbrock4<double>>(1.0e-8,1.0e-8), 
//                                               make_pair( m_computeRHS(), m_computeJacobian() ) ,
//                                               &lp_mesh_options->getMesh()[0], &lp_mesh_options->getMesh()[ &lp_mesh_options->getNumberMeshPoints - 1 ],
//                                               &lp_mesh_options->getMeshdx);

}

//=============================================================================================================

RHS::RHS( Mutation::Mixture* lp_mixture, PostShockConditionsColdGas& l_post_shock_conditions ) 
         : mp_mixture( lp_mixture ),
           m_ns( lp_mixture->nSpecies() ){ }

RHS::operator(  const vector_type& x , vector_type& dxdt , double t ){


    // use state vector x to set the state in M++
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    // Compute a, b, c, d, e, f
    //
    // Get dxdt

    for ( int i_ns = 0; i_ns < m_ns; ++i_ns ){
        dxdt[ i_ns ] = 0.E0;
    }

    dxdt[ i_ns + 1 ] = 0.0E0;    
    dxdt[ i_ns + 2 ] = 0.0E0;

    //size_t iUIndex = ns;
    //size_t iTIndex = ns+1;

    double* Wdot = new double [ns];
    double* Yi = new double [ns];
    double rhoi [ns];
    double* Cp = new double [ns];
    double* H  = new double [ns];
    
    double T = x[iTIndex];
    double u = x[iUIndex];
    double rho = mdot/u;
    
    for (size_t iSpecies = 0; iSpecies<ns; iSpecies++) {
        Yi[iSpecies] = x[iSpecies];
        rhoi[iSpecies] = rho*Yi[iSpecies];
    }
    double Mm = mix->mixtureMwMass(Yi);
    
    mix->speciesCpOverR(T, Cp);
    mix->speciesHOverRT(T, H);
    
    double p = rho*RU*T/Mm;
    
    mix->setState(rhoi,&T,1);
    mix->netProductionRates(Wdot);
    
    double a = mdot*u/(RU*T) - p/(RU*T);
    double b = 0.0;
    double c = 0.0;
    double d = mdot*u;
    double e = 0.0;
    double f = 0.0;
    
    for (int iSpecies = 0; iSpecies<ns; iSpecies++){
        double Mmi = mix->speciesMw(iSpecies);
        
        b += (mdot/T)*(x[iSpecies]/Mmi);
        c -= Wdot[iSpecies]/Mmi;
        e += mdot*(RU/Mmi)*(x[iSpecies]*(Cp[iSpecies]));
        f -= Wdot[iSpecies]*H[iSpecies]*(RU/Mmi)*T;
        
        dxdt[ iSpecies ] = Wdot[iSpecies]*overmdot;
    }
    
    double overdetA = 1.0/(a*e-b*d);
    dxdt[iUIndex] = (e*c-b*f)*overdetA;
    dxdt[iTIndex] = (a*f-c*d)*overdetA;
    
}

//=============================================================================================================


