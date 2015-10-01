#include "PostShockConditions.h"

//=============================================================================================================

PostShockConditionsColdGas::PostShockConditionsColdGas( Mutation::Mixture* const lp_mixture, InputData* const lp_input_data ) 
                                                      : m_ns( lp_mixture->nSpecies() ),
                                                        m_nEnergyEqns( lp_mixture->nEnergyEqns() ),
                                                        v_rhoi_post( m_ns, 0.E0 ),
                                                        v_X_post( m_ns, 0.E0 ),
                                                        v_temp_post( m_nEnergyEqns, 0.E0 ),
                                                        set_state_rhoi_T( 1 ),
                                                        mp_mixture( lp_mixture ),
                                                        mp_input_data( lp_input_data ) { 

    computePostShockConditions();

}

//=============================================================================================================

void PostShockConditionsColdGas::computePostShockConditions(){

    mp_mixture->setState( &mp_input_data->getPreShockPartialDensities()[0], &mp_input_data->getPreShockTemperature()[0], set_state_rhoi_T );

    double R_mix = Mutation::RU / mp_mixture->mixtureMw();
    double gamma = mp_mixture->mixtureFrozenGamma();
    double gamma_p1 = gamma + 1.0;

    double c1 = mp_mixture->frozenSoundSpeed();
    double M1 = mp_input_data->getPreShockVelocity() / c1;
    double M1_s = h1 * M1;

    // Mole Fractions do not change across the shock
    v_X_post = mp_input_data->getPreShockMoleFrac();
    
    m_P_post = mp_input_data->getPreShockPressure() * ( 2.0 * gamma * M1_s - gamma + 1.0 ) / gamma_p1 ;
    m_V_post = mp_input_data->getPreShockVelocity() - c1 * 2.0 / gamma_p1 * ( M1 - 1.0 / M1 );
    v_temp_post[0] = mp_input_data->getPreShockTemperature()[0] * ( 2.0 * gamma * M1_s - gamma + 1.0 ) * ( gamma - 1.0 + 2.0 / M1_s ) / ( gamma_p1 * gamma_p1 );

    m_VsmV2 =  mp_input_data->getPreShockVelocity() - m_V_post;

    for ( int i_nEn = 1; i_nEn < m_nEnergyEqns; ++i_nEn ){
        v_temp_post[i_nEn] = mp_input_data->getPreShockTemperature()[ i_nEn ];
    }

    // Density computed from mass conservation
    m_rho_post = getMomentumDensity() / m_V_post;
    p_mixture->convert<Mutation::Thermodynamics::X_TO_Y>( &v_X[0], &v_rhoi_post[0] );

    for ( int i_ns = 0; i_ns < m_ns; +i_ns ){
        v_rhoi_post[i_ns] *= m_rho_post;
    }

}

//=============================================================================================================
