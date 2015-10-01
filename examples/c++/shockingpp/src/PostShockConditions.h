#ifndef POSTSHOCKCONDITIONS_H
#define POSTSHOCKCONDITIONS_H

#include "mutation++.h"

#include "InputData.h"

class PostShockConditions {
public:
    virtual ~PostShockConditions(){}
    void computationalError(){}

    virtual std::vector<double> getPostShockMoleFrac() = 0;

protected:
    virtual void computePostShockConditions() = 0;

};

class PostShockConditionsColdGas : public PostShockConditions{
public:
    PostShockConditionsColdGas( Mutation::Mixture* const lp_mixture, InputData* const lp_input_data );

    std::vector<double> getPostShockMoleFrac(){ return v_X_post; }
    std::vector<double> getPostShockTemperature(){ return v_temp_post; }
    double getVsminusV2(){ return m_VsmV2; }
    double getPostShockDensity(){ return m_rho_post; }
    double getPostMomentumDensity(){ return getMomentumDensity(); }

protected:
    void computePostShockConditions();
    Mutation::Mixture* mp_mixture;
    InputData* mp_input_data;

private:
    int m_ns;
    int m_nEnergyEqns;
    
    double m_P_post;
    double m_V_post;
    double m_VsmV2;
    std::vector<double> v_rhoi_post;
    std::vector<double> v_temp_post;

    std::vector<double> v_X_post;
    double m_rho_post;

    int set_state_rhoi_T;
};

#endif // POSTSHOCKCONDITIONS_H
