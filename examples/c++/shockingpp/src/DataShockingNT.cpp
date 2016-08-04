#include "DataShockingNT.h"

DataShockingNT::DataShockingNT(Mutation::Mixture& l_mix)
                        : m_mix(l_mix),
                          n_sp(l_mix.nSpecies()),
                          n_meq(1), // An 1D solver can have 1 component of velocity
                          n_eneq(l_mix.nEnergyEqns()),
                          n_eq(n_sp+n_meq+n_eneq),
                          pos_V(n_sp),
                          pos_T(n_sp+n_meq),
                          m_P(0.0),
                          m_V(0.0),
                          m_rho(0.0),
                          m_mdot(0.0),
                          v_rhoi(n_sp, 0.0),
                          v_xi(n_sp, 0.0),
                          v_yi(n_sp, 0.0),
                          v_T(n_eneq, 0.0),
                          v_X(n_eq, 0.0)
{}

DataShockingNT::DataShockingNT(Mutation::Mixture& l_mix, const std::string& l_init)
                        : m_mix(l_mix),
                          n_sp(l_mix.nSpecies()),
                          n_meq(1), // An 1D solver can have 1 component of velocity
                          n_eneq(l_mix.nEnergyEqns()),
                          n_eq(n_sp+n_meq+n_eneq),
                          pos_V(n_sp),
                          pos_T(n_sp+n_meq),
                          m_P(0.0),
                          m_V(0.0),
                          m_rho(0.0),
                          m_mdot(0.0),
                          v_rhoi(n_sp, 0.0),
                          v_xi(n_sp, 0.0),
                          v_yi(n_sp, 0.0),
                          v_T(n_eneq, 0.0),
                          v_X(n_eq, 0.0)
{
    // Parse info

    // Or for the time being set them arbitrary!
    m_P = 5.2;      // Pa
    for (int i_en = 0; i_en < n_eneq; i_en++){
        v_T[i_en] = 210.0; // K Free stream in equilibrium
    }
    m_V = 11310.0;  // m/s

    errorStateNotSetProperly();
    buildState();
}

DataShockingNT::~DataShockingNT() {}

void DataShockingNT::buildState(){

    m_mix.equilibriumComposition(v_T[0], m_P, &v_xi[0]);
    m_mix.convert<Mutation::Thermodynamics::X_TO_Y>(&v_xi[0], &v_yi[0]);
    m_rho = m_mix.density(v_T[0], m_P, &v_xi[0]);

    // Convert yi to rhoi
    for (int i_sp = 0; i_sp < n_sp; ++i_sp){
        v_rhoi[i_sp] = m_rho *v_yi[i_sp];
    }

    // compute rho*u
    m_mdot = m_rho * m_V;

    fillStateVectorX();
}

void DataShockingNT::buildStatePostShock(){

    m_mix.convert<Mutation::Thermodynamics::Y_TO_X>(&v_yi[0], &v_xi[0]);
    m_rho = m_mix.density(v_T[0], m_P, &v_xi[0]);

    // compute rho*u
    // Important Note:
    // m_mdot is slightly changed after the shock...
    // Taking it from the pre shock and computing rho that way might be more consistent
    m_mdot = m_rho * m_V;

    fillStateVectorX();
}

void DataShockingNT::fillStateVectorX(){

    for (int i_sp = 0; i_sp < n_sp; ++i_sp){
        v_X[i_sp] = v_yi[i_sp];
    }
    v_X[pos_V] = m_V;
    for (int i_en = 0; i_en < n_eneq; ++i_en){
        v_X[pos_T+i_en] = v_T[i_en];
    }
}

void DataShockingNT::errorStateNotSetProperly(){
    // if (){
    //     std::cerr << " not set in the input file!" << std::endl;
    // }

}
