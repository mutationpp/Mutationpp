#ifndef INPUTDATA_H
#define INPUTDATA_H

#include <fstream>
#include <string>
#include <vector>

#include "mutation++.h"

class InputData{
public:
    InputData( Mutation::Mixture* const p_mixture, std::ifstream& l_input_file );
    ~InputData(){}

    double getPreShockPressure() { return m_P; }
    double getPreShockVelocity() { return m_Vs; }
    double getMomentumDensity(){ return m_rho * m_Vs; }

    std::vector<double> getPreShockTemperature() const { return v_temp; }
    std::vector<double> getPreShockPartialDensities() { return v_rhoi; }
    std::vector<double> getPreShockMoleFrac() { return v_X; }
    
private:
    std::string line;

    int m_ns;
    int m_nEnergyEqns;
    
    bool set_pressure;
    bool set_temperature;
    bool set_velocity;

    double m_P;
    double m_Vs;
    std::vector<double> v_rhoi;
    std::vector<double> v_temp;

    std::vector<double> v_X;
    double m_rho;

    int set_state_rhoi_T;

};

#endif // INPUTDATA_H
