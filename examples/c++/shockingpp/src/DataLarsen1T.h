#ifndef DATALARSEN1T_H
#define DATALARSEN1T_H

#include <string>
#include "mutation++.h"

#include "Data.h"

class DataLarsen1T: public Data {
public:
    DataLarsen1T(Mutation::Mixture& l_mix);
    ~DataLarsen1T();
    
    vector_type getInitialState();
    double mdot(){ return 0.0; }
    int nEquations() const { return n_eq; }
    
    std::vector<double> getPartialDensities() const { return v_rhoi; }
    double getPressure() const { return m_P; }
    double getVelocity() const { return 0.0; }  // dummy function
    double getTTrans() const { return v_T[0]; }
    std::vector<double> getTemperatures() const { return v_T; }
    std::vector<double> getMassFractions() const { return v_yi; }
    
    void setPressure(const double& l_P){ m_P = l_P; }
    void setTTrans(const double& l_T){ v_T[0] = l_T;}
    void setTemperatures(const std::vector<double>& lv_T){ v_T = lv_T; }
    void setMassFractions(const std::vector<double>& l_yi){ v_yi = l_yi; }
    
    void readDataFile();          // reads external data file and fills internal vectors
    void setCurrentStep(int pos); // sets an internal variable to the value 'pos'
    void getStepValues(double x_now[], double u_now[], double rho_now[]);  // returns values at the set step
    int  getStrNele();

private:
    Mutation::Mixture& m_mix;
    const size_t n_sp;
    const size_t n_eneq;
    const size_t n_eq;

    double m_P;
    double m_rho;
    std::vector<double> v_rhoi;
    std::vector<double> v_xi;
    std::vector<double> v_yi;
    std::vector<double> v_T;

    vector_type v_X;

    void errorStateNotSetProperly();
    inline void fillStateVectorX();

    // ------------>  SPECIFICALLY FOR LARSEN  <------------------
    // Those vectors store the whole streamline
    std::vector<double> x_ext, u_ext, rho_ext, T_ext;

    // Initial condition
    std::vector<double> Y0_ext;

    // Current step
    int currentStep;

};

#endif /* DATALARSEN1T_H */
