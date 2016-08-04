#ifndef DATALARSENTTv_H
#define DATALARSENTTv_H

#include <string>
#include "mutation++.h"

#include "Data.h"

class DataLarsenTTv: public Data {
public:
    DataLarsenTTv(Mutation::Mixture& l_mix);
    ~DataLarsenTTv();
    
    vector_type getInitialState();
    double mdot(){ return m_mdot; }
    int nEquations() const { return n_eq; }
    
    std::vector<double> getPartialDensities() const { return v_rhoi; }
    double getPressure() const { return m_P; }
    double getVelocity() const { return m_V; }
    double getTTrans() const { return v_T[0]; }
    std::vector<double> getTemperatures() const { return v_T; }
    std::vector<double> getMassFractions() const { return v_yi; }
    
    void setPressure(const double& l_P){ m_P = l_P; }
    void setVelocity(const double& l_V){ m_V = l_V; }
    void setTTrans(const double& l_T){ v_T[0] = l_T;}
    void setTemperatures(const std::vector<double>& lv_T){ v_T = lv_T; }
    void setMassFractions(const std::vector<double>& l_yi){ v_yi = l_yi; }
    
    void readDataFile();     // reads external data file and fills internal vectors
    void setCurrentStep(int pos);
    void getStepValues(double x_now[], double u_now[], double rho_now[]);  // returns values at the set step
    int  getStrNele();

private:
    Mutation::Mixture& m_mix;
    const size_t n_sp;
    const size_t n_meq;
    const size_t n_eneq;
    const size_t n_eq;
    const size_t pos_V;
    const size_t pos_T;

    double m_P;
    double m_V;
    double m_rho;
    double m_mdot;
    std::vector<double> v_rhoi;
    std::vector<double> v_xi;
    std::vector<double> v_yi;
    std::vector<double> v_T;
    std::vector<double> v_Tv;

    vector_type v_X;

    void errorStateNotSetProperly();
    inline void fillStateVectorX();

    // ------------>  SPECIFICALLY FOR LARSEN  <------------------
    // Those vectors store the whole streamline
    std::vector<double> x_ext, u_ext, rho_ext, T_ext, Tv_ext;

    // Initial condition
    std::vector<double> Y0_ext;

    // Current step
    int currentStep;

};

#endif /* DATALARSENTTv_H */
