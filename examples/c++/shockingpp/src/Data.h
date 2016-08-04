#ifndef DATA_H
#define DATA_H

#include <boost/numeric/odeint.hpp>
#include <vector>

typedef boost::numeric::ublas::vector<double> vector_type;

class Data {
public:
    virtual ~Data(){}
    
    virtual vector_type getInitialState() = 0;
    virtual double mdot(){ return 0.0; }
    virtual int nEquations() const = 0;
    
    virtual std::vector<double> getPartialDensities() const = 0;
    virtual double getPressure() const = 0;
    virtual double getVelocity() const = 0;
    virtual double getTTrans() const = 0;
    virtual std::vector<double> getTemperatures() const = 0;
    virtual std::vector<double> getMassFractions() const = 0;
    
    virtual void setPressure(const double& l_P) = 0;
    virtual void setTTrans(const double& l_T) = 0;
    virtual void setTemperatures(const std::vector<double>& l_T) = 0;
    virtual void setMassFractions(const std::vector<double>& l_yi) = 0;
   
    // S H O C K I N G
    virtual void setVelocity(const double& l_V) { errorFunctionCall("setVelocity()"); };
    virtual void buildState() { errorFunctionCall("buildState()"); };
    virtual void buildStatePostShock() { errorFunctionCall("buildStatePostShock()"); };

    // L A R S E N
    virtual void readDataFile() { errorFunctionCall("readDataFile()"); };
    virtual void setCurrentStep(int pos) { errorFunctionCall("setCurrentStep()"); };
    virtual void getStepValues(double x[], double u[], double rho[]) { errorFunctionCall("getStepValues()"); };
    virtual int  getStrNele() { errorFunctionCall("getStrNele()"); };
   
private:
    void errorFunctionCall(const char * funcName){
      std::cerr << " In: Data.h - " << funcName << " was called, however it is not defined for ";
      std::cerr << "the current problem_type!\n"; 
      exit(1); 
    }
};  

#endif /* DATA_H */
