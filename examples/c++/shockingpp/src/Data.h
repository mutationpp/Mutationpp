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
	virtual void setVelocity(const double& l_V) = 0;
	virtual void setTTrans(const double& l_T) = 0;
    virtual void setTemperatures(const std::vector<double>& l_T) = 0;
	virtual void setMassFractions(const std::vector<double>& l_yi) = 0;

	virtual void buildState() = 0;
	virtual void buildStatePostShock() = 0;
};

#endif /* DATA_H */
