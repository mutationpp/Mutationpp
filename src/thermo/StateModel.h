#ifndef THERMO_STATE_MODEL_H
#define THERMO_STATE_MODEL_H

/**
 * Base class for all state models.  A mixture state is completely determined
 * when enough thermodynamic values are combined with the mixture composition
 * such that all other mixture quantities can be successfully determined.  For 
 * example, with a simple mixture modelled in thermal equilibrium, the mixture
 * temperature, pressure, and species mole fractions are enough to fully
 * determine all other mixture quantities.
 *
 * A state model describes how a particular mixture is to be modelled.  For
 * example, a mixture could be modelled using a single temperature or multiple
 * temperatures that represent the various energy modes in the mixture.  This
 * base class provides the framework for defining new state models above the 
 * simple T,P model.
 */
class StateModel
{
public:

    typedef int ARGS;
    
    /**
     * Constructor.
     *
     * @param ns - number of species
     */
    StateModel(int ns)
        : m_ns(ns)
    {
        init();
    }
    
    /**
     * Destructor.
     */
    virtual ~StateModel()
    {
        delete [] mp_T;
        delete [] mp_P;
        delete [] mp_X;
    }
    
    /**
     * Sets the current mixture state via temperature(s), pressure(s), and 
     * species mole fractions.
     *
     * @param T - temperature array
     * @param P - pressure array
     * @param X - species mole fraction array
     */
    virtual void setStateTPX(
        const double* const T, const double* const P, const double* const X);
    
    /**
     * Sets the current mixture state using thermal equilibrium (ie: all the 
     * temperatures are equal.
     */
    void setStateTPX(const double T, const double P, const double* X);
    
    /**
     * Returns the mixture translational temperature.
     */
    virtual double T() const {
        return mp_T[0];
    }
    
    /**
     * Returns the mixture vibrational temperature.
     */
    virtual double Tv() const {
        return mp_T[0];
    }
    
    /**
     * Returns the mixture electron temperature.
     */
    virtual double Te() const {
        return mp_T[0];
    }
    
    /**
     * Returns the mixture rotational temperature.
     */
    virtual double Tr() const {
        return mp_T[0];
    }
    
    /**
     * Returns the mixture electronic temperature.
     */
    virtual double Tel() const {
        return mp_T[0];
    }
    
    /**
     * Returns the mixture static pressure.
     */
    virtual double P() const {
        return mp_P[0];
    }
    
    /**
     * Returns the species mole fractions.
     */
    const double* const X() const {
        return mp_X;
    }

protected:

    /**
     * Allocates storage for the state model and initializes the memory to zero.
     */
    virtual void init();
    
    /**
     * Returns the number of temperatures used in this model.
     */
    virtual int nT() const {
        return 1;
    }
    
    /**
     * Returns the number of pressures used in this model.
     */
    virtual int nP() const {
        return 1;
    }

protected:

    int m_ns;

    double* mp_T;
    double* mp_P;
    double* mp_X;
    
}; // class StateModel


#endif // THERMO_STATE_MODEL_H

