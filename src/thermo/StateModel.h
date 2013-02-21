#ifndef THERMO_STATE_MODEL_H
#define THERMO_STATE_MODEL_H

#include <vector>
#include <cstdlib>

namespace Mutation {
    namespace Thermodynamics {

class StateModelUpdateHandler;

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
 * base class provides the framework for defining new state models beyond the 
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
     * Registers a StateModelUpdateHandler object with this state model so
     * that it will be notified when the mixture state is updated.
     */
    void notifyOnUpdate(StateModelUpdateHandler* p_handler) 
    {
        if (p_handler != NULL)
            m_updateHandlers.push_back(p_handler);
    }
    
    /**
     * Sets the current mixture state via temperature(s), pressure(s), and 
     * species mole fractions.
     *
     * @param T - temperature array
     * @param P - pressure array
     * @param X - species mole fraction array
     */
    void setStateTPX(
        const double* const T, const double* const P, const double* const X)
    {
        // Update the state based on the state model
        stateModelTPX(T, P, X);
        
        // Notify update handlers that the state has been updated
        notifyHandlers();
    }
    
    /**
     * Sets the current mixture state using thermal equilibrium (ie: all the 
     * temperatures are equal.
     */
    void setStateTPX(const double T, const double P, const double* X)
    {
        // Update the state based on the state model
        stateModelTPX(T, P, X);
        
        // Notify update handlers that the state has been updated
        notifyHandlers();
    }
    
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
     * Applies this state model to set the state.
     */
    virtual void stateModelTPX(
        const double* const T, const double* const P, const double* const X);
    
    /**
     * Applies this state model to set the state.
     */
    virtual void stateModelTPX(
        const double T, const double P, const double* const X);
    
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

private:

    void notifyHandlers();
    
private:
    
    std::vector<StateModelUpdateHandler*> m_updateHandlers;
    
}; // class StateModel


/**
 * Interface class which should be implemented when a user class wishes
 * to be told when the StateModel has had its state updated.
 */
class StateModelUpdateHandler
{
public:
    
    /**
     * Destructor.
     */
    virtual ~StateModelUpdateHandler() {};
    
    /**
     * This method is called by the StateModel when the mixture state has been
     * changed.
     */
    virtual void stateUpdated() = 0;
};

    } // namespace Thermodynamics
} // namespace Mutation

#endif // THERMO_STATE_MODEL_H

