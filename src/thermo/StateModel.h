#ifndef THERMO_STATE_MODEL_H
#define THERMO_STATE_MODEL_H

#include "Thermodynamics.h"
#include "Kinetics.h"
#include "TransferModel.h"
#include "Transport.h"

namespace Mutation {
    namespace Thermodynamics {


/**
 * Simple interface for any object who wants to be updated whenever the a
 * StateModel changes its state.
 */
//class StateModelWatcher
//{
//protected:
//    void stateChanged() = 0;
//}


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

    typedef const Thermodynamics& ARGS;
    
    /**
     * Constructor initializes state.
     *
     * @param ns - number of species
     */
    StateModel(ARGS thermo, const int nenergy, const int nmass)
        : m_thermo(thermo), m_nenergy(nenergy), m_nmass(nmass)
    {
        m_T = m_Tr = m_Tv = m_Tel = m_Te = 300.0;
        m_P = 0.0;
        mp_X = new double [m_thermo.nSpecies()];
        for (int i = 0; i < thermo.nSpecies(); ++i)
            mp_X[i] = 0.0;
    }
    
    /**
     * Destructor.
     */
    virtual ~StateModel()
    {
        delete [] mp_X;
    }
    
//    void addWatcher(StateModelWatcher* p_watcher);
    
    inline int nEnergyEqns() const { return m_nenergy; }
    inline int nMassEqns() const { return m_nmass; }

    /**
     * Abstract method to set the state of the mixture which is dependent on the
     * type of the concrete StateModel class.
     *
     * @param p_mass - The "mass" vector (depends on concrete class)
     * @param p_energy - The "energy" vector (depends on concrete class)
     * @param vars - Index representing which variable set is given in the mass
     * and energy vectors
     */
    virtual void setState(
        const double* const p_mass, const double* const p_energy,
        const int vars = 0) = 0;
    
    /**
     * Returns the mixture translational temperature.
     */
    inline double T() const {
        return m_T;
    }
    
    /**
     * Returns the mixture vibrational temperature.
     */
    inline double Tv() const {
        return m_Tv;
    }
    
    /**
     * Returns the mixture electron temperature.
     */
    inline double Te() const {
        return m_Te;
    }
    
    /**
     * Returns the mixture rotational temperature.
     */
    inline double Tr() const {
        return m_Tr;
    }
    
    /**
     * Returns the mixture electronic temperature.
     */
    inline double Tel() const {
        return m_Tel;
    }
    
    /**
     * Returns the mixture static pressure.
     */
    inline double P() const {
        return m_P;
    }
    
    /**
     * Returns the species mole fractions.
     */
    inline const double* const X() const {
        return mp_X;
    }
    
    /**
     * Creates a new TransferModel object which can compute the energy transfer
     * source terms for this state model.
     */
    virtual Mutation::Transfer::TransferModel* createTransferModel(
        Thermodynamics& thermo,
        Mutation::Transport::Transport& transport,
        Mutation::Kinetics::Kinetics& kinetics)
    {
        return NULL;
    }
    
//protected:

//    void notifyWatchers();
    
protected:

    const Thermodynamics& m_thermo;
    const int m_nenergy;
    const int m_nmass;
    
    double m_T;
    double m_Tv;
    double m_Tr;
    double m_Tel;
    double m_Te;
    double m_P;
    
    double* mp_X;
    
//private:
    
//    std::list<StateModelWatcher*> m_watcher_list;
    
}; // class StateModel


    } // namespace Thermodynamics
} // namespace Mutation

#endif // THERMO_STATE_MODEL_H

