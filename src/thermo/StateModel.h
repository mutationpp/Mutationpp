#ifndef THERMO_STATE_MODEL_H
#define THERMO_STATE_MODEL_H

#include "Thermodynamics.h"
#include "Kinetics.h"
#include "TransferModel.h"
#include "Transport.h"

namespace Mutation {
    namespace Thermodynamics {

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
     * Constructor.
     *
     * @param ns - number of species
     */
    StateModel(ARGS thermo)
        : m_thermo(thermo)
    {
        mp_X = new double [m_thermo.nSpecies()];
    }
    
    /**
     * Destructor.
     */
    virtual ~StateModel()
    {
        delete [] mp_X;
    }
    
    /**
     * Sets the current mixture state.
     */
    virtual void setState(const double* const v1, const double* const v2) = 0;
    
    /**
     * Returns the mixture translational temperature.
     */
    double T() const {
        return m_T;
    }
    
    /**
     * Returns the mixture vibrational temperature.
     */
    double Tv() const {
        return m_Tv;
    }
    
    /**
     * Returns the mixture electron temperature.
     */
    double Te() const {
        return m_Te;
    }
    
    /**
     * Returns the mixture rotational temperature.
     */
    double Tr() const {
        return m_Tr;
    }
    
    /**
     * Returns the mixture electronic temperature.
     */
    double Tel() const {
        return m_Tel;
    }
    
    /**
     * Returns the mixture static pressure.
     */
    double P() const {
        return m_P;
    }
    
    /**
     * Returns the species mole fractions.
     */
    const double* const X() const {
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
    
protected:

    const Thermodynamics& m_thermo;
    
    double m_T;
    double m_Tv;
    double m_Tr;
    double m_Tel;
    double m_Te;
    double m_P;
    
    double* mp_X;
    
}; // class StateModel


    } // namespace Thermodynamics
} // namespace Mutation

#endif // THERMO_STATE_MODEL_H

