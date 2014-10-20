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
 * example, with a simple mixture modeled in thermal equilibrium, the mixture
 * temperature, pressure, and species mole fractions are enough to fully
 * determine all other mixture quantities.
 *
 * A state model describes how a particular mixture is to be modeled.  For
 * example, a mixture could be modeled using a single temperature or multiple
 * temperatures that represent the various energy modes in the mixture.  This
 * base class provides the framework for defining new state models.
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
    
    /**
     * Returns the number of energy equations represented by this StateModel.
     */
    inline int nEnergyEqns() const { return m_nenergy; }

    /**
     * Returns the number of mass equations represented by this StateModel.
     */
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
     * Fills an array of the temperatures represented by this StateModel.
     */
    virtual void getTemperatures(double* const p_T) const {
        std::cerr << "getTemperatures()"
                  << " not implemented by this StateModel!" << std::endl;
        std::exit(1);
    }

    /**
     * Returns a vector of length n_species times n_energies with each corresponding
	 * energy per unit mass.  The first n_species values correspond to the total energy
	 * vector.  The remaining n_species vectors correspond to the additional energy modes.
     */
    virtual void getEnergiesMass(double* const p_e) {
        std::cerr << "getEnergiesMass()"
                  << " not implemented by this StateModel!" << std::endl;
        std::exit(1);
    }
    
    /**
     * Returns a vector of length n_species times n_energies with each corresponding
	 * enthalpy per unit mass.  The first n_species values correspond to the total enthalpy
	 * vector.  The remaining n_species vectors correspond to the additional energy modes.
     */
    virtual void getEnthalpiesMass(double* const p_h) {
        std::cerr << "getEnthalpiesMass()"
                  << " not implemented by this StateModel!" << std::endl;
        std::exit(1);
    }
    
    /**
     * Returns a vector of length n_species times n_energies with each corresponding
	 * cp per unit mass.  Each n_species vector corresponds to a temperature in the state
     * model. The first one is associated with the heavy particle translational temperature. 
     */
    
    virtual void getCpsMass(double* const p_Cp){
        std::cerr << "getCpsMass()"
                  << " not implemented by this StateModel!" << std::endl;
        std::exit(1);
    }

    /**
     * Returns a vector of length n_species times n_energies with each corresponding
	 * cp per unit mass.  Each n_species vector corresponds to a temperature in the state
     * model. The first one is associated with the heavy particle translational temperature. 
     */
    
    virtual void getCvsMass(double* const p_Cp){
        std::cerr << "getCvsMass()"
                  << " not implemented by this StateModel!" << std::endl;
        std::exit(1);
    }

	/**
	 * Assigns a unique temperature to each energy mode.
     */
    virtual void getTagModes(int* const p_tag) {
        std::cerr << "getTagModes()"
                  << " not implemented by this StateModel!" << std::endl;
        std::exit(1);
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
/*    virtual Mutation::Transfer::TransferModel* createTransferModel(
        Thermodynamics& thermo,
        Mutation::Transport::Transport& transport,
        Mutation::Kinetics::Kinetics& kinetics)
    {
        return NULL;
    } */ // Probably will not be used.
    
    /**
     * Initializes the energy transfer terms that will be used by each State Model.
     */
    virtual void initializeTransferModel(
        Thermodynamics& thermo,
        Mutation::Transport::Transport& transport,
        Mutation::Kinetics::Kinetics& kinetics) {}
    
    /**
     * This functions provides the total energy transfer source terms
     */
    virtual void energyTransferSource(double* const omega){}
    
protected:

    /**
     * Solves the general form of an energy equation
     * \f[ f(T) = T \left[\sum_i \tilde{\rho}_i \left(\frac{H_i(T)}{R_uT}\right)
     *            + \alpha\right] - \frac{\rho e}{R_u} = 0 \f]
     * for the temperature given the \f$C_{p,i}/R_u\f$ and \f$H_i/R_uT\f$
     * function providers using a Newton-Rhapson iterative procedure.  This
     * function assumes that the mp_X array is filled with the species molar
     * densities (concentrations) prior to being called. The convergence
     * criteria for the Newton iterations is taken to be
     * \f[ |f(T)| < \epsilon_r\left|\frac{\rho e}{R_u}\right| + \epsilon_a \f]
     * where \f$\epsilon_a\f$ and \f$\epsilon_r\f$ are absolute and relative
     * tolerances respectively.  Note that if the total energy equation is
     * being solved, \f$\alpha\f$ should be set to the negative mixture
     * concentration, \f$-\sum \tilde{\rho}_i\f$. For internal energy equations,
     * \f$\alpha = 0\f$.
     *
     * @param cp        class providing species specific heats
     * @param h         class providing species enthalpies
     * @param rhoe      value of \f$\rho e\f$
     * @param T         initial guess on input, output is the solution
     * @param p_work    work array, at least number of species long
     * @param alpha     arbitrary temperature coefficient in energy equation
     * @param atol      absolute tolerance, \f$\epsilon_a\f$, on \f$f\f$
     * @param rtol      relative tolerance, \f$\epsilon_r\f$, on \f$f\f$
     * @param max_iters maximum number of iterations (warns user if exceeded)
     *
     * @return false if max iterations exceeded, true otherwise
     */
    template <typename CpProvider, typename HProvider>
    bool getTFromRhoE(
        const CpProvider& cp,
        const HProvider& h,
        const double rhoe,
        double& T,
        double* const p_work,
        const double alpha = 0.0,
        const double atol = 1.0e-12,
        const double rtol = 1.0e-12,
        const int max_iters = 100)
    {
        const int ns = m_thermo.nSpecies();
        const double rhoe_over_Ru = rhoe/RU;
        const double tol = rtol*std::abs(rhoe_over_Ru) + atol;

        double f, fp, dT;

        // Compute initial value of f
        h(T, p_work);
        f = alpha;
        for (int i = 0; i < ns; ++i)
            f += mp_X[i]*p_work[i];
        f = T*f - rhoe_over_Ru;

        int iter = 0;
        //cout << iter << " " << f << " " << T << endl;
        while (std::abs(f) > tol) {
            // Check for max iterations
            if (iter++ == max_iters) {
                using std::cerr;
                cerr << "Exceeded max iterations when computing temperature!\n";
                cerr << "res = " << f / rhoe_over_Ru << ", T = " << T << endl;
                return false;
            }

            // Compute df/dT
            cp(T, p_work);
            fp = alpha;
            for (int i = 0; i < ns; ++i)
                fp += mp_X[i]*p_work[i];

            // Update T
            dT = f/fp;
            while (dT >= T) dT *= 0.5; // prevent non-positive T
            T -= dT;

            // Recompute f
            h(T, p_work);
            f = alpha;
            for (int i = 0; i < ns; ++i)
                f += mp_X[i]*p_work[i];
            f = T*f - rhoe_over_Ru;
            //cout << iter << " " << f << " " << T << endl;
        }

        // Let the user know if we converged or not
        return true;
    }

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
    
}; // class StateModel

    } // namespace Thermodynamics
} // namespace Mutation

#endif // THERMO_STATE_MODEL_H

