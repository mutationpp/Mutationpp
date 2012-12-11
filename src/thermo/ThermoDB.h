#ifndef THERMO_THERMODB_H
#define THERMO_THERMODB_H

#include <vector>

// Forward declare the species class.
class Species;

/**
 * Abstract base class for all thermodynamic databases.  This provides the 
 * construct for self registering thermodynamic database types so that it is 
 * possible to dynamically select which database a user would like to use.  All
 * thermodynamic database classes must provide the functionality listed below as
 * well as a constructor which takes only a std::vector<Species>& as arguments.
 * Note that the cv() is implemented using the relation \f$Cv/Ru = Cp/Ru - 1\f$,
 * therefore it does not need to be implemented for concrete types of this
 * class. 
 */
class ThermoDB
{
public:
    
    /**
     * Type of arguments required in constructor for all classes derived from
     * ThermoDB.
     */
    typedef const std::vector<Species>& ARGS;
    
    /**
     * Constructor for all ThermoDB types takes the list of species for which
     * the the thermodynamic database should be created.
     */
    ThermoDB(ARGS species);
    
    /**
     * Desctructor.
     */
    virtual ~ThermoDB() {};
    
    /**
     * Returns the standard state temperature in K employeed in this database.
     */
    virtual double standardTemperature() const = 0;
    
    /**
     * Returns the standard state pressure in Pa employeed in this database.
     */
    virtual double standardPressure() const = 0;
    
    /**
     * Computes the unitless species specific heats at constant pressure
     * \f$ C_{p,i}/R_U\f$ of each species in thermal nonequilibrium.
     *
     * @param Th  - heavy particle translational temperatures
     * @param Te  - free electron temperature
     * @param Tr  - mixture rotational temperature
     * @param Tv  - mixture vibrational temperature
     * @param Tel - mixture electronic temperature
     * @param cp   - if not NULL, the array of species non-dimensional cp
     * @param cpt  - if not NULL, the array of species translational cp
     * @param cpr  - if not NULL, the array of species rotational cp
     * @param cpv  - if not NULL, the array of species vibrational cp
     * @param cpel - if not NULL, the array of species electronic cp
     */
    virtual void cp(
        double Th, double Te, double Tr, double Tv, double Tel, 
        double* const cp, double* const cpt, double* const cpr,
        double* const cpv, double* const cpel) = 0;
    
    /**
     * Computes the unitless species specific heats at constant volume
     * \f$ C_{v,i} / R_U\f$ of each species in thermal nonequilibrium.
     *
     * @param Th  - heavy particle translational temperatures
     * @param Te  - free electron temperature
     * @param Tr  - mixture rotational temperature
     * @param Tv  - mixture vibrational temperature
     * @param Tel - mixture electronic temperature
     * @param cv   - if not NULL, the array of species non-dimensional cv
     * @param cvt  - if not NULL, the array of species translational cv
     * @param cvr  - if not NULL, the array of species rotational cv
     * @param cvv  - if not NULL, the array of species vibrational cv
     * @param cvel - if not NULL, the array of species electronic cv
     */
    virtual void cv(
        double Th, double Te, double Tr, double Tv, double Tel, 
        double* const cv, double* const cvt, 
        double* const cvr, double* const cvv, 
        double* const cvel);
    
    /**
     * Computes the unitless species enthalpy \f$ h_i/R_U T_h\f$ of each
     * species in thermal nonequilibrium, which is non-dimensionalized by the 
     * heavy particle translational temperature. 
     *
     * @param Th  - heavy particle translational temperature
     * @param Te  - free electron temperature
     * @param Tr  - mixture rotational temperature
     * @param Tv  - mixture vibrational temperature
     * @param Tel - mixture electronic temperature
     * @param h   - on return, the array of species non-dimensional enthalpies
     * @param ht  - if not NULL, the array of species translational enthalpies
     * @param hr  - if not NULL, the array of species rotational enthalpies
     * @param hv  - if not NULL, the array of species vibrational enthalpies
     * @param hel - if not NULL, the array of species electronic enthalpies
     * @param hf  - if not NULL, the array of the species formation enthalpies
     */
    virtual void enthalpy(
        double Th, double Te, double Tr, double Tv, double Tel, double* const h,
        double* const ht, double* const hr, double* const hv, double* const hel, 
        double* const hf) = 0;
    
    /**
     * Computes the unitless species entropy \f$s_i/R_u\f$ allowing for thermal 
     * nonequilibrium.
     *
     * @param Th  - heavy particle translational temperature
     * @param Te  - free electron temperature
     * @param Tr  - mixture rotational temperature
     * @param Tv  - mixture vibrational temperature
     * @param Tel - mixture electronic temperature
     * @param P   - mixture static pressure
     * @param s   - on return, the array of species non-dimensional entropies
     * @param st  - if not NULL, the array of species translational entropies
     * @param sr  - if not NULL, the array of species rotational entropies
     * @param sv  - if not NULL, the array of species vibrational entropies
     * @param sel - if not NULL, the array of species electronic entropies
     */    
    virtual void entropy(
        double Th, double Te, double Tr, double Tv, double Tel, double P,
        double* const s, double* const st, double* const sr, double* const sv, 
        double* const sel) = 0;
    
    /**
     * Computes the unitless Gibbs free energy of each species i, 
     * \f$G_i / R_u T_h\f$ where \f$G_i\f$ is non-dimensionalized by the heavy
     * particle translational temperature.
     *
     * @param Th  - heavy particle translational temperature
     * @param Te  - free electron temperature
     * @param Tr  - mixture rotational temperature
     * @param Tv  - mixture vibrational temperature
     * @param Tel - mixture electronic temperature
     * @param P   - mixture static pressure
     * @param g   - on return, the array of species non-dimensional energies
     * @param gt  - if not NULL, the array of species translational energies
     * @param gr  - if not NULL, the array of species rotational energies
     * @param gv  - if not NULL, the array of species vibrational energies
     * @param gel - if not NULL, the array of species electronic energies
     */
    virtual void gibbs(
        double Th, double Te, double Tr, double Tv, double Tel, double P,
        double* const g, double* const gt, double* const gr, double* const gv,
        double* const gel) = 0;
    
protected:

    bool m_ns;
    
}; // class ThermoDB

#endif // THERMO_THERMODB_H

