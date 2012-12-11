#ifndef THERMO_THERMODYNAMICS_H
#define THERMO_THERMODYNAMICS_H

#include <map>
#include <numeric>
#include <string>
#include <vector>

#include "Species.h"
#include "Numerics.h"

class ThermoDB;
class StateModel;

/**
 * Possible conversion methods that can be used with the 
 * Thermodynamics::convert() function.
 *
 * @see Thermodynamics::convert()
 */
enum ConversionType {
    // Species quantities
    RHO_TO_CONC, ///< species density in kg/m^3 to molar density in mol/m^3
    RHO_TO_Y,    ///< species density in kg/m^3 to mass fraction
    RHO_TO_X,    ///< species density in kg/m^3 to mole fraction
    CONC_TO_RHO, ///< species molar density in mol/m^3 to density in kg/m^3
    CONC_TO_Y,   ///< species molar density in mol/m^3 to mass fraction
    CONC_TO_X,   ///< species molar density in mol/m^3 to mole fraction
    X_TO_Y,      ///< species mole fraction to mass fraction
    Y_TO_X,      ///< species mass fraction to mole fraction
    
    // Elemental quantities
    XE_TO_YE,    ///< elemental mole fraction to mass fraction
    YE_TO_XE,    ///< elemental mass fraction to mole fraction
};

// Forward declaration of GfcEquilSolver
class GfcEquilSolver;

/**
 * @todo Fix the gibbsOverRT funtion to take a pointer instead of vector.
 */
class Thermodynamics
{
public:
    
    /**
     * Constructs a Thermodynamics object given a vector of species names.  The
     * names must correspond to defined species in the thermo database of 
     * mutation++.
     */
    Thermodynamics(
        const std::vector<std::string> &species_names, 
        const std::string& database, const std::string& state_model);
    
    /**
     * Destructor.
     */
    ~Thermodynamics();
    
    /**
     * Returns the number of species considered in the mixture.
     */
    int nSpecies() const { 
        return m_species.size();
    }
    
    /**
     * Returns the number of elements considered in the mixture.
     */
    int nElements() const {
        return m_elements.size();
    }
    
    /**
     * Returns true if this mixture includes electrons, false otherwise.
     */
    bool hasElectrons() const {
        return m_has_electrons;
    }
    
    /**
     * Sets the default elemental composition of the mixture from a given
     * vector of (element name, mole fraction) pairs.  Note that an error will
     * result if each element is not included exactly once in the list or if an
     * element that does not exist in the mixture is included.
     */
    void setDefaultComposition(
        const std::vector<std::pair<std::string, double> >& composition);
    
    /**
     * Returns the default elemental fraction for the element with a given
     * index.  Note that this index should be in 0 <= i < nElements().
     */
    double getDefaultComposition(const int index) {
        return mp_default_composition[index];
    }
    
    /**
     * Returns the index in the species array assigned to the species with the
     * given name.
     */
    int speciesIndex(const std::string &name) const {
        std::map<std::string, int>::const_iterator iter = 
            m_species_indices.find(name);
        return (iter != m_species_indices.end() ? iter->second : -1);
    }
    
    /**
     * Returns the index in the element array assigned to the element with the 
     * given name.
     */
    int elementIndex(const std::string &name) const {
        std::map<std::string, int>::const_iterator iter = 
            m_element_indices.find(name);
        return (iter != m_element_indices.end() ? iter->second : -1);
    }
    
    /**
     * Returns the name of the species with the given index in the species 
     * array.
     */
    const std::string &speciesName(const int &index) const {
        return m_species[index].name();
    }
    
    /**
     * Returns the name of the element with the given index in the element
     * array.
     */
    const std::string &elementName(const int &index) const {
        return m_elements[index].name();
    }
    
    /**
     * Sets the current state of the mixture using temperatures, pressure, and 
     * species mole fractions.
     */
    void setStateTPX(
        const double* const T, const double* const P, const double* const X);
    
    /**
     * Sets the current state of the mixture using temperatures pressures, and 
     * species mass fractions.
     */
    void setStateTPY(
        const double* const T, const double* const P, const double* const Y);
        
    /**
     * Computes the equilibrium mole fractions of the mixture given the
     * elemental composition.
     */
    void equilibrate(
        double T, double P, const double* const p_c, double* const p_X,
        bool set_state = true)
        const;
    
    /**
     * Computes the equilibrium mole fractions of the mixture using the default
     * elemental composition.
     */
    void equilibrate(double T, double P) const;
    
    /**
     * Returns the mixture translational temperature.
     */
    double T() const;
    
    /**
     * Returns the mixture rotational temperature.
     */
    double Tr() const;
    
    /**
     * Returns the mixture vibrational temperature.
     */
    double Tv() const;
    
    /**
     * Returns the electron temperature.
     */
    double Te() const;
    
    /**
     * Returns the mixture electronic temperature.
     */
    double Tel() const;
    
    /**
     * Returns the current mixture static pressure in Pa.
     */
    double P() const;
    
    /**
     * Returns the current species mass fractions.
     */
    const double* const X() const;
    
    /**
     * Returns the standard state temperature for the thermodynamic database 
     * begin used in K.
     */
    double standardStateT() const;
    
    /**
     * Returns the standard state pressure for the thermodynamic database being
     * used in Pa.
     */
    double standardStateP() const;
    
    /**
     * Returns the average molecular weight of the mixture given the species
     * mass fractions.
     */
    double mixtureMwMass(const double *const Y) const {
        return 1.0 / ((Numerics::asVector(Y, nSpecies()) / m_species_mw).sum());
    }
    
    /**
     * Returns the average molecular weight of the mixture given the species
     * mole fractions.
     */
    double mixtureMwMole(const double *const X) const {
        return (Numerics::asVector(X, nSpecies()) * m_species_mw).sum();
    }
    
    /**
     * Returns the average molecular weight of the mixture.
     */
    double mixtureMw() const;
    
    /**
     * Returns the molecular weight in kg/mol of the species with the given 
     * index in the species array. 
     */
    double speciesMw(const int index) const { 
        return m_species_mw(index);
    }
    
    /**
     * Returns the atomic weight in kg/mol of the element with the given index
     * in the element array.
     */
    double atomicMass(const int index) const {
        return m_elements[index].atomicMass();
    }
    
    /**
     * Returns the species vector.
     */
    const std::vector<Species> &species() const {
        return m_species;
    }
    
    /**
     * Returns the element matrix \f$ E_{i,j} \f$ in which the (i,j)
     * component is equal to the number of atoms of the jth element in one
     * molecule of the ith species.
     */
    const Numerics::RealMatrix &elementMatrix() const {
        return m_element_matrix;
    }
    
    /**
     * Returns the number density of the mixture given the mixture temperature
     * and pressure.
     */
    double numberDensity(const double T, const double P) const;
    
    /**
     * Returns the current number density of the mixture.
     */
    double numberDensity() const;
    
    /**
     * Returns the pressure of the mixture in Pa given the mixture temperature 
     * and density and species mass fractions.
     */
    double pressure(
        const double T, const double rho, const double *const Y) const;
    
    /**
     * Returns the density of the mixture in kg/m^3 given the mixture 
     * temperature and pressure and species mole fractions.
     */
    double density(
        const double T, const double P, const double *const X) const;
        
    /**
     * Returns the current mixture density in kg/m^3.
     */
    double density() const;

    /**
     * Returns the unitless vector of species specific heats at constant
     * pressure \f$ C_{p,i} / R_u \f$.
     */
    void speciesCpOverR(double *const p_cp) const;
    
    /**
     * Returns the unitless vector of species specific heats at constant volume
     * \f$ C_{v,i} / R_u \f$.
     */
    void speciesCvOverR(
        double *const p_cv, double *const p_cvt, double *const p_cvr,
        double *const p_cvv, double *const p_cvel) const;
    
    /**
     * Returns the frozen mixture averaged specific heat at constant pressure in
     * J/mol-K.
     */
    double mixtureFrozenCpMole() const;
    
    /**
     * Returns the frozen mixture averaged specific heat at constant pressure in 
     * J/kg-K.
     */
    double mixtureFrozenCpMass() const;
    
    /**
     * Returns the equilibrium mixture averaged specific heat at constant
     * pressure in J/mol-K given the temperature, pressure and equilibrium mole
     * fractions.
     */
    double mixtureEquilibriumCpMole(
        double T, double P, const double* const Xeq) const;
    
    /**
     * Returns the current equilibrium mixture averaged specific heat at
     * constant pressure in J/mol-K.  This method assumes that the current state
     * of the mixture was already set using the equilibrate() method.
     */
    double mixtureEquilibriumCpMole() const;
    
    /**
     * Returns the equilibrium mixture averaged specific heat at constant
     * pressure in J/kg-K given the temperature, pressure and equilibrium mole
     * fractions.
     */
    double mixtureEquilibriumCpMass(
        double T, double P, const double* const Xeq) const;
    
    /**
     * Returns the current equilibrium mixture averaged specific heat at
     * constant pressure in J/kg-K.  This method assumes that the current state
     * of the mixture was already set using the equilibrate() method.
     */
    double mixtureEquilibriumCpMass() const;
    
    /**
     * Returns the mixture averaged specific heat at constant volume in J/mol-K.
     */
    double mixtureFrozenCvMole() const;
    
    /**
     * Returns the mixture averaged specific heat at constant volume in J/kg-K.
     */
    double mixtureFrozenCvMass() const;
    
    /**
     * Returns the equilibrium mixture averaged specific heat at constant
     * volume in J/mol-K.
     */
    double mixtureEquilibriumCvMole(
        double T, double P, const double* const Xeq) const;
    
    /**
     * Returns the current equilibrium mixture averaged specific heat at
     * constant volume in J/mol-K.  This method assumes that the current state of
     * the mixture was already set using the equilibrate() method.
     */
    double mixtureEquilibriumCvMole() const;
    
    /**
     * Returns the equilibrium mixture averaged specific heat at constant
     * volume in J/kg-K.
     */
    double mixtureEquilibriumCvMass(
        double T, double P, const double* const Xeq) const;
    
    /**
     * Returns the current equilibrium mixture averaged specific heat at
     * constant volume in J/kg-K.  This method assumes that the current state of
     * the mixture was already set using the equilibrate() method.
     */
    double mixtureEquilibriumCvMass() const;
    
    /**
     * Returns the ratio of mixture averaged specific heats which is a unitless
     * quantity.
     */
    double mixtureFrozenGamma() const;
    
    /**
     * Returns the equilibrium ratio of mixture averaged specific heats which is
     * a unitless quantity.
     */
    double mixtureEquilibriumGamma(
        double T, double P, const double* const Xeq) const;
    
    /**
     * Returns the current equilibrium ratio of mixture averaged specific heats
     * which is a unitless quantity.
     */
    double mixtureEquilibriumGamma() const;
    
    /**
     * Returns the unitless vector of species enthalpies \f$ H_i / R_u T \f$.
     */
    void speciesHOverRT(
        double* const h, double* const ht = NULL, 
        double* const hr = NULL, double* const hv = NULL,
        double* const hel = NULL, double* const hf = NULL) const;
    
    /**
     * Returns the mixture averaged enthalpy in J/mol.
     */
    double mixtureHMole() const;
    
    /**
     * Returns the mixture averaged enthalpy in J/kg.
     */
    double mixtureHMass() const;
    
    /**
     * Returns the unitless vector of species entropies \f$ S_i / R_u \f$.
     */
    void speciesSOverR(double* const p_s) const;
    
    /**
     * Returns the mixture averaged entropy in J/mol-K.
     */
    double mixtureSMole() const;
    
    /**
     * Returns the mixture averaged entropy in J/kg-K.
     */
    double mixtureSMass() const;
    
    /**
     * Returns the unitless vector of species Gibbs free energies
     * \f$ G_i / R_u T = H_i / R_u T - S_i / R_u \f$.
     */
    void speciesGOverRT(double* const p_g) const;
    void speciesGOverRT(double T, double P, double* const p_g) const;
    
    /**
     * Returns the number of moles of each element in a mixture with a given
     * set of species moles.
     */
    void elementMoles(
        const double *const species_N, double *const element_N) const;
    
    /**
     * Returns the element fractions given the species mole fractions.
     */
    void elementFractions(
        const double* const Xs, double* const Xe) const;
    
    /**
     * Converts species properties from one value to another such as mole to 
     * mass fractions.  The different conversion methods implemented are
     * enumerated in ConversionType.
     *
     * @see enum ConversionType
     */
    template <ConversionType>
    inline void convert(const double *const a, double *const b) const { }
    
private:
    
    /**
     * Loads the elements and species given by the species names vector from the
     * thermo database.
     */
    void loadSpeciesFromList(const std::vector<std::string> &species_names);

private:
  
    ThermoDB* mp_thermodb;
    GfcEquilSolver* mp_equil;
    StateModel* mp_state;
    
    std::vector<Species> m_species;    
    std::map<std::string, int> m_species_indices;
    
    std::vector<Element> m_elements;
    std::map<std::string, int> m_element_indices;
    
    Numerics::RealMatrix m_element_matrix;
    Numerics::RealVector m_species_mw;
    
    double* mp_work1;
    double* mp_work2;
    double* mp_default_composition;
    
    bool m_has_electrons;
    
}; // class Thermodynamics

/// @cond CONVERT_SPECIALIZATIONS
template <>
inline void Thermodynamics::convert<RHO_TO_CONC>(
    const double *const a, double *const b) const {
    for (int i = 0; i < nSpecies(); ++i)
        b[i] = a[i] / speciesMw(i);
}

template <>
inline void Thermodynamics::convert<CONC_TO_RHO>(
    const double *const a, double *const b) const {
    for (int i = 0; i < nSpecies(); ++i)
        b[i] = a[i] * speciesMw(i);
}

template <>
inline void Thermodynamics::convert<RHO_TO_Y>(
    const double *const a, double *const b) const {
    const double sum = std::accumulate(a, a + nSpecies(), 0.0);
    for (int i = 0; i < nSpecies(); ++i)
        b[i] = a[i] / sum;
}

template <>
inline void Thermodynamics::convert<CONC_TO_X>(
    const double *const a, double *const b) const {
    convert<RHO_TO_Y>(a, b);
}

template <>
inline void Thermodynamics::convert<RHO_TO_X>(
    const double *const a, double *const b) const {
    convert<RHO_TO_CONC>(a, b);
    convert<CONC_TO_X>(b, b);        
}

template <>
inline void Thermodynamics::convert<CONC_TO_Y>(
    const double *const a, double *const b) const {
    convert<CONC_TO_RHO>(a, b);
    convert<RHO_TO_Y>(b, b);
}

template <>
inline void Thermodynamics::convert<X_TO_Y>(
    const double *const a, double *const b) const {
    convert<CONC_TO_RHO>(a, b);
    convert<RHO_TO_Y>(b, b);
}

template <>
inline void Thermodynamics::convert<Y_TO_X>(
    const double *const a, double *const b) const {
    convert<RHO_TO_CONC>(a, b);
    convert<CONC_TO_X>(b, b);
}

template <>
inline void Thermodynamics::convert<XE_TO_YE>(
    const double *const a, double *const b) const {
    for (int i = 0; i < nElements(); ++i)
        b[i] = a[i] * atomicMass(i);
    const double sum = std::accumulate(b, b + nElements(), 0.0);
    for (int i = 0; i < nElements(); ++i)
        b[i] /= sum;
}

template <>
inline void Thermodynamics::convert<YE_TO_XE>(
    const double *const a, double *const b) const {
    for (int i = 0; i < nElements(); ++i)
        b[i] = a[i] / atomicMass(i);
    const double sum = std::accumulate(b, b + nElements(), 0.0);
    for (int i = 0; i < nElements(); ++i)
        b[i] /= sum;
}
/// @endcond

#endif // THERMO_THERMODYNAMICS_H
