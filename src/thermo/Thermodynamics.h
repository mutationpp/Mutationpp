#ifndef THERMO_THERMODYNAMICS_H
#define THERMO_THERMODYNAMICS_H

#include <map>
#include <numeric>
#include <string>
#include <vector>

//#include "StateModel.h"
#include "Species.h"
#include "Numerics.h"
#include "Constants.h"
#include "ThermoDB.h"
#include "MultiPhaseEquilSolver.h"

namespace Mutation {
    namespace Thermodynamics {

//class GfcEquilSolver;
class MultiPhaseEquilSolver;
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
    
    // Mixed quantities
    Y_TO_XE,     ///< species mass fractions to elemental mole fractions
    X_TO_XE,     ///< species mole fractions to elemental mole fractions
};

/**
 * @todo Fix the gibbsOverRT funtion to take a pointer instead of vector.
 */
class Thermodynamics //: public StateModelUpdateHandler
{
public:
    
    /**
     * Constructs a Thermodynamics object given a vector of species names.  The
     * names must correspond to defined species in the thermo database of 
     * mutation++.
     */
    Thermodynamics(
        const std::string& species_descriptor,
        const std::string& database,
        const std::string& state_model);
    
    /**
     * Destructor.
     */
    ~Thermodynamics();
    
    /*void stateUpdated() {
        std::cout << "stateUpdated: Thermodynamics" << std::endl;
    }*/
    
    /**
     * Returns a pointer to the StateModel object owned by this thermodynamics
     * object.
     */
    StateModel* const state() const {
        return mp_state;
    }
    
    /**
     * Returns a pointer to the Equilibrium solver object owned by this
     * Thermodynamics object.
     */
    MultiPhaseEquilSolver* const  equilSolver() const {
        return mp_equil;
    }

    /**
     * Returns a pointer to the ThermoDB object owned by this Thermodynamics
     * object.
     */
    ThermoDB* thermoDB() const {
        return mp_thermodb;
    }

    /**
     * Returns the number of species considered in the mixture.
     */
    int nSpecies() const { 
        return mp_thermodb->species().size();
    }
    
    /**
     * Returns the number of heavy particles (non electrons) in the mixture.
     */
    int nHeavy() const {
        return m_natoms + m_nmolecules;
    }
    
    /**
     * Returns the number of atomic species in the mixture.
     */
    int nAtoms() const {
        return m_natoms;
    }
    
    /**
     * Returns the number of molecules in the mixture.
     */
    int nMolecules() const {
        return m_nmolecules;
    }
    
    /**
     * Returns the number of elements considered in the mixture.
     */
    int nElements() const {
        return mp_thermodb->elements().size();
    }
    
    /**
     * Returns the number of phases belonging to this mixture.
     */
    int nPhases() const;
    
    /**
     *  Returns number of gas species in the mixture.
     */
    int nGas() const { return m_ngas; }

    /**
     * Returns the number of condensed phase species in the mixture.
     */
    int nCondensed() const { return nSpecies()-m_ngas; }

    /**
     * Returns the number of energy equations associated with the mixture
     * StateModel.
     */
    int nEnergyEqns() const;

    /**
     * Returns the number of mass equations associated with the mixture
     * StateModel.
     */
    int nMassEqns() const;

    /**
     * Returns true if this mixture includes electrons, false otherwise.
     */
    bool hasElectrons() const {
        return m_has_electrons;
    }
    
    /**
     * Returns the species at index i.
     */
    const Species& species(int i) const {
        assert(i >= 0);
        assert(i < nSpecies());
        return mp_thermodb->species()[i];
    }
    
    /**
     * Returns the species with a given name.
     */
    const Species& species(std::string name) const {
        return species(speciesIndex(name));
    }
    
    const Element& element(int i) const {
        assert(i >= 0);
        assert(i < nElements());
        return mp_thermodb->elements()[i];
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
     * given name.  If no species exists with that name, then the index is < 0.
     */
    int speciesIndex(const std::string &name) const {
        std::map<std::string, int>::const_iterator iter = 
            m_species_indices.find(name);
        return (iter != m_species_indices.end() ? iter->second : -1);
    }
    
    /**
     * Returns the index in the element array assigned to the element with the 
     * given name.  If no element exists with that name, then the index is < 0.
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
        return species(index).name();
    }
    
    /**
     * Returns the charge of the species with the given index.
     */
    double speciesCharge(const int &index) const {
        return species(index).charge()*QE;
    }
    
    /**
     * Returns the name of the element with the given index in the element
     * array.
     */
    const std::string &elementName(const int &index) const {
        return element(index).name();
    }
    
    /**
     * Returns true if the thermodynamic data in this database is valid at the
     * given temperature for the given species index.
     */
    bool speciesThermoValidAtT(const size_t i, const double T) const;
    
    /**
     * Sets the state of the mixture using the StateModel belonging to the
     * mixture.  The input variables depend on the type of StateModel being
     * used.
     */
    void setState(
        const double* const p_v1, const double* const p_v2, const int vars = 0);
        
    /**
     * Computes the equilibrium composition of the mixture at the given fixed
     * temperature and pressure using the default elemental composition.
     */
    void equilibriumComposition(
        double T, double P, double* const p_X, MoleFracDef mdf = GLOBAL) const {
        equilibriumComposition(T, P, mp_default_composition, p_X, mdf);
    }

    /**
     * Computes the equilibrium composition of the mixture at the given fixed
     * temperature, pressure, and elemental moles/mole fractions.
     */
    void equilibriumComposition(
        double T, double P, const double* const p_Xe, double* const p_X,
        MoleFracDef mdf = GLOBAL) const;
    
    /**
     * Returns the element potentials calculated by the most recent call to 
     * equilibrate().
     */
    void elementPotentials(double* const p_lambda);
    
    /**
     * Returns the phase moles calculated by the most recent call to
     * equilibrate().
     */
    void phaseMoles(double* const p_moles);
    
    /**
     * Returns the number of continuation steps on the last call to equilibrate.
     */
    int nEquilibriumSteps() const;
    
    /**
     * Returns the total number of newton iterations on the last call to
     * equilibrate.
     */
    int nEquilibriumNewtons() const;
    
    /**
     * Adds a new set of linear constraints to the equilibrium solver.
     * @see clearEquilibriumConstraints()
     */
    void addEquilibriumConstraint(const double* const p_A);
    
    /**
     * Removes all of the equilibrium constraints.
     * @see addEquilibriumConstraint()
     */
    void clearEquilibriumContraints();
    
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
     * Fills temperature array with tempertures according to the used
     * StateModel.
     */
    void getTemperatures(double* const p_T) const;

    /**
     * Fills energy per mass array with energies according to the used
     * StateModel (total + internal for each species).
     */
    void getEnergiesMass(double* const p_e) const;
    
    /**
     * Fills enthalpy per mass array with enthalpy according to the used
     * StateModel (total + internal for each species).
     */

    void getEnthalpiesMass(double* const p_h) const;
    
    /**
     * Fills the constant pressure specific heat according to the used
     * StateModel
     */
    
    void getCpsMass(double* const p_cp) const;

	/**
     * Fills the constant volume specific heat according to the used
     * StateModel
     */
    
    void getCvsMass(double* const p_cv) const;

    /**
     * Fills the tag of the modes according toi the used StateModel
     */
    void getTagModes(int* const p_tag) const;
    
    /**
     * Returns the current mixture static pressure in Pa.
     */
    double P() const;
    
    /**
     * Returns the current species mole fractions.
     */
    const double* const X() const;
    
    /**
     * Returns the current species mass fractions.
     */
    const double* const Y() const;
    
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
        return element(index).atomicMass();
    }
    
    /**
     * Returns the species vector.
     */
    const std::vector<Species>& species() const {
        return mp_thermodb->species();
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
    
    void speciesCpOverR(double T, double* const p_cp) const;
    
    void speciesCpOverR(
        double Th, double Te, double Tr, double Tv, double Tel,
        double *const p_cp, double *const p_cpt, double *const p_cpr,
        double *const p_cpv, double *const p_cpel) const;
    
    /**
     * Returns the unitless vector of species specific heats at constant volume
     * \f$ C_{v,i} / R_u \f$.
     */
    void speciesCvOverR(
        double Th, double Te, double Tr, double Tv, double Tel,
        double *const p_cv, double *const p_cvt, double *const p_cvr,
        double *const p_cvv, double *const p_cvel) const;
    
    void speciesCvOverR(
        double *const p_cv, double *const p_cvt, double *const p_cvr,
        double *const p_cvv, double *const p_cvel) const
    {
        speciesCvOverR(
            T(), Te(), Tr(), Tv(), Tel(), p_cv, p_cvt, p_cvr, p_cvv, p_cvel);        
    }
    
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
     * Returns the current equilibrium mixture averaged specific heat at
     * constant pressure in J/mol-K.  This method assumes that the current state
     * of the mixture was already set using the equilibrate() method.
     */
    double mixtureEquilibriumCpMole();
    
    /**
     * Returns the current equilibrium mixture averaged specific heat at
     * constant pressure in J/kg-K.  This method assumes that the current state
     * of the mixture was already set using the equilibrate() method.
     */
    double mixtureEquilibriumCpMass();
    
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
    //double mixtureEquilibriumCvMole(
    //    double T, double P, const double* const Xeq) const;
    
    /**
     * Returns the current equilibrium mixture averaged specific heat at
     * constant volume in J/mol-K.  This method assumes that the current state of
     * the mixture was already set using the equilibrate() method.
     */
    //double mixtureEquilibriumCvMole() const;
    
    /**
     * Returns the equilibrium mixture averaged specific heat at constant
     * volume in J/kg-K.
     */
    //double mixtureEquilibriumCvMass(
    //    double T, double P, const double* const Xeq) const;
    
    /**
     * Returns the current equilibrium mixture averaged specific heat at
     * constant volume in J/kg-K.  This method assumes that the current state of
     * the mixture was already set using the equilibrate() method.
     */
    double mixtureEquilibriumCvMass();
    
    /**
     * Returns the ratio of mixture averaged specific heats which is a unitless
     * quantity.
     */
    double mixtureFrozenGamma() const;
    
    /**
     * Returns the equilibrium ratio of mixture averaged specific heats which is
     * a unitless quantity.
     *
     * T   - equilibrium temperature in K
     * P   - equilibrium pressure in Pa
     * Xeq - equilibrium species mole fractions
     */
    double mixtureEquilibriumGamma(
        double T, double P, const double* const Xeq);
    
    /**
     * Returns the current equilibrium ratio of mixture averaged specific heats
     * which is a unitless quantity.
     */
    double mixtureEquilibriumGamma();
    
    /**
     * Returns the derivative of mole fraction w.t.r temperature for the given
     * equilibrium mixture.
     * dxdt - on return, the dX_i/dT derivatives
     */
    void dXidT(double* const dxdt) const;
    
    /**
     * Returns the species derivatives of mole fraction w.r.t. pressure for the
     * given equilibrium mixture.  Note that it is assumed the state model is an
     * equilibrium one.
     */
    void dXidP(double* const dxdp) const;

    /**
     * Returns the density derivative with respect to pressure for the current
     * equilibrium state.
     */
    double dRhodP();
    
    /**
     * Returns the unitless vector of species enthalpies \f$ H_i / R_u T \f$.
     */
    void speciesHOverRT(double T, double* const h) const; 
    
    /**
     * Computes the unitless species enthalpies and can optionally fill vectors
     * for each energy mode.  
     */
    void speciesHOverRT(
        double* const h, double* const ht = NULL, 
        double* const hr = NULL, double* const hv = NULL,
        double* const hel = NULL, double* const hf = NULL) const;

     /**
     * Computes the unitless species enthalpies and can optionally fill vectors
     * for each energy mode by explicitly passing each individual temperature .  
     */
    void speciesHOverRT(
        double T, double Te, double Tr, double Tv, double Tel,
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
     * Returns the mixture energy in J/mol.
     */
    double mixtureEnergyMole() const {
        return mixtureHMole() - P() / density() * mixtureMw();
    }
    
    /**
     * Returns the mixture energy in J/kg.
     */
    double mixtureEnergyMass() const {
        return mixtureHMass() - P() / density();
    }
    
    /**
     * Returns the frozen sound speed of the mixture in m/s.
     */
    double frozenSoundSpeed() const {
        return std::sqrt(mixtureFrozenGamma() * P() / density());
    }
    
    /**
     * Returns the current equilibrium sound speed of the mixture in m/s.
     */
    double equilibriumSoundSpeed();// {
    //    return std::sqrt(mixtureEquilibriumGamma() / dRhodP());
    //}
    
    /**
     * Returns the unitless vector of species entropies \f$ S_i / R_u \f$.
     */
    void speciesSOverR(double* const p_s) const;
    
    /**
     * Returns the mixture averaged entropy in J/mol-K including the entropy of
     * mixing.
     */
    double mixtureSMole() const;
    
    /**
     * Returns the mixture averaged entropy in J/kg-K including the entropy of
     * mixing.
     */
    double mixtureSMass() const;
    
    /**
     * Returns the unitless vector of species Gibbs free energies
     * \f$ G_i / R_u T = H_i / R_u T - S_i / R_u \f$.
     */
    void speciesGOverRT(double* const p_g) const;
    void speciesGOverRT(double T, double P, double* const p_g) const;
    void speciesSTGOverRT(double T, double* const p_g) const;
    
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
    
    /**
     * Solves the surface mass balance equation for a pure Carbon char and
     * returns the non-dimensional ablation mass flux and wall enthalpy.
     */
    void surfaceMassBalance(
        const double *const p_Yke, const double *const p_Ykg, const double T, 
        const double P, const double Bg, double &Bc, double &hw,
        double *const p_Xs = NULL);
    
private:
    
    /**
     * Loads the elements and species given by the species names vector from the
     * thermo database.
     */
    void loadSpeciesFromList(const std::vector<std::string> &species_names);
    
    /**
     * Adds a single species to the thermodynamic list.
     */
    void addSpecies(const Species& species);

private:
  
    ThermoDB* mp_thermodb;
    MultiPhaseEquilSolver* mp_equil;
    StateModel* mp_state;
    
    std::map<std::string, int> m_species_indices;
    std::map<std::string, int> m_element_indices;
    
    Numerics::RealMatrix m_element_matrix;
    Numerics::RealVector m_species_mw;
    
    double* mp_work1;
    double* mp_work2;
    double* mp_wrkcp;
    double* mp_y;
    double* mp_default_composition;
    
    bool m_has_electrons;
    int  m_natoms;
    int  m_nmolecules;
    int  m_ngas;
    
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

//template <>
//inline void Thermodynamics::convert<Y_TO_XE>(
//    const double *const a, double *const b) const {
//
//}
/// @endcond

    } // namespace Thermodynamics
} // namespace Mutation

#endif // THERMO_THERMODYNAMICS_H

