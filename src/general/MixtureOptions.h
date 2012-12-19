#ifndef GENERAL_MIXTURE_OPTIONS_H
#define GENERAL_MIXTURE_OPTIONS_H

#include <string>
#include <vector>

/**
 * Keeps track of the set of mixture options in a convient way.  Also, enforces
 * the default options so that a user will not have to set any options unless
 * they wish to change from the default.
 */
class MixtureOptions
{
public:

    /**
     * Constructs a new MixtureOptions object with default options and empty
     * species list.
     */
    MixtureOptions();
    
    /**
     * Constructs a new MixtureOptions object from a mixture file.
     */
    MixtureOptions(const std::string& mixture);
    
    /**
     * Constructs a new MixtureOptions object from a mixture file.
     */
    MixtureOptions(const char* mixture);
    
    /**
     * Destructor.
     */
    ~MixtureOptions() {}
    
    /**
     * Sets the options back to a default state.
     */
    void setDefaultOptions();
    
    /**
     * Gets the list of species names.
     */
    const std::vector<std::string>& getSpeciesNames() const {
        return m_species_names;
    }
    
    /**
     * Adds a species name to the list.
     */
    void addSpeciesName(const std::string& species) {
        m_species_names.push_back(species);
    }
    
    /**
     * Sets the list of species names.
     */
    void setSpeciesNames(const std::vector<std::string>& species) {
        m_species_names = species;
    }
    
    /**
     * Gets the mixture state model to be used.
     */
    const std::string& getStateModel() const {
        return m_state_model;
    }
    
    /**
     * Sets the mixture state model to be used.
     */
    void setStateModel(const std::string& state_model) {
        m_state_model = state_model;
    }
    
    /**
     * Gets the thermodynamic database type to use.
     */
    const std::string& getThermodynamicDatabase() const {
        return m_thermo_db;
    }
    
    /**
     * Sets which thermodynamic database is to be used.
     */
    void setThermodynamicDatabase(const std::string& thermo_db) {
        m_thermo_db = thermo_db;
    }
    
    /**
     * Gets the name of the reaction mechanism to use.
     */
    const std::string& getMechanism() const {
        return m_mechanism;
    }
    
    /**
     * Sets the name of the reaction mechanism.
     */
    void setMechanism(const std::string& mechanism) {
        m_mechanism = mechanism;
    }
    
    /**
     * Gets the viscosity algorithm to use.
     */
    const std::string& getViscosityAlgorithm() const {
        return m_viscosity;
    }
    
    /**
     * Sets which viscosity algorithm to use.
     */
    void setViscosityAlgorithm(const std::string& viscosity) {
        m_viscosity = viscosity;
    }
    
    /**
     * Gets the thermal conductivity algorithm to use.
     */
    const std::string& getThermalConductivityAlgorithm() const {
        return m_thermal_conductivity;
    }
    
    /**
     * Sets the thermal conductivity algorithm to use.
     */
    void setThermalConductivityAlgorithm(
        const std::string& thermal_conductivity) 
    {
        m_thermal_conductivity = thermal_conductivity;
    }
    
    /**
     * Sets the default mixture composition in elemental mole fractions.
     */
    void setDefaultComposition(
        const std::vector<std::pair<std::string, double> >& composition) {
        m_default_composition.assign(composition.begin(), composition.end());
        m_has_default_composition = true;
    }
    
    /**
     * Sets the default mole fraction for a single element.
     */
    void setDefaultComposition(const std::string& element, const double X) {
        m_composition_setter(element, X);
        m_has_default_composition = true;
    }
    
    /**
     * Simple class that allows the user to use a simple syntax for setting the
     * default element composition.
     * 
     * @see MixtureOptions.setDefaultComposition()
     */
    class CompositionSetter
    {
    public:
        CompositionSetter(
            std::vector<std::pair<std::string, double> >& composition)
            : m_composition(composition)
        { }
        
        CompositionSetter& operator () (
            const std::string& element, const double X)
        {
            m_composition.push_back(std::make_pair(element, X));
            return *this;
        }
        
    private:
        
        std::vector<std::pair<std::string, double> >& m_composition;
    };
    
    /**
     * Allows the user the set the default mole fractions using the simplified
     * syntax.  For example:
     *
     * @code
     * options.setDefaultComposition()
     *     ("N",  0.79)
     *     ("O",  0.21)
     *     ("e-", 0.0);
     * @endcode
     */
    CompositionSetter& setDefaultComposition() {
        m_default_composition.clear();
        m_has_default_composition = true;
        return m_composition_setter;
    }
    
    /**
     * Returns the default composition as a vector of (element name, elemental
     * fraction) pairs.
     */
    const std::vector<std::pair<std::string, double> >& getDefaultComposition()
        const
    {
        return m_default_composition;
    }
    
    /**
     * Returns true if a default composition has been set either in a mixture
     * file or explicitly by the user.  Otherwise, using getDefaultComposition()
     * is meaningless.
     */
    bool hasDefaultComposition() const {
        return m_has_default_composition;
    }

private:

    /**
     * Loads the mixture options from a mixture input file.
     */
    void loadFromFile(const std::string& mixture);

private:

    std::vector<std::string> m_species_names;
    
    std::vector<std::pair<std::string, double> > m_default_composition;
    CompositionSetter m_composition_setter;
    bool m_has_default_composition;

    std::string m_state_model;
    std::string m_thermo_db;
    std::string m_mechanism;
    std::string m_viscosity;
    std::string m_thermal_conductivity;

}; // class MixtureOptions

#endif // GENERAL_MIXTURE_OPTIONS_H

