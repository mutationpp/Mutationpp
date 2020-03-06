/**
 * @file MixtureOptions.h
 *
 * @brief Provides MixtureOptions class declaration.
 * @see Mutation::MixtureOptions
 */

/*
 * Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Mutation++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Mutation++.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef GENERAL_MIXTURE_OPTIONS_H
#define GENERAL_MIXTURE_OPTIONS_H

#include <string>
#include <vector>

#include "XMLite.h"
#include "Composition.h"

namespace Mutation {

/**
 * Keeps track of the set of mixture options in a convenient way.  Also,
 * enforces the default options so that a user will not have to set any options
 * unless they wish to change from the default behavior.
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
     * Copy constructor.
     */
    MixtureOptions(const MixtureOptions& options)
        : m_species_descriptor(options.m_species_descriptor),
          m_compositions(options.m_compositions),
          m_default_composition(options.m_default_composition),
          m_load_transport(options.m_load_transport),
          m_source(options.m_source),
          m_state_model(options.m_state_model),
          m_thermo_db(options.m_thermo_db),
          m_mechanism(options.m_mechanism),
          m_viscosity(options.m_viscosity),
          m_thermal_conductivity(options.m_thermal_conductivity),
          m_gsi_mechanism(options.m_gsi_mechanism)
    { }
    
    /**
	 * Constructs a new MixtureOptions object from a mixture file (without the
	 * .xml extension).  The file is first assumed to reside in the relative
	 * path from the local directory.  If it does not exist there, then
	 * MPP_DATA_DIRECTORY/mixtures is checked.
	 */
	MixtureOptions(const std::string& mixture);

	/**
	 * Constructs a new MixtureOptions object from a mixture file (without the
	 * .xml extension).  The file is first assumed to reside in the relative
	 * path from the local directory.  If it does not exist there, then
	 * MPP_DATA_DIRECTORY/mixtures is checked.
	 */
	MixtureOptions(const char* mixture);

	/**
	 * Constructs a new MixtureOptions from an XmlElement.
	 */
	MixtureOptions(Utilities::IO::XmlElement& element);

	/**
	 * Destructor.
	 */
	~MixtureOptions() {}

    /**
     * Assignment operator.
     */
    MixtureOptions& operator=(MixtureOptions options)
    {
        swap(*this, options);
        return *this;
    }
    
    /**
     * Loads the mixture options from a mixture input file.
     */
    void loadFromFile(const std::string& mixture);

    /**
     * Loads the mixture options from an XmlElement.
     */
    void loadFromXmlElement(Utilities::IO::XmlElement& element);

    /**
     * Gets the source of this mixture options (mixture file name).
     */
    const std::string& getSource() const {
		return m_source;
    }

    /**
     * Sets the options back to a default state.
     */
    void setDefaultOptions();
    
    /**
     * Gets the list of species names.
     */
    const std::string& getSpeciesDescriptor() const {
        return m_species_descriptor;
    }
    
    /**
     * Sets the list of species names.
     */
    void setSpeciesDescriptor(const std::string& descriptor) {
        m_species_descriptor = descriptor;
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
        m_state_model.assign(state_model);
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
     * Gets the Gas-Surface Interaction mechanism to use.
     */
    const std::string& getGSIMechanism() const {
        return m_gsi_mechanism;
    }

//    /**
//     * Sets the default mixture composition in elemental mole fractions.
//     */
//    void setDefaultComposition(
//        const std::vector<std::pair<std::string, double> >& composition) {
//        m_default_composition.assign(composition.begin(), composition.end());
//        m_has_default_composition = true;
//    }
//
//    /**
//     * Sets the default mole fraction for a single element.
//     */
//    void setDefaultComposition(const std::string& element, const double X) {
//        m_default_composition.push_back(std::make_pair(element, X));
//        m_has_default_composition = true;
//    }
//
//    /**
//     * Simple class that allows the user to use a simple syntax for setting the
//     * default element composition.
//     *
//     * @see MixtureOptions.setDefaultComposition()
//     */
//    class CompositionSetter
//    {
//    public:
//        CompositionSetter(
//            std::vector<std::pair<std::string, double> >& composition)
//            : m_composition(composition)
//        { }
//
//        CompositionSetter(const CompositionSetter& setter)
//            : m_composition(setter.m_composition)
//        { }
//
//        CompositionSetter& operator=(CompositionSetter setter)
//        {
//            std::swap(m_composition, setter.m_composition);
//            return *this;
//        }
//
//        CompositionSetter& operator () (
//            const std::string& element, const double X)
//        {
//            // Just replace the current composition for this element if it has
//            // already been set once
//            for (int i = 0; i < m_composition.size(); ++i) {
//                if (m_composition[i].first == element) {
//                    m_composition[i].second = X;
//                    return *this;
//                }
//            }
//
//            // Otherwise add this element to the default composition set
//            m_composition.push_back(std::make_pair(element, X));
//            return *this;
//        }
//
//    private:
//
//        std::vector<std::pair<std::string, double> >& m_composition;
//    };
//
//    /**
//     * Allows the user the set the default mole fractions using the simplified
//     * syntax.  For example:
//     *
//     * @code
//     * options.setDefaultComposition()
//     *     ("N",  0.79)
//     *     ("O",  0.21)
//     *     ("e-", 0.0);
//     * @endcode
//     */
//    CompositionSetter setDefaultComposition() {
//        m_default_composition.clear();
//        m_has_default_composition = true;
//        return CompositionSetter(m_default_composition);
//    }
//
    bool setDefaultComposition(const std::string& name);

    /**
     * Returns true if a default composition has been set either in a mixture
     * file or explicitly by the user.  Otherwise, using getDefaultComposition()
     * is meaningless.
     */
    bool hasDefaultComposition() const {
        return m_default_composition >= 0;
    }

    /**
     * Returns index of default composition in composition list if defined,
     * otherwise the index is < 0.
     */
    int getDefaultComposition() const {
        return m_default_composition;
    }

    bool addComposition(
        const Thermodynamics::Composition& c, bool make_default = false);
    
    const std::vector<Thermodynamics::Composition>& compositions() const {
        return m_compositions;
    }

    friend void swap(MixtureOptions&, MixtureOptions&);

private:

    void loadElementCompositions(const Utilities::IO::XmlElement& element);

private:

    std::string m_species_descriptor;
    
    std::vector<Thermodynamics::Composition> m_compositions;
    int m_default_composition;
    //std::vector<std::pair<std::string, double> > m_default_composition;
    //bool m_has_default_composition;

    bool m_load_transport;

	std::string m_source;
    std::string m_state_model;
    std::string m_thermo_db;
    std::string m_mechanism;
    std::string m_viscosity;
    std::string m_thermal_conductivity;
    std::string m_gsi_mechanism;

}; // class MixtureOptions

/**
 * Performs a swap on two MixtureOption objects.
 */
void swap(MixtureOptions& opt1, MixtureOptions& opt2);

} // namespace Mutation

#endif // GENERAL_MIXTURE_OPTIONS_H

