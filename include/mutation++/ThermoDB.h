/**
 * @file ThermoDB.h
 *
 * @brief Declaration of the ThermoDB class.
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

#ifndef THERMO_THERMODB_H
#define THERMO_THERMODB_H

#include "SpeciesListDescriptor.h"
#include "Species.h"

#include <vector>
#include <list>

namespace Mutation {
    namespace Thermodynamics {

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
    typedef int ARGS;
    
    /// Returns name of this type.
    static std::string typeName() { return "ThermoDB"; }

    /**
     * Constructs an empty database with a given standard state temperature and 
     * pressure.
     */
    ThermoDB(double sst, double ssp);
    
    /**
     * Loads the species from the database represented by the concrete class
     * which extends this ThermoDB using the given SpeciesListDescriptor.
     */
    bool load(const SpeciesListDescriptor& descriptor);
    
    /**
     * Desctructor.
     */
    virtual ~ThermoDB() {};
    
    /**
     * Returns the element list.
     */
    const std::vector<Element>& elements() const {
        return m_elements;
    }
    
    /**
     * Returns the species list.
     */
    const std::vector<Species>& species() const {
        return m_species;
    }
    
    /**
     * Returns the standard state temperature in K employeed in this database.
     */
    double standardTemperature() const {
        return m_sst;
    }
    
    /**
     * Returns the standard state pressure in Pa employeed in this database.
     */
    double standardPressure() const {
        return m_ssp;
    }
    
    /**
     * Returns true if the thermodynaic data in this database is valid at the
     * given temperature for the given species index.
     */
    virtual bool speciesThermoValidAtT(const size_t i, const double T) const {
        return true;
    }
    
    /**
     * Computes the unitless internal specific heat.  The default behavior is to
     * simply provide \f$c_{p,\text{int}} = c_p(T) - 2.5\f$ for all species.
     * @param T temperature assigned to internal energy modes
     */
    virtual void cpint(double T, double* const cp);

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
        double* const cp = NULL, double* const cpt = NULL,
        double* const cpr = NULL, double* const cpv = NULL,
        double* const cpel = NULL) = 0;
    
    /**
     * Computes the species vibrational specific heats at the given temperature
     * nondimensionalized by the universal gas constant.
     */
    virtual void cpv(double Tv, double* const p_cpv);

    /**
     * Returns the species electronic specific heats at the given temperature
     * nondimensionalized by the universal gas constant.
     */
    virtual void cpel(double Tel, double* const p_cpel);

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
     * Returns the species vibrational enthalpy at the given temperature
     * nondimensionalized by the universal gas constant times the temperature.
     */
    virtual void hv(double Tv, double* const p_hv);

    /**
     * Returns the species electronic enthalpy at the given temperature
     * nondimensionalized by the universal gas constant times the temperature.
     */
    virtual void hel(double Tel, double* const p_hel);

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

    /**
     * Implemented by the concrete class which fills the species vector with all
     * of the species available in the database.  Note that this is only called
     * once in the constructor.
     */
    virtual void loadAvailableSpecies(std::list<Species>& species) = 0;
    
    /**
     * Implemented by the concrete class which loads any data from the database
     * required to compute the necessary thermodynamic quantities.  This is 
     * function is called during a call to load() after the list of species
     * that are needed has been determined.
     */
    virtual void loadThermodynamicData() = 0;
    
private:

    double m_sst;
    double m_ssp;
    
    std::vector<Element> m_elements;
    std::vector<Species> m_species;
    
}; // class ThermoDB

    } // namespace Thermodynamics
} // namespace Mutation

#endif // THERMO_THERMODB_H

