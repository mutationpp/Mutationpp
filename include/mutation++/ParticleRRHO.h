/**
 * @file ParticleRRHO.h
 *
 * @brief Declaration of ParticleRRHO class.
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

#ifndef THERMO_PARTICLE_RRHO_H
#define THERMO_PARTICLE_RRHO_H

#include <vector>
#include <cstdlib>

// Forward declartion of XmlElement
namespace Mutation { 
    namespace Utilities {
        namespace IO {
class XmlElement;
        }
    }
}

namespace Mutation {
    namespace Thermodynamics {

/**
 * Stores the model parameters for a particle using the Rigid-Rotator Harmonic-
 * Oscillator thermodynamic model.  Note that this class is used to load the 
 * parameters from the Species.xml input file but only to provide the 
 * information to the corresponding database model.
 */
class ParticleRRHO
{
public:
    
    /**
     * Loads the RRHO parameters from an XmlElement.
     *
     * @todo Implement error checking, units conversion, and ParticleType
     * checking.
     */
    ParticleRRHO(const Mutation::Utilities::IO::XmlElement& xml_element);
    
    /**
     * Constructs a new set of RRHO parameters using the parameters from the 
     * given object but with only one electronic energy level which is 
     * specified.
     */
    ParticleRRHO(const ParticleRRHO& p_rrho, const size_t level);
    
    /**
     * Returns the formation enthalpy in J/mol.
     */
    double formationEnthalpy() const {
        return m_hform;
    }
    
    /**
     * Returns the characteristic rotational temperature in K.
     */
    double rotationalTemperature() const {
        return m_rotational_t;
    }
    
    /**
     * Returns the steric factor.
     */
    int stericFactor() const {
        return m_steric;
    }
    
    /**
     * Returns the linearity of the molecule.
     */
    int linearity() const {
        return m_linearity;
    }
    
    /**
     * Returns the number of electronic energy levels for the heavy particle.
     */
    int nElectronicLevels() const {
        return m_electronic_energies.size();
    }
    
    /**
     * Returns the (degeneracy, energy) pair for the electronic level i where
     * the energy is specified in K.
     */
    const std::pair<int, double>& electronicEnergy(const int i) const {
        return m_electronic_energies[i];
    }
    
    /**
     * Returns the number of vibrational energy levels for the molecule.
     */
    int nVibrationalLevels() const {
        return m_vibrational_energies.size();
    }
    
    /**
     * Returns the vibrational energy of level i in K.
     */
    double vibrationalEnergy(const int i) const {
        return m_vibrational_energies[i];
    }
    
private:
    
    double m_hform;
    int    m_steric;
    int    m_linearity;
    double m_rotational_t;
    
    std::vector< std::pair<int, double> > m_electronic_energies;
    std::vector<double> m_vibrational_energies;
    
}; // class ParticleRRHO

    } // namespace Thermodynamics
} // namespace Mutation

#endif // THERMO_PARTICLE_RRHO_H


