/**
 * @file ParticleRRaHO.h
 *
 * @brief Declaration of ParticleRRaHO class.
 */

/*
 * Copyright 2014-2018 von Karman Institute for Fluid Dynamics (VKI)
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

#ifndef THERMO_PARTICLE_RRaHO_H
#define THERMO_PARTICLE_RRaHO_H

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
 * Stores the model parameters for a particle using the Rigid-Rotator anHarmonic-
 * Oscillator thermodynamic model.  Note that this class is used to load the 
 * parameters from the Species.xml input file but only to provide the 
 * information to the corresponding database model.
 */
class ParticleRRaHO
{
public:
    
    /**
     * Loads the RRaHO parameters from an XmlElement.
     *
     * @todo Implement error checking, units conversion, and ParticleType
     * checking.
     */
    ParticleRRaHO(const Mutation::Utilities::IO::XmlElement& xml_element);
    
    /**
     * Constructs a new set of RRaHO parameters using the parameters from the 
     * given object but with only one electronic energy level which is 
     * specified.
     */
    ParticleRRaHO(const ParticleRRaHO& p_rraho, const size_t level);
    
    /**
     * Returns the formation enthalpy in J/mol.
     */
    double formationEnthalpy() const {
        return m_hform;
    }
    
    /**
     * Returns the formation energy in J/mol.
     */
    double formationEnergy() const {
        return m_eform;
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
     * Returns the number of vibrational energy levels for the molecule
     * at the electronic level e.
     */
    int nVibrationalLevels(const int e) const {
        return vib_energy[e].size();
    }

    /**
     * Returns the vibrational energy of level i in K.
     */
    double vibrationalEnergy(const int i) const {
        return m_vibrational_energies[i];
    }
    
    /**
     * Returns the vibrational energy of electronic level e
     * and vibrational level v in K.
     */
    double vibrationalEnergy(const int e, const int v) const {
        return vib_energy[e][v];
    }

    // Returns the number of rotational energy levels for the molecule.
    int nRotationalLevels() const { return m_rotational_energies.size(); }

    /**
     * Returns the number of rotational energy levels for the molecule
     * at electronic level e and vibrational level v.
     */
    int nRotationalLevels(const int e, const int v) const {
        return rot_energy[e][v].size();
    }

    // Returns the rotational energy at the generic index i.
    double rotationalEnergy(const int i) const {
        return m_rotational_energies[i];
    }

    /**
     * Returns the rotational energy at the electronic level e
     * vibrational level v and rotational level r.
     */
    double rotationalEnergy(const int e, const int v, const int r) const {
        return rot_energy[e][v][r];
    }

private:
    
    double m_hform;
    double m_eform;
    int    m_steric;
    int    m_linearity;
    double m_rotational_t;

    double      m_mass;
    double      m_diameter;
    int         m_charge;
    double      m_formation_energy;
    double      m_eps_LJ;
    double      m_inertia_moment;
    double      m_parker_zeta_inf;
    double      m_reduced_oscillator_mass;
    double      m_r_e;
    double      m_mA_mAB_ratio;
    double      m_mB_mAB_ratio;
    double      m_ionization_potential;
    
    std::vector< std::pair<int, double> > m_electronic_energies;
    std::vector<double> m_electronical_energies;

    std::vector<double> m_vibrational_energies;
    std::vector<int> m_vibrational_levels;
    std::vector<std::vector<double>> vib_energy;
    std::vector<double> m_dissociation_energies;
    std::vector<double> m_vibrational_temperatures;

    std::vector<int> m_rotational_levels;
    std::vector<double> m_rotational_energies;
    std::vector<std::vector<std::vector<double>>> rot_energy;

    std::vector<double> m_we;   double we;
    std::vector<double> m_wexe; double wexe;
    std::vector<double> m_weye; double weye;
    std::vector<double> m_weze; double weze;
    std::vector<double> m_ae;   double ae;
    std::vector<double> m_be;   double be;
    std::vector<double> m_de;   double de;
    
}; // class ParticleRRaHO

    } // namespace Thermodynamics
} // namespace Mutation

#endif // THERMO_PARTICLE_RRaHO_H
