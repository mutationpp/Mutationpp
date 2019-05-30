/**
 * @file VSS.h
 *
 * @brief Declaration of classes related to Variable Soft Sphere (VSS) model.
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

#ifndef VSS_H
#define VSS_H

#include <string>
#include <vector>
#include <cmath>
#include "Thermodynamics.h"

namespace Mutation {
    namespace Transfer {

/** 
 * Data storage for single VSS collision pair.
 */
class VSSPartner
{
public:

    /**
     * Construct using given values of dref, omega
     */
    VSSPartner(double dref, double omega, double mu)
        : m_dref(dref), m_omega(omega), m_mu(mu)
    {}
    
    /**
     * Construct using default values of dref and omega.
     */
    VSSPartner()
        : m_dref(0.),
          m_omega(0.),
	  m_mu(0.)
    { }
    
    double dref() const {return m_dref; }
    double omega() const {return m_omega; }
    double mu() const {return m_mu;}

private:

    double m_dref;
    double m_omega;
    double m_mu;
};

/**
 * Data storage for all collision partners of a vibrating species.
 */
class VSSVibrator
{
public:

    /**
     * Load available data from VT_VSS.xml and use defaults when data does not 
     * exist.
     */
    VSSVibrator(
        const Mutation::Utilities::IO::XmlElement& node,
        const Mutation::Thermodynamics::Thermodynamics& thermo);
    
    /**
     * Use default data for all collision pairs with this vibrator.
     */
    VSSVibrator(
        const std::string& name,
        const Mutation::Thermodynamics::Thermodynamics& thermo);
    
    /**
     * Returns the (dref, omega) data for the i'th heavy species colliding with
     * this vibrator.
     */
    const VSSPartner& operator[](int i) const {
        assert(i >= 0);
        assert(i < m_partners.size());
        return m_partners[i];
    }
    
    /**
     * Returns the limiting cross-section for vibrational excitation at 50,000 K
     * for this vibrator. Defaults to 3.0E-21 unless given differently in the 
     * data file.
     */
    double omega() const {
        return m_omegav;
    }
    
    /**
     * Returns the index of this vibrator in the mixture species list.
     */
    int index() const {
        return m_index;
    }

    /**
     * Returns the oscillation mass of this vibrator in the mixture species list.
     */
    double osc_mass() const {
        return m_osc_massv;
    }

private:

    /**
     * Finds the characteristic vibrational temperature for this vibrator.
     */
    static double loadThetaV(const std::string& name);

    /**
     * Finds the dissociation temperature for this vibrator.
     */
    static double loadDissV(const std::string& name);

    /**
     * Finds the doscillator mass for this vibrator.
     */
    static double loadOscMassV(const std::string& name);

    std::vector<VSSPartner> m_partners;
    double m_omegav;
    double m_thetav; // not used since vibr. spectrum is hardcoded
    double m_dissv; // not used since vibr. spectrum is hardcoded
    double m_osc_massv; // oscillator mass
    int m_index;
};

/**
 * Data storage for all of the information needed to evaluate the VSS
 * vibration-translation energy transfer.
 */
class VSS
{
public:

    /**
     * Loads all of the data associated with the VSS model.
     */
    VSS(const Mutation::Thermodynamics::Thermodynamics& thermo);
    
    /**
     * Returns the number of vibrators being kept in this model.
     */
    int nVibrators() const {
        return m_vibrators.size();
    }
    
    /**
     * Returns the i'th vibrator data.
     */
    const VSSVibrator& operator[](int i) const {
        assert(i >= 0);
        assert(i < m_vibrators.size());
        return m_vibrators[i];
    }

private:

    std::vector<VSSVibrator> m_vibrators;
};

    } // namespace Kinetics
} // namespace Mutation

#endif // VSS_H
