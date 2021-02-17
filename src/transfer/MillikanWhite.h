/**
 * @file MillikanWhite.h
 *
 * @brief Declaration of classes related to Millikan and White model.
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


#ifndef MILLIKAN_WHITE_H
#define MILLIKAN_WHITE_H

#include <string>
#include <vector>
#include <cmath>
#include "Thermodynamics.h"

namespace Mutation {
    namespace Transfer {

/** 
 * Data storage for single Millikan-White collision pair.
 */
class MillikanWhitePartner
{
public:

    /**
     * Construct using given values of a, b, and mu.
     */
    MillikanWhitePartner(double a, double b, double mu)
        : m_a(a), m_b(b), m_mu(mu)
    {}
    
    /**
     * Construct using default values of a and b based on mu and theta.
     */
    MillikanWhitePartner(double mu, double theta)
        : m_a(1.16E-3*std::sqrt(mu*1000.0)*std::pow(theta, 4.0/3.0)),
          m_b(0.015*std::pow(mu*1000.0, 0.25)),
          m_mu(mu)
    { }
    
    double a()  const {return m_a; }
    double b()  const {return m_b; }
    double mu() const {return m_mu;}

private:

    double m_a;
    double m_b;
    double m_mu;
};

/**
 * Data storage for all collision partners of a vibrating species.
 */
class MillikanWhiteVibrator
{
public:

    /**
     * Load available data from VT.xml and use defaults when data does not 
     * exist.
     */
    MillikanWhiteVibrator(
        const Mutation::Utilities::IO::XmlElement& node,
        const Mutation::Thermodynamics::Thermodynamics& thermo);
    
    /**
     * Use default data for all collision pairs with this vibrator.
     */
    MillikanWhiteVibrator(
        const std::string& name,
        const Mutation::Thermodynamics::Thermodynamics& thermo);
    
    /**
     * Returns the (a,b,mu) data for the i'th heavy species colliding with this
     * vibrator.
     */
    const MillikanWhitePartner& operator[](int i) const {
        assert(i >= 0);
        assert(i < m_partners.size());
        return m_partners[i];
    }
    
    /**
     * Returns the limiting cross-section for vibrational excitation at 50,000 K
     * for this vibrator.  Defaults to 3.0E-21 unless given differently in the 
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

private:

    /**
     * Finds the characteristic vibrational temperature for this vibrator.
     */
    static double loadThetaV(const std::string& name);

    std::vector<MillikanWhitePartner> m_partners;
    double m_omegav;
    double m_thetav;
    int m_index;
};


/**
 * Data storage for all of the information needed to evaluate the Millikan-White
 * vibration-translation energy transfer.
 */
class MillikanWhite
{
public:
    /**
     * Loads all of the data associated with the Millikan-White model.
     */
    MillikanWhite(const Mutation::Thermodynamics::Thermodynamics& thermo);
    
    /**
     * Returns the number of vibrators being kept in this model.
     */
    int nVibrators() const {
        return m_vibrators.size();
    }
    
    /**
     * Returns the i'th vibrator data.
     */
    const MillikanWhiteVibrator& operator[](int i) const {
        assert(i >= 0);
        assert(i < m_vibrators.size());
        return m_vibrators[i];
    }

private:

    std::vector<MillikanWhiteVibrator> m_vibrators;
};

    } // namespace Kinetics
} // namespace Mutation

#endif // MILLIKAN_WHITE_H
