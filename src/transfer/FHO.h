/**
 * @file FHO.h
 *
 * @brief Declaration of classes related to the Forced Harmonic Oscillator
 * (FHO) model.
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

#ifndef FHO_H
#define FHO_H

#include <string>
#include <vector>
#include <cmath>
#include "Thermodynamics.h"

namespace Mutation {
    namespace Transfer {

/** 
 * Data storage for single FHO collision pair.
 */
class FHOPartner
{
public:

    /**
     * Construct using given values of a, b, and mu.
     */
    FHOPartner(double beta, double E_Morse, double svt)
        : m_beta(beta), m_E_Morse(E_Morse), m_svt(svt)
    {}
    
    /**
     * Construct using default values of beta, E_Morse and svt.
     */
    FHOPartner()
        : m_beta(0.),
          m_E_Morse(0.),
          m_svt(0.)
    { }
    
    double beta() const {return m_beta; }
    double E_Morse() const {return m_E_Morse; }
    double svt() const {return m_svt; }

private:

    double m_beta;
    double m_E_Morse;
    double m_svt;
};

/**
 * Data storage for all collision partners of a vibrating species.
 */
class FHOVibrator
{
public:

    /**
     * Load available data from VT_FHO.xml and use defaults when data does not 
     * exist.
     */
    FHOVibrator(
        const Mutation::Utilities::IO::XmlElement& node,
        const Mutation::Thermodynamics::Thermodynamics& thermo);
    
    /**
     * Use default data for all collision pairs with this vibrator.
     */
    FHOVibrator(
        const std::string& name,
        const Mutation::Thermodynamics::Thermodynamics& thermo);
    
    /**
     * Returns the (beta, E_Morse, svt) data for the i'th heavy species colliding
     * with this vibrator.
     */
    const FHOPartner& operator[](int i) const {
        assert(i >= 0);
        assert(i < m_partners.size());
        return m_partners[i];
    }
    
private:

    std::vector<FHOPartner> m_partners;
};

/**
 * Data storage for all of the information needed to evaluate the Millikan-White
 * vibration-translation energy transfer.
 */
class FHO
{
public:

    /**
     * Loads all of the data associated with the Millikan-White model.
     */
    FHO(const Mutation::Thermodynamics::Thermodynamics& thermo);
    
    /**
     * Returns the number of vibrators being kept in this model.
     */
    int nVibrators() const {
        return m_vibrators.size();
    }
    
    /**
     * Returns the i'th vibrator data.
     */
    const FHOVibrator& operator[](int i) const {
        assert(i >= 0);
        assert(i < m_vibrators.size());
        return m_vibrators[i];
    }

private:

    std::vector<FHOVibrator> m_vibrators;
};

    } // namespace Kinetics
} // namespace Mutation

#endif // FHO_H
