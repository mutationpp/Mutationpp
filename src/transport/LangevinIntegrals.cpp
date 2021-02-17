/**
 * @file LangevinIntegrals.cpp
 *
 * Provides Langevin collision integral types.
 *
 * @author J.B. Scoggins
 */

/*
 * Copyright 2015-2020 von Karman Institute for Fluid Dynamics (VKI)
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

#include "AutoRegistration.h"
#include "Constants.h"
#include "CollisionIntegral.h"
#include "CollisionPair.h"
#include "Thermodynamics.h"
#include "XMLite.h"
using namespace Mutation::Thermodynamics;

using namespace Mutation::Utilities;
using namespace Mutation::Utilities::Config;
using namespace Mutation::Utilities::IO;

#include <cassert>
using namespace std;

namespace Mutation {
    namespace Transport {

/**
 * Represents a Langevin (Polarization) collision integral.  Supported integrals
 * are Q(l,s) for l = {1,2,3} and s = {1,2,3,4,5}.  Requires the dipole-
 * polarizability to be known for the neutral species.
 */
class LangevinColInt : public CollisionIntegral
{
public:

    LangevinColInt(CollisionIntegral::ARGS args) :
        CollisionIntegral(args), m_fac(0.0), m_alpha(0.0)
    {
        // Make sure this is an ion-neutral interaction
        args.xml.parseCheck(args.pair.type() == ION_NEUTRAL,
            "must be ion-neutral interaction to use the Langevin integral.");

        // Determine the charge of the ion and the name of the neutral
        double z = args.pair.sp1().charge();
        string neutral = args.pair.sp2Name();

        if (z == 0.0) {
            z = args.pair.sp2().charge();
            neutral = args.pair.sp1Name();
        }

        // Determine (l,s)
        string kind = args.kind;
        bool format_check = (
            kind.size() == 3 && kind[0] == 'Q' &&
            kind[1] > '0' && kind[1] < '4'     &&
            kind[2] > '0' && kind[2] < '6');
        args.xml.parseCheck(format_check,
            "integral must be 'Qls' with 0 < l < 4 and 0 < s < 6.");
        int l = int(kind[1])-48;
        int s = int(kind[2])-48;

        // Load the dipole polarizability of the neutral in m^3
        m_alpha = loadSpeciesParameter(
            args.xml.document()->root(), "dipole-polarizabilities", neutral);

        // Compute the factor that is multiplied by sqrt(alpha/T)
        const double sfac [5] = {0.66467, 0.55389, 0.48466, 0.43619, 0.39984};
        const double lfac [3] = {1.1046, 1.1538, 1.2754};
        m_fac = PI*lfac[l-1]*sfac[s-1]*std::sqrt(z*z*QE*QE/(TWOPI*EPS0*KB));
    }

    // Allow tabulation of this integral type
    bool canTabulate() const { return true; }

    // Make sure we could correctly load the data
    bool loaded() const { return m_alpha * m_fac > 0.0; }

private:

    double compute_(double T) { return m_fac * std::sqrt(m_alpha / T); }

    /**
     * Returns true if the constant value is the same.
     */
    bool isEqual(const CollisionIntegral& ci) const {
        const LangevinColInt& compare = dynamic_cast<const LangevinColInt&>(ci);
        return (m_fac == compare.m_fac && m_alpha == compare.m_alpha);
    }

private:
    double m_fac;
    double m_alpha;
};

// Register the "Langevin" CollisionIntegral
ObjectProvider<LangevinColInt, CollisionIntegral> Langevin_ci("Langevin");

    } // Mutation
} // Transport
