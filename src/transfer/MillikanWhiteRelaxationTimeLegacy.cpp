/**
 * @file MillikanWhiteRelaxationTimeLegacy.cpp
 *
 * @brief Relaxation time for VT source term
 * legacy relaxation time
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


#include "MillikanWhiteRelaxationTime.h"
#include "MillikanWhite.h"
#include "Thermodynamics.h"
#include "Utilities.h"
#include "AutoRegistration.h"

#include <Eigen/Dense>

using namespace Mutation::Utilities;

namespace Mutation {
namespace Transfer {

/**
 * @brief Millikan-White mixing rule with Park high-temperature correction.
 */
class MillikanWhiteRelaxationTimeLegacy : public MillikanWhiteRelaxationTime
{
public:
    MillikanWhiteRelaxationTimeLegacy(ARGS) {}

    double relaxationTime(
        const Mutation::Thermodynamics::Thermodynamics& thermo,
        const MillikanWhiteModelData& data) const override
    {    
        const double T_fac = std::pow(thermo.T(), -1.0/3.0);
        const double p_atm = thermo.P() / ONEATM;

        const Eigen::Map<const Eigen::ArrayXd> Xh(thermo.X() + (thermo.hasElectrons() ? 1 : 0), thermo.nHeavy());
        const Eigen::ArrayXd tau_tau_mw = (data.a()*(T_fac - data.b()) - 18.421).exp()/p_atm;
        const Eigen::ArrayXd tau_park = ((PI*data.mu()*KB*thermo.T()) /
                                        (8.0*NA)).sqrt()/thermo.P()/data.limitingCrossSection(thermo.T());
        return Xh.sum() / (Xh/(tau_tau_mw+tau_park)).sum();
    }
};

// Register this algorithm
Config::ObjectProvider<
    MillikanWhiteRelaxationTimeLegacy, MillikanWhiteRelaxationTime>
    mw_Legacy("Legacy");

} // namespace Transfer
} // namespace Mutation
