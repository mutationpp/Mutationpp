/**
 * @file MillikanWhiteRelaxationTimeGnoffo.cpp
 *
 * @brief Relaxation time for VT source term 
 * according to Gnoffo mixing rule
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
 * @brief Gnoffo  mixing rule with Park high-temperature correction.
 */
class MillikanWhiteRelaxationTimeGnoffo : public MillikanWhiteRelaxationTime
{
public:
    MillikanWhiteRelaxationTimeGnoffo(ARGS) {}

    double relaxationTime(
        const Mutation::Thermodynamics::Thermodynamics& thermo,
        const MillikanWhiteModelData& data) const override
    {
        const double T_fac = std::pow(thermo.T(), -1.0/3.0);
        const double p_atm = thermo.P() / ONEATM;

        const Eigen::Map<const Eigen::ArrayXd> Xh(thermo.X() + (thermo.hasElectrons() ? 1 : 0), thermo.nHeavy());
        const double tau_mw = (Xh*(data.a()*(T_fac - data.b()) - 18.421).exp()).sum()/(Xh.sum()*p_atm);
        const double ni = thermo.numberDensity() * thermo.X()[data.speciesIndex()];
        const double ci = std::sqrt(8*RU*thermo.T()/(PI*data.molecularWeight()));
        const double tau_park = 1.0/(ni * ci * data.limitingCrossSection(thermo.T()));

        return tau_mw + tau_park;
    }
};

// Register this algorithm
    Config::ObjectProvider<MillikanWhiteRelaxationTimeGnoffo, MillikanWhiteRelaxationTime> 
    mw_gnoffo("Gnoffo");

} // namespace Transfer
} // namespace Mutation
