/**
 * @file OmegaET.cpp
 *
 * @brief Implementation of OmegaET.
 */

/*
 * Copyright 2014-2017 von Karman Institute for Fluid Dynamics (VKI)
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

#include "Mixture.h"
#include "TransferModel.h"

#include <cmath>

#include <eigen3/Eigen/Dense>
using namespace Eigen;

namespace Mutation {
    namespace Transfer {

/**
 * Computes the source term of the electron-heavy translational energy
 * transfer in \f$J/(m^3 s)\f$, using a Landau-Teller formula:
 * \f[
 * \Omega^{ET} = \rho_e \frac{e^T_e\left(T\right)-e^T_e\left(T_{e}\right)}
 *     {\tau^{ET}},
 * \f]
 * where relaxation time \f$\tau^{ET}\f$ is given by
 * \f[
 * \frac{1}{\tau^{ET}} = \sum_{j\in\mathcal{H}} \frac{2m_e}{m_i}\nu_{ei}.
 * \f]
 *
 * @todo Check that the equations in the documentation match the implementation.
 * @todo Merge compute_tau_ET() into the source() function directly.
 */
class OmegaET : public TransferModel
{
public:
	OmegaET(TransferModel::ARGS mix) :
	    TransferModel(mix), m_collisions(mix.collisionDB()),
		 m_has_electrons(mix.hasElectrons())
	{ }

	double source()
	{
	    if (!m_has_electrons)
		    return 0.0;

	    const double * p_X = m_mixture.X();
        double nd = m_mixture.numberDensity();
        double T = m_mixture.T();
        double Te = m_mixture.Te();

        double tau = compute_tau_ET();
        return 1.5*KB*nd*p_X[0]*(T-Te)/tau;
	}

private:
	Transport::CollisionDB& m_collisions;
	bool m_has_electrons;

	double const compute_tau_ET();
};


double const OmegaET::compute_tau_ET()
{
	const double T = m_mixture.T();
    const double Te = m_mixture.Te();
    const double nd = m_mixture.numberDensity();
    const double ns = m_mixture.nSpecies();
    Map<const ArrayXd> X(m_mixture.X(),ns);

    // CollisionDB data
    const ArrayXd& Q11ei = m_collisions.Q11ei();
    const ArrayXd& mass  = m_collisions.mass();

    // Electron velocity
    double ve = sqrt(KB*8.*Te/(PI*mass(0)));

    // Tau
    return 1.0/(mass(0)*8./3.*ve*nd*(X*Q11ei/mass).tail(ns-1).sum());
}

// Register the transfer model
Mutation::Utilities::Config::ObjectProvider<
    OmegaET, TransferModel> omegaET("OmegaET");

    } // namespace Transfer
} // namespace Mutation
