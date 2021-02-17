/**
 * @file OmegaET.cpp
 *
 * @brief Implementation of OmegaET.
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

#include "Mixture.h"
#include "TransferModel.h"

#include <cmath>

#include <Eigen/Dense>
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
 *
 * \\f]
 * and the collision frequency \f$\nu_{ei}\f$ is given by
 * \f[
 * \nu_{ei} = (8/3)v_e n_j \overline{\Omega}_{ej}^{11}
 *
 *\\f]
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

        const double ns = m_mixture.nSpecies();
        Map<const ArrayXd> X(m_mixture.X(),ns);

        // CollisionDB data
        const ArrayXd& Q11ei = m_collisions.Q11ei();
        const ArrayXd& mass  = m_collisions.mass();

        // Electron velocity
        double ve = sqrt(KB*8.*Te/(PI*mass(0)));
        double tau = 1.0/(mass(0)*8./3.*ve*nd*(X*Q11ei/mass).tail(ns-1).sum());

        return 1.5*KB*nd*p_X[0]*(T-Te)/tau;
	}

private:
	Transport::CollisionDB& m_collisions;
	bool m_has_electrons;

};


// Register the transfer model
Mutation::Utilities::Config::ObjectProvider<
    OmegaET, TransferModel> omegaET("OmegaET");

    } // namespace Transfer
} // namespace Mutation
