/**
 * @file OmegaET.cpp
 *
 * @brief Implementation of OmegaET.
 */

/*
 * Copyright 2014 von Karman Institute for Fluid Dynamics (VKI)
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

namespace Mutation {
    namespace Transfer {
      
class OmegaET : public TransferModel
{
public:
	OmegaET(Mutation::Mixture& mix)
		: TransferModel(mix), mp_collisions(mix.collisionData())
	{ }

	/**
	 * Computes the source term of the Electron-Heavy Translational energy transfer in \f$ [J/(m^3\cdot s)] \f$
	 * using a Landau-Teller formula:
	 * \f[ \Omega^{ET} = \rho_e \frac{e^T_e\left(T\right)-e^T_e\left(T_{e}\right)}
	 *      {\tau^{ET}}  \f]
	 *
	 * The relaxation time \f$ \tau^{ET} \f$ is given by the expression:
	 *
	 * \f[ \tau^{ET} = \left( \sum_{j \in \mathcal{H}} \frac{2 M_e}{ M_i} \nu_{ei} \right)^{-1} \f]
	 *
	 */

	double source()
	{
		const double * p_X = m_mixture.X();
        double nd = m_mixture.numberDensity();
        double T = m_mixture.T();
        double Te = m_mixture.Te();

        double tau = compute_tau_ET();
        return 1.5E0*KB*nd*p_X[0]*(T-Te)/tau;
	}

private:
	Transport::CollisionDB* mp_collisions;

	double const compute_tau_ET();
};


double const OmegaET::compute_tau_ET()
{
	const double T = m_mixture.T();
    const double Te = m_mixture.Te();
    const double nd = m_mixture.numberDensity();
    const double ns = m_mixture.nSpecies();
    const double *const p_X = m_mixture.X();
    const double mwel = m_mixture.speciesMw(0);
    double mwis, nu, sum, tau;

    // Collisional integrals
    const Eigen::MatrixXd& Q11 = mp_collisions->Q11(T, Te, nd, p_X);

    // Electron velocity
    double ve = sqrt(RU*8.0*Te/(PI*mwel));

    // Compute mass averaged collision frequency
    sum = 0.0;
    for (int is = 1; is < ns; ++is) {
        mwis = m_mixture.speciesMw(is);
        nu = ve * p_X[is] * nd * Q11(is);
        sum += nu/mwis;
    }
    sum *= mwel*8.0/3.0;

    // Return tau
    tau = 1.0/sum;
    return tau;
}

// Register the transfer model
Mutation::Utilities::Config::ObjectProvider<
    OmegaET, TransferModel> omegaET("OmegaET");

    } // namespace Transfer
} // namespace Mutation
