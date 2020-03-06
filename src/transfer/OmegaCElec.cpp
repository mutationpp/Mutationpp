/*
 * @file OmegaCElec.cpp
 *
 * @brief Implementation of OmegaCElec.
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

namespace Mutation {
    namespace Transfer {

/**
 * Represents a coupling between chemistry and electronic energy modes.
 *
 * \f[ \sum_{i \in \varepsilon } \dot{\omega}_i e^E_i \f].
 *
 * where the index \f$ \varepsilon \f$ represents the set oh heavy species
 *
 * More information about the above model can be found in @cite Panesi 2008
 */
class OmegaCElec : public TransferModel
{
public:
	OmegaCElec(Mutation::Mixture& mix)
		: TransferModel(mix)
	{
		mp_wrk1 = new double [mix.nSpecies()];
		mp_wrk2 = new double [mix.nSpecies()];
	};

	~OmegaCElec()
	{
		delete [] mp_wrk1;
		delete [] mp_wrk2;
	};

	double source()
	{
		m_mixture.speciesHOverRT(NULL, NULL, NULL, NULL, mp_wrk1, NULL);
		m_mixture.netProductionRates(mp_wrk2);

		double sum = 0.0;

		for(int i = 0 ; i < m_mixture.nSpecies(); ++i)
			sum += mp_wrk1[i]*mp_wrk2[i]/m_mixture.speciesMw(i);

		return (sum*m_mixture.T()*RU);
	}

private:
	double* mp_wrk1;
	double* mp_wrk2;
};

// Register the transfer model
Mutation::Utilities::Config::ObjectProvider<
    OmegaCElec, TransferModel> omegaCElec("OmegaCElec");

    } // namespace Transfer
} // namespace Mutation
