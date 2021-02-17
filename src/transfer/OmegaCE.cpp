/**
 * @file OmegaCE.cpp
 *
 * @brief Implementation of OmegaCE.
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

namespace Mutation {
    namespace Transfer {

/**
 * Represents a coupling between chemistry and electrons by
 *  \f$ \dot{\omega _e} e^T_e \f$.
 **/

class OmegaCE : public TransferModel
{
public:
	OmegaCE(Mutation::Mixture& mix)
		: TransferModel(mix)
	{
		mp_wrk1 = new double [mix.nSpecies()];
	}

	~OmegaCE()
	{
		delete [] mp_wrk1;
	}

	double source()
	{
		static const double cv = 1.5*RU/m_mixture.speciesMw(0);
		m_mixture.netProductionRates(mp_wrk1);
		return mp_wrk1[0]*cv*m_mixture.Te();
	}

private:
	double* mp_wrk1;
};

// Register the transfer model
Mutation::Utilities::Config::ObjectProvider<
    OmegaCE, TransferModel> omegaCE("OmegaCE");
      
    } // namespace Transfer
} // namespace Mutation
