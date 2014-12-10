/*
 * @file OmegaCElec.cpp
 *
 * @brief Implementation of OmegaCElec.
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

namespace Mutation {
    namespace Transfer {

/**
 * Represents a coupling between chemistry and electronic energy modes.
 */
class OmegaCElec : public TransferModel
{
   /**
	* The necessary functions for the computation of the  source terms for
	* Electronic-Chemistry-Electronic coupling are implemented.
	*
	*/

public:
	OmegaCElec(Mutation::Mixture& mix)
		: TransferModel(mix)
	{
		m_ns = m_mixture.nSpecies();
		mp_wrk1 = new double [m_ns];
		mp_wrk2 = new double [m_ns];
	};

	~OmegaCElec()
	{
		delete [] mp_wrk1;
		delete [] mp_wrk2;
	};

	double source()
	{
	  static int i_transfer_model = 0;
	  switch (i_transfer_model){
		  case 0:
			  return compute_source_Candler();
		  break;
		  default:
			  std::cerr << "The selected Electronic-Chemistry-Electronic model is not implemented yet";
			  return 0.0;
	  }
	};

private:
	int m_ns;
	double* mp_wrk1;
	double* mp_wrk2;

	double const compute_source_Candler();
};


double const OmegaCElec::compute_source_Candler()
{
	/**
	* Non-preferential Model with
	*
	* \f[ \Omega^{CElec}_{mi} = \Sumc_1 e^{Elec}_{mi} \dot{\omega}_i\f]
	*
	* with \f$ c_1 \f$ equal to 1 for non-preferential models.
	*
	*/

	// Getting Electronic Energy
	m_mixture.speciesHOverRT(NULL, NULL, NULL, NULL, mp_wrk1, NULL);

	// Getting Production Rate
	m_mixture.netProductionRates(mp_wrk2);

	// Inner Product
	double c1 = 1.0E0;
	double sum = 0.E0;

	for(int i = 0 ; i < m_ns; ++i)
		sum += mp_wrk1[i]*mp_wrk2[i] *m_mixture.T()*RU/m_mixture.speciesMw(i);

	return(c1*sum);
}

// Register the transfer model
Utilities::Config::ObjectProvider<
    OmegaCElec, TransferModel> omegaCElec("OmegaCElec");

    } // namespace Transfer
} // namespace Mutation
