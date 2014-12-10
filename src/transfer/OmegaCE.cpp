/**
 * @file OmegaCE.cpp
 *
 * @brief Implementation of OmegaCE.
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

namespace Mutation {
    namespace Transfer {

/**
 * Represents a coupling between chemistry and electrons.
 **/
class OmegaCE : public TransferModel
{
public:
	OmegaCE(Mutation::Mixture& mix)
		: TransferModel(mix)
	{
		mp_wrk1 = new double [mix.nSpecies()];
		mp_wrk2 = new double [mix.nSpecies()];
	};

	~OmegaCE()
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
			  std::cerr << "The selected Electron-Chemistry-Electron model is not implemented yet";
			  return 0.0;
	  }
	};

private:
	double* mp_wrk1;
	double* mp_wrk2;

	double const compute_source_Candler();
};


 double const OmegaCE::compute_source_Candler()
 {
	 m_mixture.netProductionRates(mp_wrk1);
	 return mp_wrk1[0]*1.5*RU*m_mixture.Te()/m_mixture.speciesMw(0);
 }

 // Register the transfer model
 Mutation::Utilities::Config::ObjectProvider<
     OmegaCE, TransferModel> omegaCE("OmegaCE");
      
    } // namespace Transfer
} // namespace Mutation
