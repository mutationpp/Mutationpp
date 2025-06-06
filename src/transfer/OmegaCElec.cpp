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
                m_ns = mix.nSpecies();
                m_nt = mix.nEnergyEqns();
		mp_wrk1 = new double [mix.nSpecies()];
		mp_wrk2 = new double [mix.nSpecies()];
		mp_wrk3 = new double [mix.nSpecies()*mix.nSpecies()];
		mp_wrk4 = new double [mix.nSpecies()];
		mp_wrk5 = new double [mix.nSpecies()];
		mp_wrk6 = new double [mix.nSpecies()];
	};

	~OmegaCElec()
	{
		delete [] mp_wrk1;
		delete [] mp_wrk2;
		delete [] mp_wrk3;
		delete [] mp_wrk4;
		delete [] mp_wrk5;
		delete [] mp_wrk6;
	};

	double source()
	{
		mixture().speciesHOverRT(NULL, NULL, NULL, NULL, mp_wrk1, NULL);
		mixture().netProductionRates(mp_wrk2);

		double sum = 0.0;

		for(int i = 0 ; i < mixture().nSpecies(); ++i)
			sum += mp_wrk1[i]*mp_wrk2[i]/mixture().speciesMw(i);

		return (sum*mixture().T()*RU);
	}

        void jacobianRho(double* const p_jacRho)
        {
                std::fill(p_jacRho, p_jacRho + m_ns, 0.);

		mixture().speciesHOverRT(NULL, NULL, NULL, NULL, mp_wrk1, NULL);
                mixture().jacobianRho(mp_wrk3);

                double aux = mixture().T()*RU;

                for(int i = 0; i < m_ns; ++i)
                        for(int j = 0; j < m_ns; ++j)
                                p_jacRho[i] += aux*mp_wrk1[i]*mp_wrk3[i+j*m_ns]/mixture().speciesMw(i);

        }

        void jacobianTTv(double* const p_jacTTv)
        {
		std::fill(p_jacTTv, p_jacTTv + m_nt, 0.);

        	double T=mixture().T();
        	double Tv=mixture().Tv();

		mixture().speciesHOverRT(NULL, NULL, NULL, NULL, mp_wrk1, NULL);
		mixture().netProductionRates(mp_wrk2);
                mixture().speciesCpOverR(T, Tv, T, Tv, Tv, NULL, NULL, NULL, NULL, mp_wrk4);
                mixture().jacobianT(mp_wrk5);
                mixture().jacobianTv(mp_wrk6);

                double aux = T*RU;
                
                for(int i = 0; i < m_ns; ++i) {
                        p_jacTTv[0] += aux*mp_wrk1[i]*mp_wrk5[i]/mixture().speciesMw(i);
                        p_jacTTv[1] += aux*mp_wrk1[i]*mp_wrk6[i]/mixture().speciesMw(i);
                        p_jacTTv[1] += RU*mp_wrk4[i]*mp_wrk2[i]/mixture().speciesMw(i);
                }

        }

private:
        int m_ns;
        int m_nt;
	double* mp_wrk1;
	double* mp_wrk2;
	double* mp_wrk3;
	double* mp_wrk4;
	double* mp_wrk5;
	double* mp_wrk6;
};

// Register the transfer model
Mutation::Utilities::Config::ObjectProvider<
    OmegaCElec, TransferModel> omegaCElec("OmegaCElec");

    } // namespace Transfer
} // namespace Mutation
