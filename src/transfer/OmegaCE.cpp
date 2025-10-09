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
        mp_wrk2 = new double [mix.nSpecies()*mix.nSpecies()];
        mp_wrk3 = new double [mix.nSpecies()];
        mp_wrk4 = new double [mix.nSpecies()];
    }

    ~OmegaCE()
    {
        delete [] mp_wrk1;
        delete [] mp_wrk2;
        delete [] mp_wrk3;
        delete [] mp_wrk4;
    }

    double source()
    {
        static const double cv = 1.5*RU/mixture().speciesMw(0);
        mixture().netProductionRates(mp_wrk1);
        return mp_wrk1[0]*cv*mixture().Te();
    }
     
    void jacobianRho(double* const p_jacRho)
    {
        const size_t ns=mixture().nSpecies();
        std::fill(p_jacRho, p_jacRho + ns, 0.);
        static const double cv = 1.5*RU/mixture().speciesMw(0);
        mixture().jacobianRho(mp_wrk2);
        for(int i = 0; i < ns; ++i)
            p_jacRho[i] = mp_wrk2[i]*cv*mixture().Te();
    }
        
    void jacobianTTv(double* const p_jacTTv)
    {
        const size_t nt=mixture().nEnergyEqns();
        std::fill(p_jacTTv, p_jacTTv + nt, 0.);
        static const double cv = 1.5*RU/mixture().speciesMw(0);
        mixture().netProductionRates(mp_wrk1);
        mixture().jacobianT(mp_wrk3);
        mixture().jacobianTv(mp_wrk4);
        p_jacTTv[0] = cv*mixture().Te()*mp_wrk3[0];
        p_jacTTv[1] = mp_wrk1[0]*cv;
        p_jacTTv[1] += cv*mixture().Te()*mp_wrk4[0];
    }

private:
    double* mp_wrk1;
    double* mp_wrk2;
    double* mp_wrk3;
    double* mp_wrk4;
};

// Register the transfer model
Mutation::Utilities::Config::ObjectProvider<
    OmegaCE, TransferModel> omegaCE("OmegaCE");
      
    } // namespace Transfer
} // namespace Mutation
