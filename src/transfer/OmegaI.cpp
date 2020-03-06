/**
 * @file OmegaI.cpp
 *
 * @brief Implementation of OmegaI.
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
#include <vector>


namespace Mutation {
    namespace Transfer {
      
class OmegaI : public TransferModel
{

public:
    OmegaI(Mutation::Mixture& mix)
        : TransferModel(mix)
    {
        m_ns = m_mixture.nSpecies();
        m_nr = m_mixture.nReactions();
        mp_hf = new double [m_ns];
        mp_h = new double [m_ns];
        mp_rate = new double [m_nr];
        mp_delta = new double [m_nr];
        for(int i=0; i<m_nr; ++i) {
            if ( (mix.reactions()[i].type() == Kinetics::IONIZATION_E)
                || (mix.reactions()[i].type() == Kinetics::DISSOCIATION_E)
                || (mix.reactions()[i].type() == Kinetics::ELECTRONIC_DETACHMENT_E)
                || (mix.reactions()[i].type() == Kinetics::EXCITATION_E)
                || (mix.reactions()[i].type() == Kinetics::ION_RECOMBINATION_E)
                || (mix.reactions()[i].type() == Kinetics::RECOMBINATION_E)
                || (mix.reactions()[i].type() == Kinetics::ELECTRONIC_ATTACHMENT_E) )
                m_rId.push_back(i);
        }
    };

    ~OmegaI() {
        delete [] mp_hf;
        delete [] mp_h;
        delete [] mp_rate;
        delete [] mp_delta;
    };

    /**
      * Computes the Electron-Impact reactions heat Generation in \f$ [J/(m^3\cdot s)] \f$
      * which acts as a sink to the free electron energy equation. Considers both reactions
      * such as electron impact ionization and dissociation
      *
      * \f[ \Omega^{I} = - \sum_{r \in \mathcal{R}} \Delta h_r \xi_r \f]
      *
      * \f[ \mathcal{R} \f] denotes the set of electron impact reactions.
      * \f[ \Delta h_r \f] is the reaction enthalpy \f[ [J/mol] \f]
      * \f[ \xi_r \f] is the molar rate of progress \f[ [mol/(m^3\cdot s)] \f]
      *
      * @todo consider the ionization energy for N and O from an already
      * excited state as suggested by Johnston (PhD thesis)
      */
    double source()
    {
        // Get Formation enthalpy
        m_mixture.speciesHOverRT(mp_h, NULL, NULL, NULL, NULL, mp_hf);

        // Get reaction enthalpies
        std::fill(mp_delta, mp_delta+m_nr, 0.0);
        m_mixture.getReactionDelta(mp_hf,mp_delta);

        // Get molar rates of progress
        m_mixture.netRatesOfProgress(mp_rate);

        double src = 0.0;
        int j;
        for (int i = 0; i < m_rId.size(); ++i) {
            j = m_rId[i];
            src += mp_delta[j]*mp_rate[j];
        }

        return (-src*RU*m_mixture.T());

    }

private:
    int m_ns;
    int m_nr;
    std::vector<int> m_rId;
    double* mp_hf;
    double* mp_h;
    double* mp_rate;
    double* mp_delta;
};
  
// Register the transfer model
Mutation::Utilities::Config::ObjectProvider<
    OmegaI, TransferModel> omegaI("OmegaI");

    } // namespace Transfer
} // namespace Mutation
