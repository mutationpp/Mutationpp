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
        m_ns = mixture().nSpecies();
        m_nr = mixture().nReactions();
        m_nt = mixture().nEnergyEqns();
        mp_hf = new double [m_ns];
        mp_h = new double [m_ns];
        mp_cp = new double [m_ns];
        mp_rate = new double [m_nr];
        mp_delta = new double [m_nr];
        mp_delta1 = new double [m_nr];
        mp_dropRho = new double [m_nr*m_ns];
        mp_dropT = new double [m_nr];
        mp_dropTv = new double [m_nr];
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
        delete [] mp_cp;
        delete [] mp_rate;
        delete [] mp_delta;
        delete [] mp_delta1;
        delete [] mp_dropRho;
        delete [] mp_dropT;
        delete [] mp_dropTv;
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
        mixture().speciesHOverRT(mp_h, NULL, NULL, NULL, NULL, mp_hf);

        // Get reaction enthalpies
        std::fill(mp_delta, mp_delta+m_nr, 0.0);
        mixture().getReactionDelta(mp_hf,mp_delta);
        
        // Get molar rates of progress
        mixture().netRatesOfProgress(mp_rate);

        double src = 0.0;
        int j;
        for (int i = 0; i < m_rId.size(); ++i) {
            j = m_rId[i];
            src += mp_delta[j]*mp_rate[j];
        }

        return (-src*RU*mixture().T());

    }

    /**
      * Computes the partial density-based derivatives for Electron-Impact reactions 
      * heat Generation.
      *
      * \f[ \frac{\partial \Omega^{I}}{\partial \rhoi_{i}} = 
      *       - \sum_{r \in \mathcal{R}} \Delta h_r \frac{\partial \xi_r}{\partial \rho_{i}} \f]
      *
      * \f[ \mathcal{R} \f] denotes the set of electron impact reactions.
      * \f[ \Delta h_r \f] is the reaction enthalpy \f[ [J/mol] \f].
      * \f[ \frac{\partial \xi_r}{\partial \rho_{i}} \f] is the density-derivative of molar rate of progress.
      *
      */
    void jacobianRho(double* const p_jacRho)
    {

	std::fill(p_jacRho, p_jacRho + m_ns, 0.);

        // Get Formation enthalpy
        mixture().speciesHOverRT(mp_h, NULL, NULL, NULL, NULL, mp_hf);

        // Get reaction enthalpies
        std::fill(mp_delta, mp_delta+m_nr, 0.0);
        mixture().getReactionDelta(mp_hf,mp_delta);
        
        // Get molar rates of progress
        mixture().netRateProgressRhoJacobian(mp_dropRho);

        int j;
        for (int i = 0; i < m_rId.size(); ++i) {
            j = m_rId[i];
            for (int k = 0; k < m_ns; ++k) {
            	p_jacRho[k] += mp_delta[j]*mp_dropRho[k+j*m_ns];
            }
        }

        for (int i = 0; i < m_ns; ++i)
            p_jacRho[i] *= -RU*mixture().T();

    }

    /**
      * Computes the partial temperature-based derivatives for Electron-Impact reactions 
      * heat Generation.
      *
      * \f[ \frac{\partial \Omega^{I}}{\partial T} = 
      *     - \sum_{r \in \mathcal{R}} \Delta h_r \frac{\partial \xi_r}{\partial T} 
      *     - \sum_(r \in \mathcal{R}} \frac{\partial \Delta h_r}{\partial T} \xi_r \f]
      *
      * \f[ \mathcal{R} \f] denotes the set of electron impact reactions.
      * \f[ \Delta h_r \f] is the reaction enthalpy \f[ [J/mol] \f].
      * \f[ \frac{\partial \Delta h_r}{\partial T} \f] is the temperature-derivative of reaction enthalpy.
      * \f[ \xi_r \f] is the molar rate of progress \f[ [mol/(m^3\cdot s)] \f].
      * \f[ \frac{\partial \xi_r}{\partial T} \f] is the temperature-derivative of molar rate of progress. 
      *
      */
    void jacobianTTv(double* const p_jacTTv)
    {

	std::fill(p_jacTTv, p_jacTTv + m_nt, 0.);

        // Get Formation enthalpy
        mixture().speciesHOverRT(mp_h, NULL, NULL, NULL, NULL, mp_hf);

        // Get derivative of formation enthalpy
        double Tss = mixture().standardStateT();
        mixture().speciesCpOverR(Tss,mp_cp);
        for (int i = 0; i < m_ns; ++i)
            mp_cp[i] *= -1.;

        // Get reaction enthalpies
        std::fill(mp_delta, mp_delta+m_nr, 0.0);
        mixture().getReactionDelta(mp_hf,mp_delta);
        
        // Get reaction specific heats
        std::fill(mp_delta1, mp_delta1+m_nr, 0.0);
        mixture().getReactionDelta(mp_cp,mp_delta1);

        // Get molar rates of progress
        mixture().netRatesOfProgress(mp_rate);

        // Get molar rates of progress
        mixture().netRateProgressTJacobian(mp_dropT);
        mixture().netRateProgressTvJacobian(mp_dropTv);

        int j;
        for (int i = 0; i < m_rId.size(); ++i) {
            j = m_rId[i];
            p_jacTTv[0] += -RU*mixture().T()*mp_delta[j]*mp_dropT[j];
            p_jacTTv[0] += -RU*mp_delta1[j]*mp_rate[j];
            p_jacTTv[1] += -RU*mixture().T()*mp_delta[j]*mp_dropTv[j];
        }

    }

private:
    int m_ns;
    int m_nr;
    int m_nt;
    std::vector<int> m_rId;
    double* mp_hf;
    double* mp_h;
    double* mp_cp;
    double* mp_rate;
    double* mp_delta;
    double* mp_delta1;
    double* mp_dropRho;
    double* mp_dropT;
    double* mp_dropTv;
};
  
// Register the transfer model
Mutation::Utilities::Config::ObjectProvider<
    OmegaI, TransferModel> omegaI("OmegaI");

    } // namespace Transfer
} // namespace Mutation
