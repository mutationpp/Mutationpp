/**
 * @file ChemNonEqTTvStateModel.h
 *
 * @brief Provides ChemNonEqTTvStateModel.
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

#include "StateModel.h"
#include "AutoRegistration.h"

namespace Mutation {
    namespace Thermodynamics {

using namespace Mutation::Transfer;
using namespace Mutation::Utilities::Config;

/**
 * @ingroup statemodels
 *
 * Chemical Non-equilibrium 2-T StateModel. Thermal nonequilibrium with
 * T = Tr and Tv = Tel = Te with finite rate chemistry.
 */
class ChemNonEqTTvStateModel : public StateModel
{
public:
  
     //Mutation::Transfer::TransferModel* p_transfer_model;

    ChemNonEqTTvStateModel(const Thermodynamics& thermo)
        : StateModel(thermo, 2, thermo.nSpecies())
    {
        mp_work1 = new double [thermo.nSpecies()];
        mp_work2 = new double [thermo.nSpecies()];
        mp_work3 = new double [thermo.nSpecies()];
        mp_work4 = new double [thermo.nSpecies()];
    }

    ~ChemNonEqTTvStateModel()
    {
        delete [] mp_work1;
        delete [] mp_work2;
        delete [] mp_work3;
        delete [] mp_work4;
    }

    /**
     * Sets the current state of the mixture.  Variable sets can be the
     * following:
     *   0: {species densities}, {total energy density, vib. energy density}
     *   1: {species densities} and {T, Tv}
     */
    void setState(
        const double* const p_mass, const double* const p_energy,
        const int vars = 0)
    {
        const int ns = m_thermo.nSpecies();

        // Compute the species concentrations which are used through out this
        // method regardless of variable set 
        double conc = 0.0;
        for (int i = 0; i < ns; ++i){
            // Check that species densities are at least positive
            assert(p_mass[i] >= 0.0);
            mp_X[i] = std::max(p_mass[i] / m_thermo.speciesMw(i), 0.0);
            conc += mp_X[i];
        }
        double elec = 0.0;
        if(m_thermo.hasElectrons())
            elec = mp_X[0];

        // Compute the temperatures and make sure the variable set is implemented
        switch (vars) {
        case 0: {
            // First step is to solve one nonlinear equation to get Tv from the
            // vibrational energy density equation
            if (!getTFromRhoE(
                Cvv(m_thermo), Ev(m_thermo), p_energy[1], m_Tv, mp_work1, 0.0)) {
                std::cout << "was trying to compute Tv..." << std::endl;
                std::cout << "rhoe_v = " << p_energy[1] << std::endl;
                for (int i = 0; i < m_thermo.nSpecies(); ++i)
                    std::cout << std::setw(20) << m_thermo.speciesName(i)
                              << std::setw(20) << p_mass[i] << std::endl;
                std::cout << std::endl;
            }

            // Next compute the temperature using the total energy density
            m_Tel = m_Te = m_Tv;

            if (!getTFromRhoE(
                Cv(m_thermo), E(m_thermo, m_Tv), p_energy[0], m_T, mp_work1, 0.0)) {
                std::cout << "was trying to compute T..." << std::endl;
                std::cout << "rhoe = " << p_energy[0] << std::endl;
                for (int i = 0; i < m_thermo.nSpecies(); ++i)
                    std::cout << std::setw(20) << m_thermo.speciesName(i)
                              << std::setw(20) << p_mass[i] << std::endl;
                std::cout << std::endl;
            }
//            getTFromRhoE(
//                Cvv(m_thermo), Ev(m_thermo), p_energy[1], m_Tv, mp_work1, 0.0);
//            m_Tel = m_Te = m_Tv;
//
//            // Next compute the temperature using the total energy density
//            getTFromRhoE(
//                Cv(m_thermo), E(m_thermo, m_Tv), p_energy[0], m_T, mp_work1, 0.0);

            break;
        }
        case 1:
            // Check that temperature is at least positive
            assert(p_energy[0] > 0.0);
            assert(p_energy[1] > 0.0);
            m_T = p_energy[0];
            m_Tv = p_energy[1];
            break;
        default:
            std::cout << "Variable-set " << vars << " not implemented in StateModel!" << std::endl;
            exit(1);
        }

        // Set Tr, Tel, and Te
        m_Tr = m_T;
        m_Tel = m_Te = m_Tv;

        // Compute the pressure and species mole fractions from T and rho_i
        for (int i = 0; i < ns; ++i)        
            mp_X[i] /= conc;
        m_P = conc * RU * (m_T + (m_Tv - m_T) * elec / conc);
    }
    
    void initializeTransferModel(Mutation::Mixture& mix)
    {
        // Heavy particle terms
        addTransferTerm(0, Factory<TransferModel>::create("OmegaVT", mix));
        addTransferTerm(0, Factory<TransferModel>::create("OmegaCV", mix));
        addTransferTerm(0, Factory<TransferModel>::create("OmegaCElec", mix));

        // Terms only included when electrons are present
        if (m_thermo.hasElectrons()) {
            addTransferTerm(0, Factory<TransferModel>::create("OmegaET", mix));
            addTransferTerm(0, Factory<TransferModel>::create("OmegaCE", mix));
            addTransferTerm(0, Factory<TransferModel>::create("OmegaI", mix));
        }
    }
    
    void getTemperatures(double* const p_T) const {
        p_T[0] = T();
        p_T[1] = Tv();
    }
    
    void getEnergiesMass(double* const p_e)
    {    
        int ns = m_thermo.nSpecies();
        m_thermo.speciesHOverRT(mp_work1, mp_work2, NULL, mp_work3, mp_work4, NULL);
        int offset = (m_thermo.hasElectrons() ? 1 : 0);

        for(int i = offset; i < ns; ++i)
            p_e[i] = (mp_work1[i] - 1.0)*m_T*RU/m_thermo.speciesMw(i);

        for(int i = offset; i < ns; ++i)
            p_e[i+ns] = (mp_work3[i] + mp_work4[i])*m_T*RU/m_thermo.speciesMw(i);

        if (m_thermo.hasElectrons()) {
            p_e[0] = (mp_work1[0]*m_T-m_Tv)*RU/m_thermo.speciesMw(0);
            p_e[ns] = (mp_work2[0]*m_T-m_Tv)*RU/m_thermo.speciesMw(0);
        }
    }

    void getEnthalpiesMass(double* const p_h)
    {    
        int ns = m_thermo.nSpecies();
        m_thermo.speciesHOverRT(mp_work1, mp_work2, NULL, mp_work3, mp_work4, NULL);
        
        for(int i = 0; i < ns; ++i)
            p_h[i] = mp_work1[i]*m_T*RU/m_thermo.speciesMw(i);
        for(int i = 0; i < ns; ++i)
            p_h[i+ns] = (mp_work3[i] + mp_work4[i])*m_T*RU/m_thermo.speciesMw(i);
        if(m_thermo.hasElectrons())
            p_h[ns] = mp_work2[0]*m_T*RU/m_thermo.speciesMw(0);
    }
    
    void getCpsMass(double* const p_Cp)
    {       
        int ns = m_thermo.nSpecies();
        int offset = (m_thermo.hasElectrons() ? 1 : 0);
        m_thermo.speciesCpOverR(
			m_T, m_Te, m_Tr, m_Tv, m_Tel, NULL, mp_work1, mp_work2, mp_work3, mp_work4);

        for(int i = offset; i < ns; ++i)
            p_Cp[i] = (mp_work1[i]+mp_work2[i])*RU/m_thermo.speciesMw(i);
        for(int i = offset; i < ns; ++i)
            p_Cp[i+ns] = (mp_work3[i]+mp_work4[i])*RU/m_thermo.speciesMw(i);
        if(m_thermo.hasElectrons()) {
            p_Cp[0] = 0.0;
            p_Cp[ns] = mp_work1[0]*RU/m_thermo.speciesMw(0);
        }
    }

    void getCvsMass(double* const p_Cv)
    {       
        int ns = m_thermo.nSpecies();
        int offset = (m_thermo.hasElectrons() ? 1 : 0);
        m_thermo.speciesCpOverR(
			m_T, m_Te, m_Tr, m_Tv, m_Tel, NULL, mp_work1, mp_work2, mp_work3, mp_work4);

        for(int i = offset; i < ns; ++i)
            p_Cv[i] = (mp_work1[i]+mp_work2[i]-1.0)*RU/m_thermo.speciesMw(i);
        for(int i = offset; i < ns; ++i)
            p_Cv[i+ns] = (mp_work3[i]+mp_work4[i])*RU/m_thermo.speciesMw(i);
        if(m_thermo.hasElectrons()) {
            p_Cv[0] = 0.0;
            p_Cv[ns] = (mp_work1[0]-1.0)*RU/m_thermo.speciesMw(0);
        }
    }

    void getTagModes(int* const p_tag)
    {
     	p_tag[0] = 1; p_tag[5] = 0; // Heavy translation
     	p_tag[1] = 0; p_tag[6] = 1; // Electron translation
     	p_tag[2] = 1; p_tag[7] = 0; // Rotation excitation
     	p_tag[3] = 0; p_tag[8] = 1; // Vibration excitation
     	p_tag[4] = 0; p_tag[9] = 1; // Electronic excitation
    }
    
    private:

    /**
     * Small helper class which provides wrapper to
     * Thermodynamics::speciesCpOverR() to be used with getTFromRhoE().
     */
    
    class Cv {
    public:
        Cv(const Thermodynamics& t) : thermo(t) {
             p_cvt = new double [thermo.nSpecies()]; 
             p_cvr = new double [thermo.nSpecies()]; 
        }
        ~Cv() {
             delete [] p_cvt;
             delete [] p_cvr;
        }
        void operator () (double T, double* const cv) const {
          thermo.speciesCpOverR(T, T, T, T, T, NULL, p_cvt, p_cvr, NULL, NULL);
          int n_species = thermo.nSpecies();
          for(int iCv = 0; iCv < n_species; ++iCv) 
              cv[iCv] = p_cvt[iCv]+p_cvr[iCv]-1.0;
          if(thermo.hasElectrons())
              cv[0] = 0.0;
        }

    private:
        const Thermodynamics& thermo;
        double* p_cvt;
        double* p_cvr;
    }; // class Cp
    
    class Cvv {
    public:
        Cvv(const Thermodynamics& t) : thermo(t) {
             p_cvt = new double [thermo.nSpecies()]; 
             p_cvv = new double [thermo.nSpecies()];
             p_cvel = new double [thermo.nSpecies()];
        }
        ~Cvv() {
             delete [] p_cvt;
             delete [] p_cvv;
             delete [] p_cvel;
        }
        void operator () (double T, double* const cv) const {
            thermo.speciesCpOverR(T, T, T, T, T, NULL, p_cvt, NULL, p_cvv, p_cvel);
            for(int iCvv = 0; iCvv < thermo.nSpecies(); ++iCvv)
                cv[iCvv] = p_cvv[iCvv] + p_cvel[iCvv];
            if(thermo.hasElectrons())
                cv[0] = p_cvt[0]-1.0;
        }
    private:
        const Thermodynamics& thermo;
	double* p_cvt;
	double* p_cvv;
	double* p_cvel;
    }; // class Cpv

    /**
     * Small helper classes which provide wrapper to
     * Thermodynamics::speciesHOverRT() to be used with getTFromRhoE().
     */
    
    class E {
    public:
        E(const Thermodynamics& t, const double& T_vibration) : thermo(t) {
            Tv = T_vibration;
        }
        void operator () (double T, double* const e) const {
          thermo.speciesHOverRT(T, Tv, T, Tv, Tv, e, NULL, NULL, NULL, NULL, NULL);

          int offset = (thermo.hasElectrons() ? 1 : 0);
          for(int i = offset; i < thermo.nSpecies(); ++i)
              e[i] -= 1.0;
          if(thermo.hasElectrons())
              e[0] -= Tv/T;

        }
    private:
        const Thermodynamics& thermo;
        double Tv;
    }; // class E
    
    class Ev {
    public:
        Ev(const Thermodynamics& t) : thermo(t) { 
             p_et  = new double [thermo.nSpecies()]; 
             p_ev  = new double [thermo.nSpecies()]; 
             p_eel = new double [thermo.nSpecies()]; 
        }
        ~Ev() {
             delete [] p_et;
             delete [] p_ev;
             delete [] p_eel;
        }
        void operator () (double T, double* const e) const {
            thermo.speciesHOverRT(T, T, T, T, T, NULL, p_et, NULL, p_ev, p_eel, NULL);
            for(int iev = 0; iev < thermo.nSpecies(); ++iev){
                e[iev] = p_ev[iev] + p_eel[iev];
            if(thermo.hasElectrons())
                e[0] = p_et[0]-1.0;
            }
        }
    private:
        const Thermodynamics& thermo;
        double* p_et;
        double* p_ev;
        double* p_eel;
    }; // class Ev

private:
    double* mp_work1;
    double* mp_work2;
    double* mp_work3;
    double* mp_work4;
    
}; // class ChemNonEqStateModel

// Register the state model
Utilities::Config::ObjectProvider<
    ChemNonEqTTvStateModel, StateModel> chem_non_eq_TTv("ChemNonEqTTv");

    } // namespace Thermodynamics
} // namespace Mutation
