/**
 * @file ChemNonEqT3TvStateModel.h
 *
 * @brief Provides ChemNonEqT3TvStateModel.
 */

/*
 * Copyright 2014-2018 von Karman Institute for Fluid Dynamics (VKI)
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
#include <eigen3/Eigen/Dense>
#include "Thermodynamics.h"
#include "Transport.h"

using namespace Eigen;

namespace Mutation {
    namespace Thermodynamics {

using namespace Mutation::Transfer;
using namespace Mutation::Utilities::Config;

/**
 * @ingroup statemodels
 *
 * Chemical Non-equilibrium 4-T StateModel. Thermal nonequilibrium with
 * T = Tr and T1 = Tv_NO, T2 = Tv_N2, T3 = Tv_O2 with finite rate chemistry.
 */
class ChemNonEqT3TvStateModel : public StateModel
{
public:
  
    ChemNonEqT3TvStateModel(const Thermodynamics& thermo)
        : StateModel(thermo, 4, thermo.nSpecies())
    {
        mp_work1 = new double [thermo.nSpecies()];
        mp_work2 = new double [thermo.nSpecies()];
        mp_work3 = new double [thermo.nSpecies()];
        mp_work4 = new double [thermo.nSpecies()];

	mp_Tvs = new double [3];
    }

    ~ChemNonEqT3TvStateModel()
    {
        delete [] mp_work1;
        delete [] mp_work2;
        delete [] mp_work3;
        delete [] mp_work4;

	delete [] mp_Tvs;
    }

    /**
     * Sets the current state of the mixture.  Variable sets can be the
     * following:
     *   0: {species densities}, {total energy density, vib. energy density}
     *   1: {species densities} and {T, Tv[3]}
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

            solveEnergies(p_mass, p_energy);

            break;
        }
        case 1:
            // Check that temperature is at least positive
            assert(p_energy[0] > 0.0);
            assert(p_energy[1] > 0.0);
            assert(p_energy[2] > 0.0);
            assert(p_energy[3] > 0.0);
            m_T      = p_energy[0];
            mp_Tvs[0] = p_energy[1];
            mp_Tvs[1] = p_energy[2];
            mp_Tvs[2] = p_energy[3];
	    std::cout << m_T       << " " 
		      << mp_Tvs[0] << " " 
		      << mp_Tvs[1] << " " 
		      << mp_Tvs[2] << std::endl;  
		   
            break;
        default:
            throw InvalidInputError("variable set", vars)
                << "This variable-set is not implemented in ChemNonEqStateModel"
                << ". Possible variable-sets are:\n"
                << "  0: (species densities, static energy density)\n"
                << "  1: (species densities, T and Tv[3])";
        }

        // Set Tr, Tel, and Te
        m_Tr = m_T;
        m_Tel = m_Te = m_T; // Not used in air5! 

        // Compute the pressure and species mole fractions from T and rho_i
        for (int i = 0; i < ns; ++i)        
            mp_X[i] /= conc;
        m_P = conc * RU * (m_T + /*for air 5 the terms below is null*/
			  (mp_Tvs[0] + mp_Tvs[1] + mp_Tvs[2] - m_T) * elec / conc);
    }
    /**
     *  Add proper
     *   0: {species densities}, {total energy density, vib. energy density}
     *   1: {species densities} and {T, Tv[3]}
     */
    void initializeTransferModel(Mutation::Mixture& mix)
    {
        typedef Factory<TransferModel> Factory;

        try { // TODO: models below need to be written
            // Heavy particle terms
            //addTransferTerm(0, Factory::create("R_VT", mix));
            //addTransferTerm(0, Factory::create("R_VV", mix));

	    // Not sure this can be here ...
	    // ... maybe should be in kinetics but here is clearer
            //addTransferTerm(0, Factory::create("R_diss_rec", mix));

            //addTransferTerm(0, Factory::create("OmegaVT", mix));
            //addTransferTerm(0, Factory::create("OmegaCV", mix));
            //addTransferTerm(0, Factory::create("OmegaCElec", mix));

            //// Terms only included when electrons are present
            //if (m_thermo.hasElectrons()) {
            //    addTransferTerm(0, Factory::create("OmegaET", mix));
            //    addTransferTerm(0, Factory::create("OmegaCE", mix));
            //    addTransferTerm(0, Factory::create("OmegaI", mix));
            //}
        } catch (Error& e) {
            e << "\nWas trying to load a energy transfer model.";
            throw;
        }
    }
    
    void getTemperatures(double* const p_T) const {
        p_T[0] = T();
	// Remember that this is a 1+3-temperature model
	for (int i=0; i<3; ++i)
            p_T[i+1] = Tvs(i);
    }
    
    void getEnergiesMass(double* const p_e)
    {    
        int ns = m_thermo.nSpecies();

	// This call raises some issues because of the function signature
	// At lower level it calls hV which use one entry T but loops over
	// all molecular species.
	// Possible solution (not dirty) is to pass Tv[] vector and the 
	// index loop will correspond to the molecular one.
	// It requires a neverending sequence of changes ... uff
	// Maybe better to call a function with different signature!
	m_thermo.speciesHOverRT(
            m_T, m_Tv, m_T, mp_Tvs, m_Tv, mp_work1, mp_work2, NULL, mp_work3, 
	    mp_work4, NULL);
        //m_thermo.speciesHOverRT(
        //    m_T, mp_Tvs, mp_work1, mp_work2, NULL, mp_work3, NULL, NULL);

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

	// This call raises some issues because of the function signature
	// At lower level it calls hV which use one entry T but loops over
	// all molecular species.
	// Possible solution (not dirty) is to pass Tv[] vector and the 
	// index loop will correspond to the molecular one.
	// It requires a neverending sequence of changes ... uff
	// Maybe better to call a function with different signature!
	m_thermo.speciesHOverRT(
            m_T, m_Tv, m_T, m_Tv, m_Tv, mp_work1, mp_work2, NULL, mp_work3, 
	    mp_work4, NULL);
        //m_thermo.speciesHOverRT(
        //    m_T, mp_Tvs, mp_work1, mp_work2, NULL, mp_work3, NULL, NULL);
        
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

	// This call raises some issues because of the function signature
	// At lower level it calls cpV which use one entry T but loops over
	// all molecular species.
	// Possible solution (not dirty) is to pass Tv[] vector and the 
	// index loop will correspond to the molecular one.
	// It requires a neverending sequence of changes ... uff
	// Maybe better to call a function with different signature!
	m_thermo.speciesHOverRT(
            m_T, m_Tv, m_T, m_Tv, m_Tv, mp_work1, mp_work2, NULL, mp_work3, 
	    mp_work4, NULL);
        //m_thermo.speciesHOverRT(
        //    m_T, mp_Tvs, mp_work1, mp_work2, NULL, mp_work3, NULL, NULL);

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

	// This call raises some issues because of the function signature
	// At lower level it calls cpV which use one entry T but loops over
	// all molecular species.
	// Possible solution (not dirty) is to pass Tv[] vector and the 
	// index loop will correspond to the molecular one.
	// It requires a neverending sequence of changes ... uff
	// Maybe better to call a function with different signature!
	m_thermo.speciesHOverRT(
            m_T, m_Tv, m_T, m_Tv, m_Tv, mp_work1, mp_work2, NULL, mp_work3, 
	    mp_work4, NULL);
        //m_thermo.speciesHOverRT(
        //    m_T, mp_Tvs, mp_work1, mp_work2, NULL, mp_work3, NULL, NULL);

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

    /**
    * @brief Get temperatures from total energy by solving a non-linear system
    *
    * @param p_rhoi - vector of the partial densities
    * @param p_rhoe - total and internal mass energy
    */
    void solveEnergies(const double* const p_rhoi, const double* const p_rhoe)
    {
        const double atol = 1.0e-12;
        const double rtol = 1.0e-12;
        const int    imax = 100;

        Map<const VectorXd> rhoi(p_rhoi, m_thermo.nSpecies());
        const double density = rhoi.sum();

        static VectorXd yi; yi = rhoi / density;
        const Vector2d emix = Map<const Vector2d>(p_rhoe) / density;

        static Matrix<double, Dynamic, Dynamic, RowMajor> ei;
        static Matrix<double, Dynamic, Dynamic, RowMajor> ci;
        ei.resize(2, m_thermo.nSpecies());
        ci.resize(2, m_thermo.nSpecies());

        getEnergiesMass(ei.data());

        Vector2d f, e, cv;
        e = ei*yi;
        f = e - emix;

        int i;
        for (i = 0; (f.norm() > rtol*emix.norm() + atol) && (i < imax); ++i) {

            // Update temperatures
            getCvsMass(ci.data());
            cv = ci*yi;

            mp_Tvs[0] = std::max(mp_Tvs[0] - f[1]/cv[1], 0.1*mp_Tvs[0]);
            mp_Tvs[1] = std::max(mp_Tvs[1] - f[1]/cv[1], 0.1*mp_Tvs[1]);
            mp_Tvs[2] = std::max(mp_Tvs[2] - f[1]/cv[1], 0.1*mp_Tvs[2]);
            m_T  = std::max(m_T + (f[1]-f[0])/cv[0], 0.1*m_T);

            // Update function evaluation
            getEnergiesMass(ei.data());
            e = ei*yi;
            f = e - emix;
        }

        if (i == imax)
            std::cout << "Warning, didn't converge temperatures: f = " 
		      << f.norm() << std::endl;
    }

private:
    double* mp_work1;
    double* mp_work2;
    double* mp_work3;
    double* mp_work4;

    double* mp_Tvs;
    
}; // class ChemNonEqStateModel

// Register the state model
Utilities::Config::ObjectProvider<
    ChemNonEqT3TvStateModel, StateModel> chem_non_eq_T3Tv("ChemNonEqT3Tv");

    } // namespace Thermodynamics
} // namespace Mutation
