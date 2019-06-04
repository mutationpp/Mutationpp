/**
 * @file R_vibr_VV.cpp
 *
 * @brief Implementation of R_vibr_VV.
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

#include "Mixture.h"
#include "TransferModel.h"
#include "Numerics.h"

#include <cmath>
#include <fstream>
#include <vector>

using namespace Mutation;

namespace Mutation {
    namespace Transfer {

/**
 * Represents a coupling between vibrational energy modes of different molecules
 */
class R_vibr_VV : public TransferModel
{
public:

    R_vibr_VV(Mixture& mix)
        : TransferModel(mix)
    {
        m_ns              = m_mixture.nSpecies();
        m_transfer_offset = m_mixture.hasElectrons() ? 1 : 0;

        mp_Mw = new double [m_ns];
        for(int i = 0; i < m_ns; ++i)
            mp_Mw[i] = m_mixture.speciesMw(i);

	mp_rVV = new double [m_mixture.nEnergyEqns()-1];
    }

    virtual ~R_vibr_VV()
    {
        delete [] mp_Mw;
	delete [] mp_rVV;
    }

    double source() {};
    void rVT(double* const mp_rVT) {};
    void R_w_22(double* const p_r_w_22) {};
    void R_w_23(double* const p_r_w_23) {};

    void rVV(double* const mp_rVV)
    {
        const double * p_X = m_mixture.X();
        const double   n   = m_mixture.numberDensity();
        std::cout << "number density: " << n << std::endl;
        double T = m_mixture.T();
        std::cout << T << std::endl;
        double Tv[3];
        for (int i=0; i<3; ++i) {
            Tv[i] = m_mixture.Tvs(i); // NO, N2, O2
            std::cout << Tv[i] << std::endl;
        }
        double KBT = KB * T;
    
        // Number of vibrational levels
	// for harmonic oscillator model
        int lNO = 28;
        int lN2 = 33;
        int lO2 = 26;
    
        // Vibrational energy arrays
        static double ve[87];
        static double veNO[28];
        static double veN2[33];
        static double veO2[26];
    
        veNO[0] = 0.00000000e+00;
        veNO[1] = 3.78297694e-20;
        veNO[2] = 7.56595389e-20;
        veNO[3] = 1.13489308e-19;
        veNO[4] = 1.51319078e-19;
        veNO[5] = 1.89148847e-19;
        veNO[6] = 2.26978617e-19;
        veNO[7] = 2.64808386e-19;
        veNO[8] = 3.02638156e-19;
        veNO[9] = 3.40467925e-19;
        veNO[10] = 3.78297694e-19;
        veNO[11] = 4.16127464e-19;
        veNO[12] = 4.53957233e-19;
        veNO[13] = 4.91787003e-19;
        veNO[14] = 5.29616772e-19;
        veNO[15] = 5.67446542e-19;
        veNO[16] = 6.05276311e-19;
        veNO[17] = 6.43106081e-19;
        veNO[18] = 6.80935850e-19;
        veNO[19] = 7.18765620e-19;
        veNO[20] = 7.56595389e-19;
        veNO[21] = 7.94425158e-19;
        veNO[22] = 8.32254928e-19;
        veNO[23] = 8.70084697e-19;
        veNO[24] = 9.07914467e-19;
        veNO[25] = 9.45744236e-19;
        veNO[26] = 9.83574006e-19;
        veNO[27] = 1.02140378e-18;
    
        veN2[0] = 0.00000000e+00;
        veN2[1] = 4.68454043e-20;
        veN2[2] = 9.36908086e-20;
        veN2[3] = 1.40536213e-19;
        veN2[4] = 1.87381617e-19;
        veN2[5] = 2.34227021e-19;
        veN2[6] = 2.81072426e-19;
        veN2[7] = 3.27917830e-19;
        veN2[8] = 3.74763234e-19;
        veN2[9] = 4.21608639e-19;
        veN2[10] = 4.68454043e-19; 
        veN2[11] = 5.15299447e-19;
        veN2[12] = 5.62144851e-19;
        veN2[13] = 6.08990256e-19;
        veN2[14] = 6.55835660e-19;
        veN2[15] = 7.02681064e-19;
        veN2[16] = 7.49526469e-19;
        veN2[17] = 7.96371873e-19;
        veN2[18] = 8.43217277e-19;
        veN2[19] = 8.90062681e-19;
        veN2[20] = 9.36908086e-19;
        veN2[21] = 9.83753490e-19;
        veN2[22] = 1.03059889e-18;
        veN2[23] = 1.07744430e-18;
        veN2[24] = 1.12428970e-18;
        veN2[25] = 1.17113511e-18;
        veN2[26] = 1.21798051e-18;
        veN2[27] = 1.26482592e-18;
        veN2[28] = 1.31167132e-18;
        veN2[29] = 1.35851672e-18;
        veN2[30] = 1.40536213e-18;
        veN2[31] = 1.45220753e-18;
        veN2[32] = 1.49905294e-18;
    
        veO2[0] = 0.00000000e+00;
        veO2[1] = 3.13821409e-20; 
        veO2[2] = 6.27642817e-20;
        veO2[3] = 9.41464226e-20;
        veO2[4] = 1.25528563e-19;
        veO2[5] = 1.56910704e-19;
        veO2[6] = 1.88292845e-19;
        veO2[7] = 2.19674986e-19;
        veO2[8] = 2.51057127e-19;
        veO2[9] = 2.82439268e-19;
        veO2[10] = 3.13821409e-19; 
        veO2[11] = 3.45203549e-19;
        veO2[12] = 3.76585690e-19;
        veO2[13] = 4.07967831e-19;
        veO2[14] = 4.39349972e-19;
        veO2[15] = 4.70732113e-19;
        veO2[16] = 5.02114254e-19;
        veO2[17] = 5.33496395e-19;
        veO2[18] = 5.64878535e-19;
        veO2[19] = 5.96260676e-19;
        veO2[20] = 6.27642817e-19;
        veO2[21] = 6.59024958e-19;
        veO2[22] = 6.90407099e-19;
        veO2[23] = 7.21789240e-19;
        veO2[24] = 7.53171381e-19;
        veO2[25] = 7.84553521e-19;
    
        // Dissociation energy
        const double de[m_mixture.nMolecules()] = {1.0409017238088574e-18,
    	1.5636156480913654e-18, 8.196091362099268e-19};
        for (int i=0; i<m_mixture.nMolecules(); ++i) {
    	std::cout << "de: " << de[i] << std::endl;
        }
    
        // Partial number densities
        double nd[m_mixture.nSpecies()] = {};
        for (int i=0; i<m_mixture.nSpecies(); ++i) {
            nd[i] = p_X[i] * n;
    	std::cout << "nd: " << nd[i] << std::endl;
        }
    
        // Calculate the non-equilibrium partition functions
        double z_vibr_T_Tv_NO = 0.;
        for (int i=0; i<lNO; ++i) {
    
    	// harmonic oscillator
            z_vibr_T_Tv_NO += exp(-veNO[i]/(KB*Tv[0]));
    
    	// anharmonic oscillator
            //z_vibr_T_Tv_NO += exp(-(veNO[i]-i*veNO[0])/(KB*T)-
    	//		                i*veNO[0]/(KB*Tv[0]));
        }

        double z_vibr_T_Tv_N2 = 0.;
        for (int i=0; i<lN2; ++i) {
            z_vibr_T_Tv_N2 += exp(-veN2[i]/(KB*Tv[1]));
            //z_vibr_T_Tv_N2 += exp(-(veN2[i]-i*veN2[0])/(KB*T)-
    	//		                i*veN2[0]/(KB*Tv[1]));
        }

        double z_vibr_T_Tv_O2 = 0.;
        for (int i=0; i<lO2; ++i) {
            z_vibr_T_Tv_O2 += exp(-veO2[i]/(KB*Tv[2]));
            //z_vibr_T_Tv_O2 += exp(-(veO2[i]-i*veO2[0])/(KB*T)-
    	//		                i*veO2[0]/(KB*Tv[2]));
        }

        std::cout << "zNO: " << z_vibr_T_Tv_NO << " " << 
    	         "zN2: " << z_vibr_T_Tv_N2 << " " << 
    		 "zO2: " << z_vibr_T_Tv_O2 << std::endl;
    
        // Compute the vibrational populations
        // according to Treanor-Marrone or Boltzmann
        double nNO[lNO];
        double nN2[lN2];
        double nO2[lO2];
    
        std::cout << " Vibrational population ... " << std::endl;
        for (int i=0; i<lNO; i++) {
            nNO[i] = nd[2] * exp(-veNO[i] / (KB * Tv[0])) / z_vibr_T_Tv_NO;
            //nNO[i] = nd[2] * exp(-(veNO[i] - i * veNO[0]) / KBT 
    	    //                 - i * veNO[0] / (KB * Tv[0])) / z_vibr_T_Tv_NO;
            std::cout << i << " " << nNO[i] << std::endl;
        }

        for (int i=0; i<lN2; i++) {
            nN2[i] = nd[3] * exp(-veN2[i] / (KB * Tv[1])) / z_vibr_T_Tv_N2;
            //nN2[i] = nd[3] * exp(-(veN2[i] - i * veNO[0]) / KBT 
    	    //                 - i * veN2[0] / (KB * Tv[1])) / z_vibr_T_Tv_N2;
            std::cout << i << " " << nN2[i] << std::endl;
        }

        for (int i=0; i<lO2; i++) {
            nO2[i] = nd[4] * exp(-veO2[i] / (KB * Tv[2])) / z_vibr_T_Tv_O2;
            //nO2[i] = nd[4] * exp(-(veO2[i] - i * veO2[0]) / KBT 
    	    //                 - i * veO2[0] / (KB * Tv[2])) / z_vibr_T_Tv_O2;
            std::cout << i << " " << nO2[i] << std::endl;
        }
    
	double deps_n2[lN2] = {};
	for (int i=0; i<lN2-1; ++i) {
	    deps_n2[i] = veN2[i] - veN2[i+1];
	}

	double deps_o2[lO2] = {};
	for (int i=0; i<lO2-1; ++i) {
	    deps_o2[i] = veO2[i] - veO2[i+1];
	}

	double deps_no[lNO] = {};
	for (int i=0; i<lNO-1; ++i) {
	    deps_no[i] = veNO[i] - veNO[i+1];
	}

        // Compute STS vibrational exchange reaction rate coefficients.
	double kvvs_down_n2_o2[lN2][lO2] = {};
        double kvvs_down_n2_no[lN2][lNO] = {};
	double kvvs_down_o2_n2[lO2][lN2] = {}; 
	double kvvs_down_o2_no[lO2][lNO] = {};
	double kvvs_down_no_n2[lNO][lN2] = {}; 
	double kvvs_down_no_o2[lNO][lO2] = {};

	double kvvs_up_n2_o2[lN2][lO2] = {}; 
	double kvvs_up_n2_no[lN2][lNO] = {};
	double kvvs_up_o2_n2[lO2][lN2] = {}; 
	double kvvs_up_o2_no[lO2][lNO] = {};
	double kvvs_up_no_n2[lNO][lN2] = {}; 
	double kvvs_up_no_o2[lNO][lO2] = {};

        double mNO = 4.9826300488143997e-26; // NO
        double mN2 = 4.6517344343135997e-26; // N2
        double mO2 = 5.3135256633152e-26;    // O2
        double coll_mass = 0.;

    	coll_mass = (mN2 * mO2) / (mN2 + mO2);
	for (int i=0; i<lN2; ++i) {
	    for (int j=0; j<lO2; ++j) {

   	        kvvs_down_n2_o2[i][j] = k_VV_SSH(T, i, j, coll_mass, 
		    /*diameter=*/3.621e-10, /*omega_e=*/1., 
		            /*epsilon=*/1.3465465015559998e-21, 
			        /*r_e=*/1.20752e-10);

   	        kvvs_up_n2_o2[i][j] = kvvs_down_n2_o2[i][j] * 
	            exp((deps_n2[i] - deps_o2[j]) / KBT);
	    }
	}

    	coll_mass = (mN2 * mNO) / (mN2 + mNO);
	for (int i=0; i<lN2; ++i) {
	    for (int j=0; j<lNO; ++j) {

   	        kvvs_down_n2_no[i][j] = k_VV_SSH(T, i, j, coll_mass, 
		    /*diameter=*/3.621e-10, /*omega_e=*/1., 
		            /*epsilon=*/1.3465465015559998e-21, 
			        /*r_e=*/1.20752e-10);

    	        kvvs_up_n2_no[i][j] = kvvs_down_n2_no[i][j] * 
		    exp((deps_n2[i] - deps_no[j]) / KBT);
	    }
	}

    	coll_mass = (mN2 * mO2) / (mN2 + mO2);
	for (int i=0; i<lO2; ++i) {
	    for (int j=0; j<lN2; ++j) {

   	        kvvs_down_o2_n2[i][j] = k_VV_SSH(T, i, j, coll_mass, 
		    /*diameter=*/3.458e-10, /*omega_e=*/1., 
		            /*epsilon=*/1.48281651048e-21, 
			        /*r_e=*/1.20752e-10);

    	        kvvs_up_o2_n2[i][j] = kvvs_down_o2_n2[i][j] * 
		    exp((deps_o2[i] - deps_n2[j]) / KBT);
	    }
	}

	// O2+NO
    	coll_mass = (mO2 * mNO) / (mO2 + mNO);
	for (int i=0; i<lO2; ++i) {
	    for (int j=0; j<lNO; ++j) {

   	        kvvs_down_o2_no[i][j] = k_VV_SSH(T, i, j, coll_mass, 
		    /*diameter=*/3.458e-10, /*omega_e=*/1., 
		            /*epsilon=*/1.48281651048e-21, 
			        /*r_e=*/1.20752e-10);

    	        kvvs_up_o2_no[i][j] = kvvs_down_o2_no[i][j] * 
		    exp((deps_o2[i] - deps_no[j]) / KBT);
	    }
	}

    	coll_mass = (mNO * mN2) / (mNO + mN2);
	for (int i=0; i<lNO; ++i) {
	    for (int j=0; j<lN2; ++j) {

   	        kvvs_down_no_n2[i][j] = k_VV_SSH(T, i, j, coll_mass, 
		    /*diameter=*/3.47e-10, /*omega_e=*/1., 
		            /*epsilon=*/1.6429717387999998e-21, 
			        /*r_e=*/1.15077e-10);

    	        kvvs_up_no_n2[i][j] = kvvs_down_no_n2[i][j] * 
		    exp((deps_no[i] - deps_n2[j]) / KBT);
	    }
	}

    	coll_mass = (mNO * mO2) / (mNO + mO2);
	for (int i=0; i<lNO; ++i) {
	    for (int j=0; j<lO2; ++j) {

   	        kvvs_down_no_n2[i][j] = k_VV_SSH(T, i, j, coll_mass, 
		    /*diameter=*/3.47e-10, /*omega_e=*/1., 
		            /*epsilon=*/1.6429717387999998e-21, 
			        /*r_e=*/1.15077e-10);

    	        kvvs_up_no_o2[i][j] = kvvs_down_no_o2[i][j] * 
		    exp((deps_no[i] - deps_o2[j]) / KBT);
	    }
	}

        double R_w_VV[m_mixture.nMolecules()] = {}; // NO N2 O2
	double fac0, fac1, fac2, fac3, fac4, fac5 = 0.;
	int i, k = 0;
    
	// N2+NO
	for (i=1; i<lNO-1; ++i) { 
	    for (k=1; k<lN2-1; ++k) {

		fac0 += nN2[k]   * kvvs_down_no_n2[i][k];
	        fac1 += nN2[k+1] * kvvs_up_no_n2[i-1][k];
	 	fac2 += nN2[k+1] * kvvs_up_no_n2[i][k] +
	            	nN2[k]   * kvvs_down_no_n2[i-1][k];
	    }

	    for (k=1; k<lO2-1; ++k) {

		fac3 += nO2[k]   * kvvs_down_no_o2[i][k];
		fac4 += nO2[k+1] * kvvs_up_no_o2[i-1][k];
		fac5 += nO2[k+1] * kvvs_up_no_o2[i][k] +
                        nO2[k]   * kvvs_down_no_o2[i-1][k];
	    }

	    R_w_VV[0] += i * (nNO[i+1] * fac0 +
			      nNO[i-1] * fac1 -
			      nNO[i]   * fac2 + 
			      nNO[i+1] * fac3 + 
			      nNO[i-1] * fac4 -
			      nNO[i]   * fac5);
	}

	fac0=0; fac1=0.; fac2=0.; fac3=0.; fac4=0.; fac5=0.;
	i=0;
	    for (k=1; k<lN2-1; ++k) {

		fac0 += nN2[k]   * kvvs_down_no_n2[i][k];
	    	fac1 += nN2[k+1] * kvvs_up_no_n2[i][k];
	    }

	    for (k=1; k<lO2-1; ++k) {

		fac2 += nO2[k]   * kvvs_down_no_o2[i][k];
		fac3 += nO2[k+1] * kvvs_up_no_o2[i][k];
	    }

	    R_w_VV[0] += i * (nNO[i+1] * fac0 -
			      nNO[i]   * fac1 +
			      nNO[i+1] * fac2 - 
			      nNO[i]   * fac3);

	fac0=0; fac1=0.; fac2=0.; fac3=0.; fac4=0.; fac5=0.;
	i=lNO;
	    for (k=1; k<lN2-1; ++k) {

		fac0 += nN2[k+1] * kvvs_down_no_n2[i-1][k];
	        fac1 += nN2[k]   * kvvs_up_no_n2[i-1][k];
	    }

	    for (k=1; k<lO2-1; ++k) {	    

		fac2 += nO2[k+1] * kvvs_down_no_o2[i-1][k];
		fac3 += nO2[k]   * kvvs_up_no_o2[i-1][k];
	    }

	    R_w_VV[0] += i * (nNO[i-1] * fac0 -
			      nNO[i]   * fac1 +
			      nNO[i-1] * fac2 -
			      nNO[i]   * fac3);

	// N2+O2
	fac0=0; fac1=0.; fac2=0.; fac3=0.; fac4=0.; fac5=0.;
	for (i=1; i<lN2-1; ++i) { 

	    for (k=1; k<lO2-1; ++k) {

		fac0 += nN2[k]   * kvvs_down_n2_o2[i][k];
	        fac1 += nN2[k+1] * kvvs_up_n2_o2[i-1][k];
	 	fac2 += nN2[k+1] * kvvs_up_n2_o2[i][k] +
	            	nN2[k]   * kvvs_down_n2_o2[i-1][k];
	    }

	    for (k=1; k<lNO-1; ++k) {

		fac3 += nNO[k]   * kvvs_down_n2_no[i][k];
		fac4 += nNO[k+1] * kvvs_up_n2_no[i-1][k];
		fac5 += nNO[k+1] * kvvs_up_n2_no[i][k] +
                        nNO[k]   * kvvs_down_n2_no[i-1][k];
	    }

	    R_w_VV[1] += i * (nNO[i+1] * fac0 +
			      nNO[i-1] * fac1 -
			      nNO[i]   * fac2 + 
			      nNO[i+1] * fac3 + 
			      nNO[i-1] * fac4 -
			      nNO[i]   * fac5);
	}

	fac0=0; fac1=0.; fac2=0.; fac3=0.; fac4=0.; fac5=0.;
	i=0;

	    for (k=1; k<lO2-1; ++k) {

		fac0 += nO2[k]   * kvvs_down_n2_o2[i][k];
	        fac1 += nO2[k+1] * kvvs_up_n2_o2[i][k];
	    }        	

	    for (k=1; k<lNO-1; ++k) {

		fac3 += nNO[k]   * kvvs_down_n2_no[i][k];
		fac4 += nNO[k+1] * kvvs_up_n2_no[i][k];
	    }

	    R_w_VV[1] += i * (nN2[i+1] * fac0 -
			      nN2[i]   * fac1 +
			      nN2[i+1] * fac2 -
			      nN2[i]   * fac3);

	fac0=0; fac1=0.; fac2=0.; fac3=0.; fac4=0.; fac5=0.;
	i=lN2;

	    for (k=1; k<lO2-1; ++k) {

		fac0 += nO2[k+1] * kvvs_down_n2_o2[i-1][k];
	        fac1 += nO2[k]   * kvvs_up_n2_o2[i-1][k];
	    }        	

	    for (k=1; k<lNO-1; ++k) {

		fac3 += nNO[k+1] * kvvs_down_n2_no[i-1][k];
		fac4 += nNO[k] * kvvs_up_n2_no[i-1][k];
	    }

	    R_w_VV[1] += i * (nN2[i+1] * fac0 -
			      nN2[i]   * fac1 +
			      nN2[i+1] * fac2 -
			      nN2[i]   * fac3);

	// O2+NO
	fac0=0; fac1=0.; fac2=0.; fac3=0.; fac4=0.; fac5=0.;
	for (i=1; i<lO2-1; ++i) { 

	    for (k=1; k<lN2-1; ++k) {

		fac0 += nN2[k]   * kvvs_down_o2_n2[i][k];
	    	fac1 += nN2[k+1] * kvvs_up_o2_n2[i-1][k];
		fac2 += nN2[k+1] * kvvs_up_o2_n2[i][k] +
			nN2[k]   * kvvs_down_o2_n2[i-1][k];
	    }

	    for (k=1; k<lNO-1; ++k) {
		    
		fac3 += nNO[k]   * kvvs_down_o2_no[i][k];
		fac4 += nNO[k+1] * kvvs_up_o2_no[i-1][k];
		fac5 += nNO[k+1] * kvvs_up_o2_no[i][k] +
			nNO[k]   * kvvs_down_o2_no[i-1][k];
	    }

	    R_w_VV[2] += i * (nO2[i+1] * fac0 +
			      nO2[i-1] * fac1 -
			      nO2[i]   * fac2 +
			      nO2[i+1] * fac3 +
			      nO2[i-1] * fac4 -
			      nO2[i]   * fac5);
	}

	fac0=0; fac1=0.; fac2=0.; fac3=0.; fac4=0.; fac5=0.;
	i=0;
	    for (k=1; k<lN2-1; ++k) {

		fac0 += nN2[k]   * kvvs_down_o2_n2[i][k];
		fac1 += nN2[k+1] * kvvs_up_o2_n2[i][k];
	    }

	    for (k=1; k<lNO-1; ++k) {

		fac2 += nNO[k]   * kvvs_down_o2_no[i][k];
		fac3 += nNO[k+1] * kvvs_up_o2_no[i][k];

	    }

	    R_w_VV[2] += i * (nO2[i+1] * fac0 -
			      nO2[i]   * fac1 +
			      nO2[i+1] * fac2 -
			      nO2[i]   * fac3);

	fac0=0; fac1=0.; fac2=0.; fac3=0.; fac4=0.; fac5=0.;
	i=lO2;
	    for (k=1; k<lN2-1; ++k) {

		fac0 += nN2[k+1] * kvvs_down_o2_n2[i-1][k];
		fac1 += nN2[k]   * kvvs_up_o2_n2[i-1][k];
	    }

	    for (k=1; k<lNO-1; ++k) {

		fac2 += nNO[k+1] * kvvs_down_o2_no[i-1][k];
	    	fac3 += nNO[k]   * kvvs_up_o2_no[i-1][k];
	    }

	    R_w_VV[2] += i * (nO2[i-1] * fac0 -
			      nO2[i]   * fac1 +
		  	      nO2[i-1] * fac2 -
			      nO2[i]   * fac3);	      

        // Compute species concentrations (mol/m^3)
        // TODO: check dimensions!
        for (int i=0; i<m_mixture.nMolecules(); ++i) {
            mp_rVV[i] = R_w_VV[i];
    	std::cout <<  " mp_rVV["<<i<<"]" << mp_rVV[i] << std::endl;
        }
    }

private:

    int m_ns;
    int m_transfer_offset;

    double* mp_wrk1;
    double* mp_wrk2;
    double* mp_Mw;
    //double* mp_hv;
    //double* mp_hveq;
    double* mp_rVV;

    // SSH model (harmonic only)
    inline double k_VV_SSH(double T, int i, int k, double coll_mass, 
	double diameter, double omega_e, double epsilon, double r_e);

    inline double P_SSH_VV_01(double T, double omega_e, double epsilon, 
	double osc_mass, double diameter, double r_e);

    inline double p_Z_coll(double T, double n, double coll_mass, double diameter);
};

//==============================================================================

   inline double R_vibr_VV::k_VV_SSH(double T, int i, int k, double coll_mass, 
       double diameter, double omega_e, double epsilon, double r_e) 
   {
       // eq. 81
       return p_Z_coll(T, 1., coll_mass, diameter) * (i+1)*(k+1) * 
           P_SSH_VV_01(T, coll_mass, omega_e, epsilon, diameter, r_e);
   }

//==============================================================================

    // Collision frequency
    inline double R_vibr_VV::p_Z_coll(double T, double n, double coll_mass, 
	double diameter) 
    {
        // eq. 82(a) assuming constant diameter
        return 4. * PI * n * diameter * diameter * 
	    sqrt(KB * T / (2. * PI * coll_mass));
    }

//==============================================================================

    inline double R_vibr_VV::P_SSH_VV_01(double T, double omega_e, 
	double epsilon, double osc_mass, double diameter, double r_e) 
    {
        // eq. 86(a)
        double alpha = 17.5 / diameter; 
    	double lambda = 0.5; // for diatomic homonuclear molecules
    	double omega = 2. * PI * C0 * omega_e; //ang. freq. of the oscillator
    
	// eq. 84
	return /*pow(lambda, 4)=*/0.0625 * 4. * KB * T / osc_mass * 
	    alpha * alpha / omega / omega;
    }

//==============================================================================

// Register the transfer model
Utilities::Config::ObjectProvider<R_vibr_VV, TransferModel> 
    r_vibr_vv("R_vibr_VV");

    } // namespace Transfer
} // namespace Mutation 
