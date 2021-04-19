/**
 * @file R_w_23.cpp
 *
 * @brief Implementation of R_w_23.
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

#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;
using namespace Mutation;

namespace Mutation {
    namespace Transfer {

/**
 * Represents the production term due to dissociation/recombination
 * which appears in the vibrational energy equations.
 */
class R_w_react_23 : public TransferModel
{
public:

    R_w_react_23(Mixture& mix)
        : TransferModel(mix)
    {
        m_ns              = m_mixture.nSpecies();
        m_transfer_offset = m_mixture.hasElectrons() ? 1 : 0;

        mp_Mw = new double [m_ns];
        for(int i = 0; i < m_ns; ++i)
            mp_Mw[i] = m_mixture.speciesMw(i);

	mp_r_w_23 = new double [m_mixture.nEnergyEqns()-1];
    }

    virtual ~R_w_react_23()
    {
        delete [] mp_Mw;
	delete [] mp_r_w_23;
    }

    double source() {};
    void rVT(double* const mp_rVT) {};
    void rVV(double* const p_rVv) {};
    void R_w_22(double* const p_r_w_22) {};

    // Compute the production term R_23 due to dissociation/recombination
    void R_w_23(double* const mp_r_w_23)
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
    
        // Number of vibrational levels (harmonic oscillator)
        int lN2 = 33; 
        int lO2 = 26;
        int lNO = 27;
    
        // Vibrational energy arrays [J]
        static double ve[87];
        static double veN2[33];
        static double veO2[26];
        static double veNO[28];
    
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
    
        // Dissociation energy [J]
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
    	//  		                      i*veNO[0]/(KB*Tv[0]));
        }

        double z_vibr_T_Tv_N2 = 0.;
        for (int i=0; i<lN2; ++i) {
            z_vibr_T_Tv_N2 += exp(-veN2[i]/(KB*Tv[1]));
            //z_vibr_T_Tv_N2 += exp(-(veN2[i]-i*veN2[0])/(KB*T)-
    	//       	                      i*veN2[0]/(KB*Tv[1]));
        }

        double z_vibr_T_Tv_O2 = 0.;
        for (int i=0; i<lO2; ++i) {
            z_vibr_T_Tv_O2 += exp(-veO2[i]/(KB*Tv[2]));
            //z_vibr_T_Tv_O2 += exp(-(veO2[i]-i*veO2[0])/(KB*T)-
    	//		                      i*veO2[0]/(KB*Tv[2]));
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
    
        // Masses 
        double m[5]; // [kg]
        m[0] = 2.3258672171567998e-26; // N
        m[1] = 2.6567628316576e-26;    // O
        m[2] = 4.9826300488143997e-26; // NO
        m[3] = 4.6517344343135997e-26; // N2
        m[4] = 5.3135256633152e-26;    // O2
    
        // Diameters
        double d[5]; // [m]
        d[0] = 3.298e-10;
        d[1] = 2.75e-10;
        d[2] = 3.4061e-10;
        d[3] = 3.40385035355259e-10;
        d[4] = 3.5155e-10;
        
        // Compute state-to-state diss/rec coefficients, k_c,diss^d(0)
        // according to the Rigid-Sphere (RS) model.
        std::cout << " Dissociation STS rates ... " << std::endl;
        double kcd_no[5][lNO]; 
        double coll_mass = 0.;
        double diameter = 0.;
        double interNO [5][3] = {}; // 1st: interaction, 2nd: property
        interNO[0][0] = 1.82659281313541e-13; 	// A for N + NO 
        interNO[0][1] = 0.;			// n
        interNO[0][2] = 1.0423910520000002e-18;	// Ea
        interNO[1][0] = 1.82659281313541e-13; 	// A for O + NO 
        interNO[1][1] = 0.;			// n
        interNO[1][2] = 1.0423910520000002e-18;	// Ea
        interNO[2][0] = 1.82659281313541e-13; 	// A for NO + NO 
        interNO[2][1] = 0.;			// n
        interNO[2][2] = 1.0423910520000002e-18;	// Ea
        interNO[3][0] = 1.16237724472253e-08; 	// A for N2 + NO 
        interNO[3][1] = -1.6;			// n
        interNO[3][2] = 1.5628962528000002e-18;	// Ea
        interNO[4][0] = 2e-10;		 	// A for O2 + NO 
        interNO[4][1] = -1.;			// n
        interNO[4][2] = 1.043e-18;		// Ea
        double A, nn, Ea = 0.;
        for (int i=0; i<m_mixture.nSpecies(); ++i) { // N O NO N2 O2
    	    coll_mass = (m[i] * m[2]) / (m[i] + m[2]);
    	    diameter = 0.5 * (d[i] + d[2]);
            for (int j=0; j<lNO; ++j) {
    
                // Treanor-Marrone, 3T model
		// *************************
    	    	A = interNO[i][0];
    	    	nn = interNO[i][1];
    	    	Ea = interNO[i][2];
    	    	kcd_no[i][j] = k_Arrhenius(T, A, nn, Ea) 
    	    	    * Z_diss(T, 3.*T, veNO, j);
    
                // Treanor-Marrone, D/6K model
		// ***************************
    	    	//kcd_no[i][j] = k_Arrhenius(T, A, nn, Ea) 
    	    	//    * Z_diss(T, de[0]/(6.*KB), veNO, j);
    	    
     	        // Rigid-sphere (RS) model	    
		// ***********************
                //kcd_no[i][j] = k_diss_RS(T, coll_mass, diameter, de[0], veNO[j],
    	        //    /*center_of_mass=*/true);
    
     	        std::cout << i << " " << j << " " << T << " " << coll_mass 
    		   << " " << diameter << " " << de[2] << " " << veNO[j] 
    		   << " " << kcd_no[i][j] << std::endl;
    	    }
        }
    
        double interN2 [5][3] = {}; 
        interN2[0][0] = 2.657e-08; 		// A for N + N2
        interN2[0][1] = -1.6;			 
        interN2[0][2] = 1.5628962528000002e-18;
        interN2[1][0] = 4.98161676309657e-08; 	// A for O + N2 
        interN2[1][1] = -1.6;			 
        interN2[1][2] = 1.5628962528000002e-18;	
        interN2[2][0] = 1.16237724472253e-08; 	// A for NO + N2
        interN2[2][1] = -1.6;			 
        interN2[2][2] = 1.5628962528000002e-18;
        interN2[3][0] = 6.144e-09; 		// A for N2 + N2
        interN2[3][1] = -1.6;		
        interN2[3][2] = 1.5628962528000002e-18;
        interN2[4][0] = 1.16237724472253e-08; 	// A for O2 + N2
        interN2[4][1] = -1.6;		
        interN2[4][2] = 1.5628962528000002e-18;	
        double kcd_n2[5][lN2]; 
        for (int i=0; i<m_mixture.nSpecies(); ++i) { // N O NO N2 O2
    	    coll_mass = (m[i] * m[3]) / (m[i] + m[3]);
    	    diameter = 0.5 * (d[i] + d[3]);
            for (int j=0; j<lN2; ++j) {
    
    	        A = interN2[i][0];
    	        nn = interN2[i][1];
    	        Ea = interN2[i][2];
    	        kcd_n2[i][j] = k_Arrhenius(T, A, nn, Ea) 
		    * Z_diss(T, 3*T, veN2, j);
    
    	        //kcd_n2[i][j] = k_Arrhenius(T, A, nn, Ea) 
    	        //    * Z_diss(T, de[1]/(6.*KB), veN2, j);
    	        
                //kcd_n2[i][j] = k_diss_RS(T, coll_mass, diameter, de[1], veN2[j],
    	        //    /*center_of_mass=*/true);
    
     	        std::cout << i << " " << j << " " << T << " " << coll_mass 
    		   << " " << diameter << " " << de[1] << " " << veN2[j] 
    		   << " " << kcd_n2[i][j] << std::endl;
    	    }
        }
    
        double interO2 [5][3] = {}; 
        interO2[0][0] = 3.32107784206438e-09;	// A for N + O2
        interO2[0][1] = -1.5;			 
        interO2[0][2] = 8.24938614e-19;
        interO2[1][0] = 1.66053892103219e-08; 	// A for O + O2 
        interO2[1][1] = -1.5;			 
        interO2[1][2] = 8.24938614e-19;	
        interO2[2][0] = 3.32107784206438e-09; 	// A for NO + O2
        interO2[2][1] = -1.5;			 
        interO2[2][2] = 8.24938614e-19;
        interO2[3][0] = 1.16237724472253e-08; 	// A for N2 + O2
        interO2[3][1] = -1.6;		
        interO2[3][2] = 1.5628962528000002e-18;
        interO2[4][0] = 3.32107784206438e-09; 	// A for O2 + O2
        interO2[4][1] = -1.5;		
        interO2[4][2] = 8.24938614e-19;	
        double kcd_o2[5][lO2]; 
        for (int i=0; i<m_mixture.nSpecies(); ++i) { // N O NO N2 O2
    	    coll_mass = (m[i] * m[4]) / (m[i] + m[4]);
    	    diameter = 0.5 * (d[i] + d[4]);
            for (int j=0; j<lO2; ++j) {
    
    	        A = interO2[i][0];
    	        nn = interO2[i][1];
    	        Ea = interO2[i][2];
    	        kcd_o2[i][j] = k_Arrhenius(T, A, nn, Ea) 
                        * Z_diss(T, 3.*T, veO2, j);
    	
    	        //kcd_o2[i][j] = k_Arrhenius(T, A, nn, Ea) 
    	        //    * Z_diss(T, de[2]/(6.*KB), veO2, j);
    
                //kcd_o2[i][j] = k_diss_RS(T, coll_mass, diameter, de[2], veO2[j],
    	        //    /*center_of_mass=*/true);
    
     	        std::cout << i << " " << j << " " << T << " " << coll_mass 
    		   << " " << diameter << " " << de[2] << " " << veO2[j] 
    		   << " " << kcd_o2[i][j] << std::endl;
    	    }
        }
    
        // Compute the rotational partition function
        // according to a simplified formula ... TODO: better formula
        double Be[3] = {199.8, 143.77, 167.20}; // m^-1
        int sigma[3] = {2, 2, 1};
        double theta_r[3] = {};
        double z_rot[3] = {};
        for (int i=0; i<3; ++i) {
            theta_r[i] = Be[i]*HP*C0/KB;
            z_rot[i] = T/(sigma[i]*theta_r[i]);
    	    std::cout << "theta & zrot:" 
    	              << theta_r[i] << " " << z_rot[i] << std::endl;
        }
    
        // Compute STS recombination coefficients, k_rec,c^d(0)
        std::cout << " Recombination STS rates ... " << std::endl;
        double fac0 = HP * HP * HP * pow(2. * PI * KB * T, -1.5);
        double fac1 = (m[2]/(m[0]*m[1])) * fac0 * z_rot[0];
    
        std::cout << " Recombination STS rates NO ... " << std::endl;
        double kcr_no[5][lNO]; 
        for (int i=0; i<m_mixture.nSpecies(); ++i) {
            for (int j=0; j<lNO; ++j) {
                kcr_no[i][j] = kcd_no[i][j] * exp(-(veNO[j]-de[0])/KBT) * fac1;
    	        std::cout << i << " " << j << " " << kcr_no[i][j] << std::endl;
    	    }
        }
    
        std::cout << " Recombination STS rates N2 ... " << std::endl;
        fac1 = (m[3]/(m[0]*m[0])) * fac0 * z_rot[1];
        double kcr_n2[5][lN2]; 
        for (int i=0; i<m_mixture.nSpecies(); ++i) {
            for (int j=0; j<lN2; ++j) {
                kcr_n2[i][j] = kcd_n2[i][j] * exp(-(veN2[j]-de[1])/KBT) * fac1;
    	        std::cout << i << " " << j << " " << kcr_n2[i][j] << std::endl;
    	    }
        }
    
        std::cout << " Recombination STS rates O2 ... " << std::endl;
        fac1 = (m[4]/(m[1]*m[1])) * fac0 * z_rot[2];
        double kcr_o2[5][lO2]; 
        for (int i=0; i<m_mixture.nSpecies(); ++i) {
            for (int j=0; j<lO2; ++j) {
                kcr_o2[i][j] = kcd_o2[i][j] * exp(-(veO2[j]-de[2])/KBT) * fac1;
    	        std::cout << i << " " << j << " " << kcr_o2[i][j] << std::endl;
    	    }
        }
    
        std::cout << " Recombination STS rates FINISH ... " << std::endl;
	// In this case, we don't need to compute such terms
	// and we can directly sum-up all indices and calculate
	// the production term, R_w_VV2.
	// eq. pp 78 of Kustova & Nagnibeda 2009.
        // Compute MT diss/rec coefficients, k_c^d(0), eq. 3.94.
        //double kd_no[m_mixture.nSpecies()];
        //double kr_no[m_mixture.nSpecies()];
        //double kd_n2[m_mixture.nSpecies()];
        //double kr_n2[m_mixture.nSpecies()];
        //double kd_o2[m_mixture.nSpecies()];
        //double kr_o2[m_mixture.nSpecies()];
    
        //// Compute the summation term in eq. 3.94.
        //// NO
        //for (int i=0; i<m_mixture.nSpecies(); ++i) {
    	//    fac0 = 0., fac1 = 0.; 
        //    for (int j=0; j<lNO; ++j) {
    
    	//        // Using the Boltzmann distribution
    	//        fac0 += nNO[j] * kcd_no[i][j];
    
    	//        // Using the generalized Treanor-Marrone distribution	
        //        //fac0 += exp(-((veNO[j]-j*veNO[0])/KBT) -
        //        //                       j*veNO[0]/(KB*Tv[0])) * kcd_no[i][j];
    	//        fac1 += kcr_no[i][j];
        //    }
    	//    kd_no[i] = fac0 / nd[2];
        //    //kd_no[i] = fac0 / z_vibr_T_Tv_NO;
        //    kr_no[i] = fac1;
    	//    std::cout << i << " kdNO: " << kd_no[i] 
    	//    	           << " krNO: " << kr_no[i] << std::endl;
        //}
    
        //// N2
        //for (int i=0; i<m_mixture.nSpecies(); ++i) {
    	//    fac0 = 0., fac1 = 0.; 
        //    for (int j=0; j<lN2; ++j) {
    
    	//        fac0 += nN2[j] * kcd_n2[i][j];
    
        //        //fac0 += exp(-((veN2[j]-j*veN2[0])/KBT) -
        //        //                       j*veN2[0]/(KB*Tv[1])) * kcd_n2[i][j];
    	//        fac1 += kcr_n2[i][j];
        //    }
    	//    kd_n2[i] = fac0 / nd[3];
        //    //kd_n2[i] = fac0 / z_vibr_T_Tv_N2;
        //    kr_n2[i] = fac1;
    	//    std::cout << i << " kdN2: " << kd_n2[i] 
    	// 	           << " krN2: " << kr_n2[i] << std::endl;
        //}
    
        //// O2
        //for (int i=0; i<m_mixture.nSpecies(); ++i) {
    	//    fac0 = 0., fac1 = 0.; 
        //    for (int j=0; j<lO2; ++j) {
    
    	//        fac0 += nO2[j] * kcd_o2[i][j];
    
        //        //fac0 += exp(-((veO2[j]-j*veO2[0])/KBT) -
        //        //                       j*veO2[0]/(KB*Tv[2])) * kcd_o2[i][j];
    	//        fac1 += kcr_o2[i][j];
        //    }
    	//    kd_o2[i] = fac0 / nd[4];
        //    //kd_o2[i] = fac0 / z_vibr_T_Tv_O2;
        //    kr_o2[i] = fac1;
    	//    std::cout << i << " kdO2: " << kd_o2[i] 
    	//	           << " krO2: " << kr_o2[i] << std::endl;
        //}
    
	fac0 = 0.;
	for (int i=0; i<lNO; ++i) {
	    fac0 += i * ( nd[0] * (nd[0]*nd[1]*kcr_no[0][i]-nd[2]*kcd_no[0][i]) + 
			  nd[1] * (nd[0]*nd[1]*kcr_no[1][i]-nd[2]*kcd_no[1][i]) +
			  nd[2] * (nd[0]*nd[1]*kcr_no[2][i]-nd[2]*kcd_no[2][i]) +
			  nd[3] * (nd[0]*nd[1]*kcr_no[3][i]-nd[2]*kcd_no[3][i]) +
			  nd[4] * (nd[0]*nd[1]*kcr_no[4][i]-nd[2]*kcd_no[4][i]) );
	}
	mp_r_w_23[0] = fac0;
        //std::cout << " mp_r_w_23[0] = " << mp_r_w_23[0] << std::endl;

	fac0 = 0.;
	for (int i=0; i<lN2; ++i) {
	    fac0 += i * ( nd[0] * (nd[0]*nd[0]*kcr_n2[0][i]-nd[3]*kcd_n2[0][i]) + 
			  nd[1] * (nd[0]*nd[0]*kcr_n2[1][i]-nd[3]*kcd_n2[1][i]) +
			  nd[2] * (nd[0]*nd[0]*kcr_n2[2][i]-nd[3]*kcd_n2[2][i]) +
			  nd[3] * (nd[0]*nd[0]*kcr_n2[3][i]-nd[3]*kcd_n2[3][i]) +
			  nd[4] * (nd[0]*nd[0]*kcr_n2[4][i]-nd[3]*kcd_n2[4][i]) );
	}
        mp_r_w_23[1] = fac0;
        //std::cout << " mp_r_w_23[1] = " << mp_r_w_23[1] << std::endl;

	fac0 = 0.;
	for (int i=0; i<lO2; ++i) {
	    fac0 += i * ( nd[0] * (nd[1]*nd[1]*kcr_o2[0][i]-nd[4]*kcd_o2[0][i]) + 
			  nd[1] * (nd[1]*nd[1]*kcr_o2[1][i]-nd[4]*kcd_o2[1][i]) +
			  nd[2] * (nd[1]*nd[1]*kcr_o2[2][i]-nd[4]*kcd_o2[2][i]) +
			  nd[3] * (nd[1]*nd[1]*kcr_o2[3][i]-nd[4]*kcd_o2[3][i]) +
			  nd[4] * (nd[1]*nd[1]*kcr_o2[4][i]-nd[4]*kcd_o2[4][i]) );
	}
        mp_r_w_23[2] = fac0;
        //std::cout << " mp_r_w_23[2] = " << mp_r_w_23[2] << std::endl;

        //for (int i=0; i<m_mixture.nEnergyEqns()-1; ++i) 
    	//    std::cout << " mp_r_w_23["<<i<<"] = " << mp_r_w_23[i] << std::endl;
    }

private:

    int m_ns;
    int m_transfer_offset;

    double* mp_wrk1;
    double* mp_wrk2;
    double* mp_Mw;
    double* mp_r_w_23;

    inline double k_diss_RS(double T, double coll_mass, double diameter,
        double diss_energy, double electron_vibr_energy,
            bool center_of_mass=true);

    inline double k_diss_VSS(double T, double coll_mass, double vss_c_cs,
        double vss_omega, double diss_energy, double electron_vibr_energy,
            bool center_of_mass=true);

    inline double integral_diss_RS(double T, int degree, double coll_mass,
        double diameter, double diss_energy, double vibr_energy,
            bool center_of_mass=true);

    inline double integral_diss_VSS(double T, int degree, double coll_mass,
            double vss_c_cs, double vss_omega, double diss_energy,
                double vibr_energy, bool center_of_mass=true);

    // Calculate the elastic cross-section sigma_{tot} using the Rigid Sphere
    // (RS) interaction potential
    inline double crosssection_elastic_RS(double diameter);

    inline double crosssection_elastic_VSS(double rel_vel, double vss_c_cs,
        double vss_omega);

    inline double crosssection_elastic_VSS(double rel_vel, double coll_mass,
        double vss_c, double vss_omega);
    
    // Calculate the dissociation cross-section (VSS model),
    // accounting for vibrational energy
    inline double crosssection_diss_VSS(double rel_vel, double coll_mass,
        double vss_c_cs, double vss_omega, double diss_energy,
            double vibr_energy, bool center_of_mass=true);

    // Calculate the minimum velocity for which a dissociation probability is
    // non-zero, accounting for vibrational energy
    inline double min_vel_diss(double coll_mass, double diss_energy,
        double vibr_energy);

    // Calculate the dissociation cross-section (RS model),
    // accounting for vibrational energy
    inline double crosssection_diss_RS(double rel_vel, double coll_mass,
        double diameter, double diss_energy, double vibr_energy,
            bool center_of_mass=true);
    
    // Calculate the probability of dissociation, accounting for vibr. energy
    inline double p_probability_diss(double rel_vel, double coll_mass,
        double diss_energy, double vibr_energy, bool center_of_mass=true);

    inline double k_Arrhenius(double T, double A, double n, double Ea);

    // SSH model
    inline double k_VT_SSH(double T, int i, double coll_mass, double diameter,
        double omega_e, double epsilon, double r_e); // harmonic

    inline double k_VT_SSH(double T, int i, double coll_mass, double diameter,
        double omega_e, double epsilon, double r_e, double Delta_E_vibr,
            double vibr_energy_1); // anharmonic

    inline double p_Z_coll(double T, double n, double coll_mass, double diameter);

    inline double P_SSH_VT_10(double T, double coll_mass, double omega_e,
        double epsilon, double diameter, double r_e);

    inline double Z_diss(double T, double U, const double *ve, int i);
    inline double p_Z_vibr_eq(double T, const double *ve);

};

//==============================================================================

    inline double R_w_react_23::k_diss_RS(double T, double coll_mass, 
	double diameter, double diss_energy, double vibr_energy, 
	    bool center_of_mass) 
    {
        return 8. * integral_diss_RS(T, 0, coll_mass, diameter, diss_energy, 
	    vibr_energy, center_of_mass);
    }

//==============================================================================

    inline double R_w_react_23::k_diss_VSS(double T, double coll_mass, 
	double vss_c_cs, double vss_omega, double diss_energy, 
	    double vibr_energy, bool center_of_mass) 
    {
        return 8. * integral_diss_VSS(T, 0, coll_mass, vss_c_cs, vss_omega, 
	    diss_energy, vibr_energy, center_of_mass);
    }

//==============================================================================

    inline double R_w_react_23::integral_diss_RS(double T, int degree,
        double coll_mass, double diameter, double diss_energy, 
	    double vibr_energy, bool center_of_mass) 
    {
        double conversion = sqrt(2. * KB * T / coll_mass);

        auto integrand = [=](double g)
        {
	    return pow(g, 2 * degree + 3) * 
	        crosssection_diss_RS(conversion * g, coll_mass, diameter, 
		    diss_energy, vibr_energy, center_of_mass) * exp(-g * g); 
	};

        return sqrt(KB * T / (2. * PI * coll_mass)) * 
	    Numerics::integrate_semi_inf(integrand, min_vel_diss(coll_mass, 
                diss_energy, vibr_energy) / conversion);
    }

//==============================================================================

    inline double R_w_react_23::integral_diss_VSS(double T, int degree, 
	double coll_mass, double vss_c_cs, double vss_omega, double diss_energy, 
	    double vibr_energy, bool center_of_mass) 
    {

        double conversion = sqrt(2. * KB * T / coll_mass);

        auto integrand = [=](double g)
	{
            return pow(g, 2 * degree + 3) * crosssection_diss_VSS(
	        conversion * g, coll_mass, vss_c_cs, vss_omega, 
	            diss_energy, vibr_energy, center_of_mass) * 
		        exp(-g * g); 
	};

        return sqrt(KB * T / (2. * PI * coll_mass)) * 
	    Numerics::integrate_semi_inf(integrand, min_vel_diss(coll_mass, 
	        diss_energy, vibr_energy) / conversion);
    }

//==============================================================================

    inline double R_w_react_23::crosssection_diss_RS(double rel_vel, 
	double coll_mass, double diameter, double diss_energy, 
	    double vibr_energy, bool center_of_mass) 
    {
        return crosssection_elastic_RS(diameter) * 
	    p_probability_diss(rel_vel, coll_mass, diss_energy, vibr_energy, 
	        center_of_mass);
    }

//==============================================================================

    inline double R_w_react_23::crosssection_diss_VSS(double rel_vel, 
	double coll_mass, double vss_c_cs, double vss_omega, double diss_energy, 
	    double vibr_energy, bool center_of_mass) 
    {   
        return crosssection_elastic_VSS(rel_vel, vss_c_cs, vss_omega) *
	    p_probability_diss(rel_vel, coll_mass, diss_energy, vibr_energy, 
	        center_of_mass);
    }

//==============================================================================

    inline double R_w_react_23::p_probability_diss(double rel_vel, 
	double coll_mass, double diss_energy, double vibr_energy, 
	    bool center_of_mass) 
    {

        double energy = vibr_energy + rel_vel * rel_vel * coll_mass / 2.;

        if (energy < diss_energy) {
            return 0.0;
        } else if (center_of_mass) {
            return 1.0 - diss_energy / energy;
        } else {
            return 1.0;
        }
    }

//==============================================================================

    inline double R_w_react_23::min_vel_diss(double coll_mass, 
	double diss_energy, double vibr_energy) 
    {
        return sqrt(2. * (diss_energy - vibr_energy) / coll_mass);
    }

//==============================================================================

    inline double R_w_react_23::crosssection_elastic_RS(double diameter) 
    {
        return PI * diameter * diameter;
    }

//==============================================================================

    inline double R_w_react_23::crosssection_elastic_VSS(double rel_vel, 
        double vss_c_cs, double vss_omega) 
    {
        return vss_c_cs * pow(rel_vel, 1. - 2. * vss_omega);
    }

//==============================================================================

    inline double R_w_react_23::crosssection_elastic_VSS(double rel_vel, 
	double coll_mass, double vss_c, double vss_omega) 
    {
        return vss_c * pow(coll_mass * rel_vel * rel_vel / (2. * KB), 
	    -vss_omega);
    }

//==============================================================================

    inline double R_w_react_23::k_Arrhenius(double T, double A, double n, 
        double Ea) 
    {
        return A * pow(T, n) * exp(-Ea / (KB * T));
    }

//==============================================================================

    // Compute the non-equilibrium factor using the Treanor-Marrone model
    inline double R_w_react_23::Z_diss(double T,double U,const double *ve,int i) 
    {
        return p_Z_vibr_eq(T, ve) * exp(ve[i] * (1./T + 1./U) / KB) 
            / p_Z_vibr_eq(-U, ve);
    }

//==============================================================================

    // Computes the equilibrium vibrational partition function
    inline double R_w_react_23::p_Z_vibr_eq(double T, const double *ve) 
    {
        double p_Z_vibr_eq = 0.;
        int size = sizeof(ve)/sizeof(ve[0]);
        for (int i=0; i<size; ++i) 
            p_Z_vibr_eq += exp(-ve[i] / (KB * T));
          
        return p_Z_vibr_eq;
    }

//==============================================================================

// Register the transfer model
Utilities::Config::ObjectProvider<R_w_react_23, TransferModel> 
    r_w_react_23("R_w_react_23");

    } // namespace Transfer
} // namespace Mutation 
