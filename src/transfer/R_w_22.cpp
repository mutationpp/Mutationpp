/**
 * @file R_w_22.cpp
 *
 * @brief Implementation of R_w_22.
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
 * Represents the production term due to chemical reactions
 * appearing in the vibrational energy equations.
 */
class R_w_react_22 : public TransferModel
{
public:

    R_w_react_22(Mixture& mix)
        : TransferModel(mix)
    {
        m_ns              = m_mixture.nSpecies();
        m_transfer_offset = m_mixture.hasElectrons() ? 1 : 0;

        mp_Mw = new double [m_ns];
        for(int i = 0; i < m_ns; ++i)
            mp_Mw[i] = m_mixture.speciesMw(i);

	mp_r_w_22 = new double [m_mixture.nEnergyEqns()-1];
    }

    virtual ~R_w_react_22()
    {
        delete [] mp_Mw;
	delete [] mp_r_w_22;
    }

    double source() {};
    void rVT(double* mp_rVT) {};
    void rVV(double* const p_rVv) {};
    void R_w_23(double* mp_r_w_23) {};

    // Compute the production term R_22 due to exchange reactions
    void R_w_22(double* mp_r_w_22)
    {
        double T = m_mixture.T();
        double KBT = KB * T;
        double Tv[3];
        for (int i=0; i<3; ++i)
            Tv[i] = m_mixture.Tvs(i); // NO, N2, O2

        const double * p_X = m_mixture.X();
        const double   n   = m_mixture.numberDensity();

        // Partial number densities
        double nd[m_mixture.nSpecies()] = {}; // [m^-3]
        std::cout << " Number densities ... " << std::endl;
        for (int i=0; i<m_mixture.nSpecies(); ++i) {
            nd[i] = p_X[i] * n; // N O NO N2 O2
            std::cout << i << " " << nd[i] << std::endl;
        }

        // Number of vibrational levels (harmonic oscillator)
        int lNO = 28;
        int lN2 = 33;
        int lO2 = 26;

        // Vibrational energy arrays [J]
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

        // Dissociation energy (ground level) [J]
        const double de[m_mixture.nMolecules()] = {1.0409017238088574e-18,
            1.5636156480913654e-18, 8.196091362099268e-19};

        // Non-equilibrium partition functions [-]
        double z_vibr_T_Tv_NO = 0.;
        for (int i=0; i<lNO; ++i) {

            // anharmonic (treanor-marrone)
            //z_vibr_T_Tv_NO += exp(-(veNO[i]-i*veNO[0])/KBT-
            //		                i*veNO[0]/(KB*Tv[0]));
	    		                
            // harmonic (boltzmann)
            z_vibr_T_Tv_NO += exp(-veNO[i]/(KB*Tv[0]));
        }

        double z_vibr_T_Tv_N2 = 0.;
        for (int i=0; i<lN2; ++i) {
            //z_vibr_T_Tv_N2 += exp(-(veN2[i]-i*veN2[0])/KBT-
            //		                i*veN2[0]/(KB*Tv[1]));
	    		                
            z_vibr_T_Tv_N2 += exp(-veN2[i]/(KB*Tv[1]));
        }

        double z_vibr_T_Tv_O2 = 0.;
        for (int i=0; i<lO2; ++i) {
            //z_vibr_T_Tv_O2 += exp(-(veO2[i]-i*veO2[0])/KBT-
            //		                i*veO2[0]/(KB*Tv[2]));
	    		                
            z_vibr_T_Tv_O2 += exp(-veO2[i]/(KB*Tv[2]));
        }

        std::cout << "Non-equilibrium partition functions ... " << std::endl;
        std::cout << z_vibr_T_Tv_NO << " " 
                  << z_vibr_T_Tv_N2 << " " 	
                  << z_vibr_T_Tv_O2 << std::endl;

        // Compute the vibrational populations according to 
        // Treanor-Marrone or Boltzmann distributions
        double nNO[lNO];
        double nN2[lN2];
        double nO2[lO2];

        std::cout << " Vibrational populations ... " << std::endl;
        for (int i=0; i<lNO; i++) {
	  
	    // Treanor-Marrone distribution	
            //nNO[i] = nd[2] * exp(-(veNO[i]-i*veNO[0])/(KB*T)+i*veNO[0]/(KB*Tv[0]))
            //    / z_vibr_T_Tv_NO;

            // Boltzmann distribution
            nNO[i] = nd[2] * exp(-veNO[i]/(KB*Tv[0])) / z_vibr_T_Tv_NO;
            std::cout << i << " " << nNO[i] << std::endl;
        }

        for (int i=0; i<lN2; i++) {

            //nN2[i] = nd[3] * exp(-(veN2[i]-i*veN2[0])/(KB*T)+i*veN2[0]/(KB*Tv[1]))
            //    / z_vibr_T_Tv_N2;

            nN2[i] = nd[3] * exp(-veN2[i]/(KB*Tv[1])) / z_vibr_T_Tv_N2;
            std::cout << i << " " << nN2[i] << std::endl;
        }

        for (int i=0; i<lO2; i++) {

            //nO2[i] = nd[4] * exp(-(veO2[i]-i*veO2[0])/(KB*T)+i*veO2[0]/(KB*Tv[2]))
            //    / z_vibr_T_Tv_O2;
	        
            nO2[i] = nd[4] * exp(-veO2[i]/(KB*Tv[2])) / z_vibr_T_Tv_O2;
            std::cout << i << " " << nO2[i] << std::endl;
        }

        // Mass and diameter
        double m[5]; // [kg]
        m[0] = 2.3258672171567998e-26; // N
        m[1] = 2.6567628316576e-26;    // O
        m[2] = 4.9826300488143997e-26; // NO
        m[3] = 4.6517344343135997e-26; // N2
        m[4] = 5.3135256633152e-26;    // O2

        double d[5]; // [m]
        d[0] = 3.298e-10;
        d[1] = 2.75e-10;
        d[2] = 3.4061e-10;
        d[3] = 3.40385035355259e-10;
        d[4] = 3.5155e-10;

        // Computation of STS rate coefficients for the exchange reactions
        // with the model of Savelev. We consider 2 reactions:
        // N2(i) + O <=> NO(k) + N
        // O2(i) + N <=> NO(k) + O
        double k_N2_O[lN2][lNO];
        double k_O2_N[lO2][lNO];

        // constants for N2,O2
        double A[2] = {8e-17, 4e-16}; // [m^3 s^-1]
        double b[2] = {0, -0.39};
        double U = 3.*T; // A/(6.*KB);
        double KBU = KB * U; // [J]

        // N2 vibr. energy
        double en2[lN2]    = {};
        double en2_eV[lN2] = {};
        for (int i=0; i<lN2; ++i) {
            en2[i] = veN2[i] + veN2[0];    // [J]
            en2_eV[i] = en2[i] * 6.242e18; // [eV]
            std::cout << " en2_eV: " << i << " " << en2_eV[i] << std::endl;
        }

        // O2 vibr. energy, J
        double eo2[lO2]    = {};
        double eo2_eV[lO2] = {};
        for (int i=0; i<lO2; ++i) {
            eo2[i] = veO2[i] + veO2[0];
            eo2_eV[i] = eo2[i] * 6.242e18; 
            std::cout << " eo2_eV: " << i << " " << eo2_eV[i] << std::endl;
        }

        double eno[lNO]      = {};
        double eno_eV[lNO]   = {};
        double Ear_n2[lNO]   = {};
        double Ear_n2_J[lNO] = {};
        double Ear_o2[lNO]   = {};
        double Ear_o2_J[lNO] = {};
        for (int i=0; i<lNO; ++i) {
            eno[i] = veNO[i] + veNO[0];
            eno_eV[i] = eno[i] * 6.242e18; 
            std::cout << " eno_eV: " << i << " " << eno_eV[i] << std::endl;
        }

        std::cout << " Ear ... " << std::endl;
        for (int i=0; i<lNO; ++i) {

            // N2+O reaction
            Ear_n2[i] = 2.8793 + 1.02227 * eno_eV[i]; // eV
            Ear_n2_J[i] = Ear_n2[i] / 6.242e18;       // J

            // O2+N reaction
            if (eno_eV[i] < 1.3706) {
                Ear_o2[i] = 0.098;
            }
            else if ((1.3706 < eno_eV[i]) && (eno_eV[i] < 2.4121)) {
                Ear_o2[i] = -0.6521 + 0.54736 * eno_eV[i];
            }
            else if (eno_eV[i] > 2.4121) {
                Ear_o2[i] = -1.8451 + 1.04189 * eno_eV[i]; // eV
            } else {
                std::cout << "Something wrong in Ear_o2" << std::endl;
            }
            Ear_o2_J[i] = Ear_o2[i] / 6.242e18; // J

            std::cout << i << " " << Ear_n2_J[i] << "  "
             	  	      << Ear_o2_J[i] << std::endl;
        }

        // Vibration partition functions
        std::cout << "Vibr. partition functions" << std::endl;
        double Zv_n2;
        for (int i=0; i<lN2; ++i)
            Zv_n2 += exp(-en2[i]/KBT);

        double Zv_o2;
        for (int i=0; i<lO2; ++i) 
            Zv_o2 += exp(-eo2[i]/KBT);

        std::cout << " Zv_n2: " << Zv_n2 << " " 
                  << " Zv_o2: " << Zv_o2 << std::endl;

        // Equilibrium coefficients
        std::cout << " N2 eq. coeffs..." << std::endl;
        double k_eq_n2[lNO] = {}; // [m^3 s^-1]
        double fac = A[0] * pow(T,b[0]);
        for (int i=0; i<lNO; ++i) {
            k_eq_n2[i] = fac * (1. + 0.333333333333333 * eno_eV[i]) 
                * exp(-Ear_n2_J[i] / KBT); 
            std::cout << i << " " << k_eq_n2[i] << std::endl;
        }

        std::cout << " O2 eq. coeffs..." << std::endl;
        double k_eq_o2[lNO] = {};
        fac = A[1] * pow(T,b[1]);
        for (int i=0; i<lNO; ++i) {
            k_eq_o2[i] = fac * (Ear_o2[i] + 0.8) * exp(-Ear_o2_J[i] / KBT);
            std::cout << i << " " << k_eq_o2[i] << std::endl;
        }

        // Energy threshold, for each k -> e_i* [-]
        std::cout << "energy threshold" << std::endl;
        double sum1_n2[lNO] = {}; //zeros(lno,1);
        double sum2_n2[lNO] = {}; 
        double sum1_o2[lNO] = {};
        double sum2_o2[lNO] = {};

        // Identify the thresholds to compute the partial sums
        int thrsN2[lNO];
        int thrsO2[lNO];
        for (int i=0; i<lNO; ++i) {

            for (int j=0; j<lN2; ++j) {
                if (en2[j] <= Ear_n2_J[i]) thrsN2[i] = j;
            }

            for (int j=0; j<lO2; ++j) {
                if (eo2[j] <= Ear_o2_J[i]) thrsO2[i] = j;
            }
            std::cout << " thrsN2["<<i<<"] = " << thrsN2[i] 
            	      << " thrsO2["<<i<<"] = " << thrsO2[i] << std::endl;
        }

        std::cout << " Partial sums " << std::endl;
        for (int i=0; i<lNO; ++i) {

            // N2 + O reaction
            for (int j=0; j<thrsN2[i]; ++j) {
                sum1_n2[i] += exp(-(Ear_n2_J[i]-en2[j])/KBU);    
            }
            for (int j=thrsN2[i]; j<lN2; ++j) {
                sum2_n2[i] += exp(-(Ear_n2_J[i]-en2[j])/KBT);    
            }

            // O2 + N reaction
            for (int j=0; j<thrsO2[i]; ++j) {
                sum1_o2[i] += exp(-(Ear_o2_J[i]-eo2[j])/KBU);    
            }
            for (int j=thrsN2[i]; j<lO2; ++j) {
                sum2_o2[i] += exp(-(Ear_o2_J[i]-eo2[j])/KBT);    
            }
            std::cout << i << " " << sum1_n2[i] << " " << sum2_n2[i] << " " 
            	       << sum1_o2[i] << " " << sum2_o2[i] << std::endl;
        }

        // Normalizing coefficients
        std::cout << "Normalizing coeffs" << std::endl;
        double C_n2[lNO]  = {}, C_o2[lNO] = {};
        double B1_n2[lNO] = {}, B2_n2[lNO] = {};
        double B1_o2[lNO] = {}, B2_o2[lNO] = {}; 

        for (int i=0; i<lNO; ++i) {
            C_n2[i] = Zv_n2 * 1./(sum1_n2[i] + sum2_n2[i]);
            B1_n2[i] = C_n2[i] * k_eq_n2[i] * exp(-Ear_n2_J[i]/KBU);
            B2_n2[i] = C_n2[i] * k_eq_n2[i] * exp(Ear_n2_J[i]/KBT);
            std::cout << " C_n2["<<i<<"] = " << C_n2[i] 
            	  << " B1_n2["<<i<<"] = " << B1_n2[i] 
            	  << " B2_n2["<<i<<"] = " << B2_n2[i] << std::endl;
        }

        for (int i=0; i<lNO; ++i) {
            C_o2[i] = Zv_o2 * 1./(sum1_o2[i] + sum2_o2[i]);
            B1_o2[i] = C_o2[i] * k_eq_o2[i] * exp(-Ear_o2_J[i]/KBU);
            B2_o2[i] = C_o2[i] * k_eq_o2[i] * exp(Ear_o2_J[i]/KBT);
            std::cout << " C_o2["<<i<<"] = " << C_o2[i] 
            	  << " B1_o2["<<i<<"] = " << B1_o2[i] 
            	  << " B2_o2["<<i<<"] = " << B2_o2[i] << std::endl;
        }

        std::cout << " k_N2_O " << std::endl;
        double kbTU =  (1./T + 1./U) / KB;
        for (int i=0; i<lN2; ++i) { // N2 rows
            for (int j = 0; j < lNO; ++j) { // NO columns
                if (en2[i] < Ear_n2_J[j]) { // Ear_n2_J has size = lNO 
                    k_N2_O[i][j] = B1_n2[i] * exp(en2[i]/kbTU);
                }
                else if (en2[i] > Ear_n2_J[j]) {
            	k_N2_O[i][j] = B2_n2[i];
                } else {
                    std::cout << " Out-of-bound error ... " << std::endl;
                }
                std::cout << " k_N2_O["<<i<<"]["<<j<<"] = " 
            	      << k_N2_O[i][j] << std::endl;
            }
        }

        std::cout << " k_O2_N " << std::endl;
        for (int i=0; i<lO2; ++i) { // O2 rows 
            for (int j = 0; j < lNO; ++j) { // NO columns
                if (eo2[i] < Ear_o2_J[j]) { // Ear_o2_J has size = lNO
                    k_O2_N[i][j] = B1_o2[i] * exp(eo2[i]/kbTU);
                }
                else if (eo2[i] > Ear_o2_J[j]) {
            	k_O2_N[i][j] = B2_o2[i];
                } else {
                    std::cout << " Out-of-bound error ... " << std::endl;
                }
                std::cout << " k_O2_N["<<i<<"]["<<j<<"] = " 
            	      << k_O2_N[i][j] << std::endl;
            }
        }

        /*
        ////////////////////////////////////////////////////////////////////////////
        // 1. Rusanov-Friedman model
        for (int i=0; i<lN2; ++i) {
            if (Ear_n2_J[i] > 0.51 * en2[i] ) {
              B2_n2[i] = A[0]*pow(T, b[0])*exp(-(Ear_n2_J[i] - 0.51*en2[i])/(KB*T));
            } else {
              B2_n2[i] = A[0] * pow(T, b[0]);
            }
        }
        ////////////////////////////////////////////////////////////////////////////
        // 2. Polak model
        for (int i=0; i<lN2; ++i) {
            if (Ear_n2_J[i] > 0.52 * en2[i] ) {
              B2_n2[i] = A[0] * pow(T, b[0]) *
                  exp(-(Ear_n2_J[i] - 0.52 * en2[i]) / (0.9 * KB * T));
            } else {
              B2_n2[i] = A[0] * pow(T, b[0]);
            }
        }
        ////////////////////////////////////////////////////////////////////////////
        // 3. Warnatz model
        for (int i=0; i<lN2; ++i) {
            if (Ear_n2_J[i] > en2[i] ) {
              B2_n2[i] = (4.17e+12 + i) * pow(T, 0) *
                  exp(-(Ear_n2_J[i] - en2[i]) / (KB * T));
            } else {
              B2_n2[i] = (4.17e+12 + i) * pow(T, 0);
            }
        }
        ////////////////////////////////////////////////////////////////////////////
        */

        // Compute the rotational partition function
        // according to a simplified formula ... TODO: better formula
        double Be[3] = {167.20, 199.8, 143.77}; // [m^-1]
        int sigma[3] = {1, 2, 2};
        double theta_r[3] = {};
        double z_rot[3] = {};
        for (int i=0; i<3; ++i) { // NO N2 O2
            theta_r[i] = Be[i]*HP*C0/KB; // [K]
            z_rot[i] = T/(sigma[i]*theta_r[i]);
        }

        std::cout << " Kz_n2 ... " << std::endl;
        double Kz_n2[lN2][lNO] = {}; 
        double fac0 = pow((m[3]*m[1]/(m[2]*m[0])),1.5) * z_rot[1]/z_rot[0] 
                        * exp((de[1]-de[0])/KBT);
        for (int i=0; i<lN2; ++i) {
            for (int j=0; j<lNO; ++j) {
                Kz_n2[i][j] = fac0 * exp((-veN2[i] + veNO[j])/KBT);
		std::cout << Kz_n2[i][j] << " \n"[j == lNO-1];
            }
        }

        std::cout << "Kz_o2" << std::endl;
        double Kz_o2[lO2][lNO] = {}; 
        double fac1 = pow((m[4]*m[0]/(m[2]*m[1])),1.5) * z_rot[2]/z_rot[0] 
                        * exp((de[2]-de[0])/KBT);
        for (int i=0; i<lO2; ++i) {
            for (int j=0; j<lNO; ++j) {
                Kz_o2[i][j] = fac1 * exp((-veNO[j] + veO2[i])/KBT);
		std::cout << Kz_o2[i][j] << " \n"[j == lNO-1];
            }
        }

        // Compute MT forward exchange coefficients, eq. 3.93.
        // N2
        double k_O_N2[lN2][lNO];
        //double fac0f, fac0b = 0.;
        for (int i=0; i<lN2; ++i) {
            for (int j=0; j<lNO; ++j) {

                // MT forward rates
                //fac0f += exp(-((veN2[i]-i*veN2[0])/KBT) -
                //                        i*veN2[0]/(KB*Tv[1])) * k_N2_O[i][j];
                // STS backward rates
                k_O_N2[i][j] = Kz_n2[i][j] * k_N2_O[i][j];
                //fac0b += exp(-((veN2[i]-i*veN2[0])/KBT) -
                //                        i*veN2[0]/(KB*Tv[1])) * k_O_N2[i][j];
            }
        }
        //double k_exch_N2_NO = fac0f / z_vibr_T_Tv_N2;
        //double k_exch_NO_N2 = fac0b / z_vibr_T_Tv_N2;
        //std::cout << " k_exch_N2_NO: " << k_exch_N2_NO << std::endl; 
        //std::cout << " k_exch_NO_N2: " << k_exch_NO_N2 << std::endl; 

        // O2
        double k_N_O2[lO2][lNO];
        for (int i=0; i<lO2; ++i) {
            //fac0f = fac0b = 0.;
            for (int j=0; j<lNO; ++j) {

                // MT forward rates
                //fac0f += exp(-((veO2[i]-i*veO2[0])/KBT) -
                //                        i*veO2[0]/(KB*Tv[2])) * k_O2_N[i][j];
                // STS backward rates
                k_N_O2[i][j] = Kz_o2[i][j] * k_O2_N[i][j];
                //fac0b += exp(-((veO2[i]-i*veO2[0])/KBT) -
                //                        i*veO2[0]/(KB*Tv[2])) * k_N_O2[i][j];
            }
        }
        //double k_exch_O2_NO = fac0f / z_vibr_T_Tv_O2;
        //double k_exch_NO_O2 = fac0b / z_vibr_T_Tv_O2;
        //std::cout << " k_exch_O2_NO: " << k_exch_O2_NO << std::endl; 
        //std::cout << " k_exch_NO_O2: " << k_exch_NO_O2 << std::endl; 

        //double R_exch_22[m_mixture.nMolecules()] = {};
	//Map<ArrayXd>(mp_r_w_22, m_mixture.nEnergyEqns()-1) = 0.;

        // NO
	fac0, fac1 = 0.;
	for (int i=0; i<lNO; ++i) {
	    for (int j=0; j<lN2; ++j) {
		fac0 += (nN2[j] * /*nO=*/nd[1] * k_O_N2[i][j] - 
			 nNO[i] * /*nN=*/nd[0] * k_N2_O[i][j]);
	    }
	    for (int j=0; j<lO2; ++j) {
		fac1 += (nO2[j] * /*nN=*/nd[0] * k_N_O2[i][j] - 
			 nNO[i] * /*nO=*/nd[1] * k_O2_N[i][j]);
	    }
	    //R_exch_22[0] += i * (fac0 + fac1);
	    mp_r_w_22[0] += i * (fac0 + fac1);
	}

        // N2
	fac0, fac1 = 0.;
	for (int i=0; i<lN2; ++i) {
	    for (int j=0; j<lNO; ++j) {
		fac0 += (nNO[j] * /*nN=*/nd[0] * k_O_N2[i][j] - 
		         nN2[i] * /*nO=*/nd[1] * k_N2_O[i][j]);
	    }
	    //R_exch_22[1] += i * fac0;
	    mp_r_w_22[1] += i * fac0;
	}

        // O2
	fac0, fac1 = 0.;
	for (int i=0; i<lO2; ++i) {
	    for (int j=0; j<lNO; ++j) {
		fac0 += (nNO[j] * /*nO=*/nd[1] * k_N_O2[i][j] - 
			 nO2[i] * /*nN=*/nd[0] * k_O2_N[i][j]);
	    }
	    //R_exch_22[2] += i * fac0;
	    mp_r_w_22[2] += i * fac0;
	}

        for (int i=0; i<m_mixture.nEnergyEqns()-1; ++i) {
            //mp_r_w_22[i] = R_exch_22[i];
            std::cout << " mp_r_w_22["<<i<<"] = " << mp_r_w_22[i] << std::endl;
        }
    }

private:

    int m_ns;
    int m_transfer_offset;
    double* mp_wrk1;
    double* mp_wrk2;
    double* mp_Mw;
    double* mp_r_w_22;

};

// Register the transfer model
Utilities::Config::ObjectProvider<R_w_react_22, TransferModel> 
    r_w_react_22("R_w_react_22");

    } // namespace Transfer
} // namespace Mutation 
