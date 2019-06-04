/**
 * @file R_vibr_VT.cpp
 *
 * @brief Implementation of R_vibr_VT.
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

#include <vector>
#include <cmath>
#include <math.h>

#include "Mixture.h"
#include "TransferModel.h"
#include "Numerics.h"
#include "VSS.h"
#include "FHO.h"

using namespace Mutation;

namespace Mutation {
    namespace Transfer {

/**
 * Represents a coupling between vibrational and translational energy modes.
 */
class R_vibr_VT : public TransferModel
{
public:

    R_vibr_VT(Mixture& mix)
        : TransferModel(mix), m_vss(mix), m_fho(mix)
    {
        m_ns              = m_mixture.nSpecies();
        m_transfer_offset = m_mixture.hasElectrons() ? 1 : 0;

        mp_rVT = new double [m_mixture.nEnergyEqns()-1];
    }

    virtual ~R_vibr_VT()
    {
        delete [] mp_rVT;
    }

    /**
     * Computes the source terms of the Vibration-Translational energy transfer
     * in \f$ [J/(m^3\cdot s)] \f$ using a Landau-Teller formula.
     *
     * \f[ 
     * \R_vibr_{VT}_m = \rho_m \frac{e^V_m\left(T\right) - 
     * 		 	             e^V_m\left(T_{vm}\right)} {\tau^{VT}_m} 
     * \f]
     *
     * More information about the above model can be found in 
     * @cite E. V. Kustova and G. P. Oblapenko 2015, 2016 and
     * @cite G. P. Oblapenko 2018.
     *
     */
    double source() {};
    void R_w_22(double* const p_r_w_22) {};
    void R_w_23(double* const p_r_w_23) {};
    void rVV(double* const p_rVv) {};

    void rVT(double* mp_rVT)
    {
        const double * p_X = m_mixture.X();
        //const double * p_Y = m_mixture.Y();
        double n = m_mixture.numberDensity();
        double T = m_mixture.T();
        double Tv[m_mixture.nMolecules()] = {};
	for (int i=0; i<m_mixture.nMolecules(); ++i)
	     Tv[i] = m_mixture.Tvs(i);

        int inv = 0; // non-vibrators
        for (int iv = 0; iv-inv < m_vss.nVibrators(); ++iv){
	    // If we write only molecular vibrators in the VSS file,
	    // the if branch is no more needed.
            //if(m_mixture.species(iv).type() != Mutation::Thermodynamics::MOLECULE){
            //    inv++;
            //} else { 
		mp_rVT[iv] = 4. * KB * (T/Tv[iv]) * (T-Tv[iv]) * (p_X[iv] * n)
		    * compute_ave_op(iv-inv);
	
	    // Alternatively, we could compute the source term, by first
	    // calculating the vibrational relaxation time tauVT ... TODO
            //}
        }
    }

private:

    VSS m_vss;
    FHO m_fho;

    /**
     * @brief Computes the averaging operator over heavy particles
     * pre-multiplied by the sum of the number density of all the
     * collisional partners.
     *
     * @param vibrator index
     *
     * @return
     */
    double compute_ave_op(int const);

    inline double VT_integral(double, double,
        std::vector<double> vibr_spectrum, double beta, 
	    double dref, double omega, double coll_mass, 
	        double osc_mass, double ram1, double ram2, 
		    double E_Morse, double Tref=273., 
		         double this_svt=0.5);

    inline double Zv(double, std::vector<double>);

    inline double vt_rate_fho(double T, double beta, double dref, 
	double omega, double coll_mass, double osc_mass, double ve_before, 
	    double ve_after, int i, int delta, double ram1, double ram2, 
	        double E_Morse, double Tref=273., double this_svt=0.5);

    inline double vt_prob_g_only_fho(double g, double mass, double beta, 
        double osc_mass, double ve_before, double ve_after, int i, int delta, 
	    double ram1, double ram2, double E_Morse, double this_svt);

    inline double vt_prob_g_only_fho_12(double g, double mass, double beta, 
        double osc_mass, double ve_before, double ve_after, int i, int delta, 
	    double ram, double E_Morse, double this_svt);

    inline double vel_avg_vt(double g, double ve_before, double ve_after, 
        double mass);

    inline double c_vibr(double Tv, std::vector<double> vibr_spectrum, 
        double mass);

    inline double cs_vss(double g, double dref, double gref, double omega);

    int m_ns;
    int m_transfer_offset;
    double* mp_rVT;
};
 
////////////////////////////////////////////////////////////////////////////////

// Compute the sum of integrals over VT transition cross-section
inline double R_vibr_VT::VT_integral(double T, double Tv,
    std::vector<double> vibr_spectrum, double beta, double dref, 
        double omega, double coll_mass, double osc_mass, double ram1, 
	    double ram2, double E_Morse, double Tref, double this_svt)
{
    double res, tmp, vtr = 0.;
    double dEsq = pow(((vibr_spectrum.at(1) - vibr_spectrum.at(0))/ (KB*T)), 2.);
    double rev_k_mult = exp((vibr_spectrum[0] - vibr_spectrum[1]) / (KB*T));
    double kTv = KB * Tv;
    for (int i=0; i<vibr_spectrum.size()-1; ++i) {

        vtr = vt_rate_fho(T, beta, dref, omega, coll_mass, osc_mass,
            vibr_spectrum[i+1], vibr_spectrum[i], i+1, -1, ram1, ram2, 
	        E_Morse, Tref, this_svt);

        tmp = dEsq * exp(-vibr_spectrum[i+1] / kTv);
        tmp *= vtr;
        res += tmp;

        // v -> v + 1
        tmp = dEsq * exp(-vibr_spectrum[i] / kTv);
        tmp *= vtr * rev_k_mult;
        res += tmp;
    }
    return res / (Zv(Tv, vibr_spectrum) * 8.);
}

///////////////////////////////////////////////////////////////////////////////////

inline double R_vibr_VT::c_vibr(double Tv, std::vector<double> vibr_spectrum, 
    double mass)
{
    double vae = 0.;
    for (int i=0; i<vibr_spectrum.size(); ++i)
        vae += vibr_spectrum[i] * exp(-vibr_spectrum[i] / (KB * Tv));
    vae /= Zv(Tv, vibr_spectrum);

    double vaesq = 0.;
    for (int i=0; i<vibr_spectrum.size(); ++i)
        vaesq += (vibr_spectrum[i] * vibr_spectrum[i] * exp(-vibr_spectrum[i]
            / (KB * Tv)));
    vaesq /= Zv(Tv, vibr_spectrum);

    return (vaesq - vae * vae) / (KB * Tv * Tv * mass);
}

//////////////////////////////////////////////////////////////////////////////////

// VSS cross-section calculation 
// g: velocity, dref and gref are reference parameters
inline double R_vibr_VT:: cs_vss(double g, double dref, double gref, double omega)
{
    return PI * dref * dref * pow((g / gref), (1. - 2. * omega))
        / tgamma(2.5 - omega);
}

//////////////////////////////////////////////////////////////////////////////////

inline double R_vibr_VT::Zv(double Tv, std::vector<double> vibr_spectrum)
{
    double zv = 0.;
    for (int i=0; i<vibr_spectrum.size(); ++i) 
       zv += std::exp(-vibr_spectrum[i] / (KB * Tv));
    return zv;
}

//////////////////////////////////////////////////////////////////////////////////    

// Compute VT transition rate
inline double R_vibr_VT::vt_rate_fho(double T, double beta, double dref, 
    double omega, double coll_mass, double osc_mass, double ve_before, 
        double ve_after, int i, int delta, double ram1, double ram2, 
	    double E_Morse, double Tref, double this_svt)
{
    double mult = 0.;	
    if (delta == 1)
        mult = (i + 1);
    else if (delta == -1)
        mult = i;
    else if (delta > 0.)
        mult = Numerics::fact_div_fact(i, i + delta) 
	    / pow(Numerics::factorial(delta), 2.);
    else
        mult = Numerics::fact_div_fact(i + delta, i) 
	    / pow(Numerics::factorial(-delta), 2.);
	
    double kT = KB *T;
    mult *= sqrt(kT / (2. * PI * coll_mass));
    double gref = sqrt(2. * KB * Tref / coll_mass);

    double min_g = 0.;
    if (ve_after <= ve_before)
        min_g = 0.;
    else
        min_g = sqrt((ve_after - ve_before) / kT);

    // capture all by reference
    auto f = [&] (double g) 
    {
        return vt_prob_g_only_fho(g * sqrt(2. * kT / coll_mass), coll_mass, 
                   beta, osc_mass, ve_before, ve_after, i, delta, ram1, ram2, 
        	       E_Morse, this_svt) * 
               cs_vss(g * sqrt(2. * kT / coll_mass), dref, gref, omega) * 
               g*g*g * std::exp(-g*g);
    };

    return 8. * mult * Numerics::integrate_semi_inf(f, min_g);
}

//////////////////////////////////////////////////////////////////////////////////    

// Compute velocity-dependent part of VT transition probability, 
// for heteronuclear molecules, compute 2 probabilities and take the average
inline double R_vibr_VT::vt_prob_g_only_fho(double g, double mass, double beta, 
    double osc_mass, double ve_before, double ve_after, int i, int delta, 
        double ram1, double ram2, double E_Morse, double this_svt)
{
    double res = 0.;
    if (ram1 == ram2) {

        return vt_prob_g_only_fho_12(g, mass, beta, osc_mass, ve_before, 
	    ve_after, i, delta, ram1, E_Morse, this_svt);

    } else {

        res = vt_prob_g_only_fho_12(g, mass, beta, osc_mass, ve_before, 
	    ve_after, i, delta, ram1, E_Morse, this_svt);

        res += vt_prob_g_only_fho_12(g, mass, beta, osc_mass, ve_before, 
	    ve_after, i, delta, ram2, E_Morse, this_svt);
    }

    return res / 2.;
}

/////////////////////////////////////////////////////////////////////////////////

// Compute the VT transition probability
// Warning: there is an additional multiplier missing here!
// It is moved outside to vt_rate, so that we do less FLOPs inside the 
// integration routine. 
inline double R_vibr_VT::vt_prob_g_only_fho_12(double g, double mass, 
    double beta, double osc_mass, double ve_before, double ve_after, int i, 
        int delta, double ram, double E_Morse, double this_svt)
{
    double res = 0.;
    double omega = 0.;

    double vel = vel_avg_vt(g, ve_before, ve_after, mass);

    if (delta == 1)
        omega = (ve_after - ve_before) / HPBAR;
    else if (delta == -1)
        omega = (ve_before - ve_after) / HPBAR;
    else 
	omega = 0.;

    double eps = 1.;
    double phi = (2. / PI) * std::atan(std::sqrt((2. * E_Morse) 
        / (mass * vel * vel)));

    eps *= std::pow((cosh((1. + phi) * PI * omega / (beta * vel))), 2.);
    eps *= 8. * ram * ram / std::pow((std::sinh(2 * PI * omega / (beta * vel))), 2.);
    eps *= this_svt * (PI * PI) * omega * mass * mass 
        / (osc_mass * beta * beta * HP);

    if (delta == 1)
        res = eps * std::exp(-(i + 1) * eps);
    else if (delta == -1)
        res = eps * std::exp(-i * eps);

    return res;
}

/////////////////////////////////////////////////////////////////////////////////    

// Average the velocities before and after a collision
inline double R_vibr_VT::vel_avg_vt(double g, double ve_before, double ve_after, 
    double mass)
{
    double gn_sq = (ve_before - ve_after) * (2. / mass) + (g * g);
    if (gn_sq < 0.)
      return -1.;
    else
      return 0.5 * (g + sqrt(gn_sq));
}

/////////////////////////////////////////////////////////////////////////////////

double R_vibr_VT::compute_ave_op(int i_vibrator)
{
    const double T = m_mixture.T();
    const double n = m_mixture.numberDensity();
    const double * p_X = m_mixture.X();
    double Tv[m_mixture.nMolecules()] = {};
    for (int i=0; i<m_mixture.nMolecules(); ++i) 
        Tv[i] = m_mixture.Tvs(i);
    double osc_mass = m_vss[i_vibrator].osc_mass();
    
    // Harmonic vibrational spectrum of NO N2 O2
    std::vector<std::vector<double> > vibr_spectrum {
    {0.00000000e+00, 3.78297694e-20, 7.56595389e-20, 1.13489308e-19, // NO
     1.51319078e-19, 1.89148847e-19, 2.26978617e-19, 2.64808386e-19,
     3.02638156e-19, 3.40467925e-19, 3.78297694e-19, 4.16127464e-19,
     4.53957233e-19, 4.91787003e-19, 5.29616772e-19, 5.67446542e-19,
     6.05276311e-19, 6.43106081e-19, 6.80935850e-19, 7.18765620e-19,
     7.56595389e-19, 7.94425158e-19, 8.32254928e-19, 8.70084697e-19,
     9.07914467e-19, 9.45744236e-19, 9.83574006e-19, 1.02140378e-18},
    {0.00000000e+00, 4.68454043e-20, 9.36908086e-20, 1.40536213e-19, // N2
     1.87381617e-19, 2.34227021e-19, 2.81072426e-19, 3.27917830e-19,
     3.74763234e-19, 4.21608639e-19, 4.68454043e-19, 5.15299447e-19,
     5.62144851e-19, 6.08990256e-19, 6.55835660e-19, 7.02681064e-19,
     7.49526469e-19, 7.96371873e-19, 8.43217277e-19, 8.90062681e-19,
     9.36908086e-19, 9.83753490e-19, 1.03059889e-18, 1.07744430e-18,
     1.12428970e-18, 1.17113511e-18, 1.21798051e-18, 1.26482592e-18,
     1.31167132e-18, 1.35851672e-18, 1.40536213e-18, 1.45220753e-18,
     1.49905294e-18},
    {0.00000000e+00, 3.13821409e-20, 6.27642817e-20, 9.41464226e-20, // O2
     1.25528563e-19, 1.56910704e-19, 1.88292845e-19, 2.19674986e-19,
     2.51057127e-19, 2.82439268e-19, 3.13821409e-19, 3.45203549e-19,
     3.76585690e-19, 4.07967831e-19, 4.39349972e-19, 4.70732113e-19,
     5.02114254e-19, 5.33496395e-19, 5.64878535e-19, 5.96260676e-19,
     6.27642817e-19, 6.59024958e-19, 6.90407099e-19, 7.21789240e-19,
     7.53171381e-19, 7.84553521e-19}};

    double ave = 0.;
    double dref, omega, mu = 0.;
    double beta, E_Morse, svt = 0.;

    // Partner offset
    for (int i_partner = m_transfer_offset; i_partner < m_ns; ++i_partner){

	// VSS stuff
	dref = m_vss[i_vibrator][i_partner].dref();
	omega = m_vss[i_vibrator][i_partner].omega();
	mu = m_vss[i_vibrator][i_partner].mu();

	// FHO stuff
	beta = m_fho[i_vibrator][i_partner].beta();
	E_Morse = m_fho[i_vibrator][i_partner].E_Morse();
	svt = m_fho[i_vibrator][i_partner].svt();

	// TODO: replace 0.5, 0.5 with computed reduced masses ...
        ave = VT_integral(T, Tv[i_vibrator], vibr_spectrum[i_vibrator],
	    beta, dref, omega, mu, osc_mass, 0.5, 0.5, E_Morse, 273., svt);

        ave += (p_X[i_partner] * n) * ave;
    }
    return ave;
}

/////////////////////////////////////////////////////////////////////////////////

// Register the transfer model
Utilities::Config::ObjectProvider<
    R_vibr_VT, TransferModel> r_vibr_VT("R_vibr_VT");

    } // namespace Transfer
} // namespace Mutation 
