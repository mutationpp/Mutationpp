/**
 * @file RrahoDB.cpp
 *
 * @brief Provides a Rigid-Rotator anHarmonic Oscillator thermodynamic database.
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

#include "Constants.h"
#include "ThermoDB.h"
#include "Species.h"
#include "ParticleRRaHO.h"
#include "AutoRegistration.h"
#include "Functors.h"
#include "LookupTable.h"
#include "Utilities.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cassert>

using namespace std;
using namespace Mutation::Numerics;
using namespace Mutation::Utilities;

namespace Mutation {
    namespace Thermodynamics {

// Simple for loop over the species just to make life a little easier
#define LOOP(__op__)\
for (int i = 0; i < m_ns; ++i) {\
    __op__ ;\
}

#define LOOP_ATOMS(__op__)\
for (int i = 0, j = 0; i < m_na; ++i) {\
    j = mp_indices[i];\
    __op__ ;\
}

// Loops over molecules. Inside the loop, the index i is zero based and indexes
// internal molecular data.  The index j is the index corresponding to the 
// original species data.
#define LOOP_MOLECULES(__op__)\
for (int i = 0, j = 0; i < m_nm; ++i) {\
    j = mp_indices[m_na+i];\
    __op__ ;\
}

// Loops over heavy particles (non electron species).  Inside the loop, index i
// is zero based and indexes internal molecular data. The index j is the index
// corresponding to the original species data.
#define LOOP_HEAVY(__op__)\
for (int i = 0, j = 0; i < m_na + m_nm; ++i) {\
    j = (m_has_electron ? i+1 : i);\
    __op__ ;\
}

typedef struct {
    double ln_omega_t;  // ln(omega^(2 / L) * theta)
    double linearity;   // L / 2
    int sym;		// symmetry/steric factor
} RotData;

typedef struct {
    double g;           // degeneracy
    double theta;       // characteristic temperature
    double energy;      // electronic energy
    double we;          // spectroscopic constant
} ElecLevel;

typedef struct {
    unsigned int offset;
    unsigned int nheavy;
    unsigned int nlevels;
    int* p_nelec;
    ElecLevel* p_levels;
} ElectronicData;


/**
 * A thermodynamic database that uses the Rigid-Rotator anHarmonic-Oscillator
 * model for computing species thermodynamic properties. See the individual 
 * thermodynamic functions for specific descriptions of the model.
 */
class RrahoDB : public ThermoDB
{
public:

    RrahoDB(int arg)
        : ThermoDB(298.15, 101325.0), m_ns(0), m_na(0), m_nm(0),
          m_has_electron(false),
          m_use_tables(true),
          m_last_bfacs_T(0.0)
    { }
    
    /**
     * Destructor.
     */
    ~RrahoDB() 
    {
        delete [] mp_lnqtmw;
        delete [] mp_hform;
        delete [] mp_eform;
        delete [] mp_indices;
        delete [] mp_rot_data;
        delete [] mp_nvib;
        delete [] mp_vib_temps;
        
        delete [] m_elec_data.p_nelec;
        delete [] m_elec_data.p_levels;
        delete [] mp_part_sst;
        delete [] mp_el_bfacs;
        
        if (m_use_tables) {
            delete mp_el_bfac_table;
        }
    }

    /**
     * Computes the unitless species specific heat at constant pressure
     * \f$ C_{P,i} / R_U\f$ in thermal nonequilibrium.
     *
     * @param Th  - heavy particle translational temperature
     * @param Te  - free electron temperature
     * @param Tr  - mixture rotational temperature
     * @param Tv  - mixture vibrational temperature
     * @param Tel - mixture electronic temperature
     * @param cp   - on return, the array of species non-dimensional \f$C_P\f$'s
     * @param cpt  - if not NULL, the array of species translational \f$C_P\f$'s
     * @param cpr  - if not NULL, the array of species rotational \f$C_P\f$'s
     * @param cpv  - if not NULL, the array of species vibrational \f$C_P\f$'s
     * @param cpel - if not NULL, the array of species electronic \f$C_P\f$'s
     *
     * @todo Compute \f$C_P\f$ directly instead of using finite-differencs.
     */
    void cp(
        double Th, double Te, double Tr, double Tv, double Tel, 
        double* const cp, double* const cpt, double* const cpr, 
        double* const cpv, double* const cpel)
    {
        // Special case if we only want total Cp
        if (cp != NULL && cpt == NULL && cpr == NULL && cpv == NULL && 
            cpel == NULL)
        {
            cpT(cp, Eq());
            cpR(cp, PlusEq());
            cpV(Tv, cp, PlusEq());
            cpE(Tel, cp, PlusEq());
            return;
        }
        
        // Otherwise we have to compute each component directly.
        // Translation
        if (cpt == NULL) {
            if (cp != NULL)
                //cpT(cp, PlusEq());
                cpT(cp, Eq());
        } else {
            cpT(cpt, Eq());
            if (cp != NULL) 
                LOOP(cp[i] = cpt[i]);
        }
        
        // Rotation
        if (cpr == NULL) {
            if (cp != NULL)
                cpR(cp, PlusEq());
        } else {
            cpR(cpr, Eq());
            if (cp != NULL) 
                LOOP_MOLECULES(cp[j] += cpr[j]);
        }
        
        // Vibration
        if (cpv == NULL) {
            if (cp != NULL)
                cpV(Tv, cp, PlusEq());
        } else {
            cpV(Tv, cpv, Eq());
            if (cp != NULL) 
                LOOP_MOLECULES(cp[j] += cpv[j]);
        }
        
        // Electronic
        if (cpel == NULL) {
            if (cp != NULL)
                cpE(Tel, cp, PlusEq());
        } else {
            cpE(Tel, cpel, Eq());
            if (cp != NULL)
                LOOP_HEAVY(cp[j] += cpel[j]);
        }
    }
    
//    /**
//     * Computes the species vibrational specific heats at the given temperature
//     * nondimensionalized by Ru.
//     */
//    void cpv(double Tv, double* const p_cpv)
//    {
//        int ilevel = 0;
//        double sum, fac1, fac2;
//        LOOP_MOLECULES(
//            sum = 0.0;
//            for (int k = 0; k < mp_nvib[i]; ++k, ilevel++) {
//                fac1 = mp_vib_temps[ilevel] / Tv;
//                fac2 = std::exp(fac1);
//                fac1 = fac1 / (fac2 - 1.0);
//                sum += fac2*fac1;
//            }
//            p_cpv[j] = sum;
//        )
//    }
//
//    void cpel(double Tel, double* const p_cpel)
//    {
//        CpelFunctor()(Tel, p_cpel, m_elec_data, Eq());
//    }

    /**
     * Computes the unitless species specific heat at constant volume
     * \f$ C_{V,i} / R_U\f$ in thermal nonequilibrium.
     *
     * @param Th  - heavy particle translational temperature
     * @param Te  - free electron temperature
     * @param Tr  - mixture rotational temperature
     * @param Tv  - mixture vibrational temperature
     * @param Tel - mixture electronic temperature
     * @param cv   - on return, the array of species non-dimensional \f$C_V\f$'s
     * @param cvt  - if not NULL, the array of species translational \f$C_V\f$'s
     * @param cvr  - if not NULL, the array of species rotational \f$C_V\f$'s
     * @param cvv  - if not NULL, the array of species vibrational \f$C_V\f$'s
     * @param cvel - if not NULL, the array of species electronic \f$C_V\f$'s
     *
     * @todo Compute \f$C_V\f$ directly instead of using finite-differencs.
     */
    void cv(
        double Th, double Te, double Tr, double Tv, double Tel, 
        double* const cv, double* const cvt, double* const cvr, 
        double* const cvv, double* const cvel)
    {
        // Special case if we only want total Cv
        if (cv != NULL && cvt == NULL && cvr == NULL && cvv == NULL && 
            cvel == NULL)
        {
            cvT(cv, Eq());
            cvR(cv, PlusEq());
	    //cvR(Tr, cv, PlusEq());
            cvV(Tv, cv, PlusEq());
            cvE(Tel, cv, PlusEq());
            return;
        }
        
        // Otherwise we have to compute each component directly.
        // Translation
        if (cvt == NULL) {
            if (cv != NULL)
                cvT(cv, Eq());
        } else {
            cvT(cvt, Eq());
            if (cv != NULL) 
                LOOP(cv[i] = cvt[i]);
        }
        
        // Rotation
        if (cvr == NULL) {
            if (cv != NULL)
                cvR(cv, PlusEq());
	        //cvR(Tr, cv, PlusEq());
        } else {
            cvR(cvr, Eq());
	    //cvR(Tr, cvr, Eq());
            if (cv != NULL) 
                LOOP_MOLECULES(cv[j] += cvr[j]);
        }
        
        // Vibration
        if (cvv == NULL) {
            if (cv != NULL)
                cvV(Tv, cv, PlusEq());
        } else {
            cvV(Tv, cvv, Eq());
            if (cv != NULL) 
                LOOP_MOLECULES(cv[j] += cvv[j]);
        }
        
        // Electronic
        if (cvel == NULL) {
            if (cv != NULL)
                cvE(Tel, cv, PlusEq());
        } else {
            cvE(Tel, cvel, Eq());
            if (cv != NULL)
                LOOP_HEAVY(cv[j] += cvel[j]);
        }
    }
    
//    /** TODO
//     * Computes the species vibrational specific heats at the given temperature
//     * nondimensionalized by Ru.
//     */
//    void cvv(double Tv, double* const p_cvv)
//    {
//    //    int ilevel = 0;
//    //    double sum, fac1, fac2;
//    //    LOOP_MOLECULES(
//    //        sum = 0.0;
//    //        for (int k = 0; k < mp_nvib[i]; ++k, ilevel++) {
//    //            fac1 = mp_vib_temps[ilevel] / Tv;
//    //            fac2 = std::exp(fac1);
//    //            fac1 = fac1 / (fac2 - 1.0);
//    //            sum += fac2*fac1;
//    //        }
//    //        p_cvv[j] = sum;
//    //    )
//    }
//  
//    // TODO
//    void cvel(double Tel, double* const p_cvel)
//    {
//        CvelFunctor()(Tel, p_cvel, m_elec_data, Eq());
//    }

    /**
     * Computes the unitless species enthalpy \f$ h_i/R_U T_h\f$ of each
     * species in thermal nonequilibrium, which is non-dimensionalized by the 
     * heavy particle translational temperature. 
     *
     * @param Th  - heavy particle translational temperature
     * @param Te  - free electron temperature
     * @param Tr  - mixture rotational temperature
     * @param Tv  - mixture vibrational temperature
     * @param Tel - mixture electronic temperature
     * @param h   - on return, the array of species non-dimensional enthalpies
     * @param ht  - if not NULL, the array of species translational enthalpies
     * @param hr  - if not NULL, the array of species rotational enthalpies
     * @param hv  - if not NULL, the array of species vibrational enthalpies
     * @param hel - if not NULL, the array of species electronic enthalpies
     * @param hf  - if not NULL, the array of the species formation enthalpies
     */
    void enthalpy(
        double Th, double Te, double Tr, double Tv, double Tel, double* const h,
        double* const ht, double* const hr, double* const hv, double* const hel,
        double* const hf)
    {
        
        // Special case where we only want the total enthalpy
        if (ht == NULL && hr == NULL && hv == NULL && hel == NULL && 
            hf == NULL && h != NULL) 
        {
            hT(Th, Te, h, Eq());
            hR(Tr, h, PlusEq());
            hV(Tv, h, PlusEq());
            hE(Tel, h, PlusEq());
            hF(h, PlusEq());
            LOOP(h[i] /= Th);
            return;
        }
        
        // Otherwise selectively choose what we want
        // Translational enthalpy
        if (ht == NULL) {
            //hT(Th, Te, h, Eq());
            if (h != NULL)
                hT(Th, Te, h, EqDiv(Th));
        } else {
            hT(Th, Te, ht, EqDiv(Th));
            if (h != NULL)
                LOOP(h[i] = ht[i]);
        }
        
        // Rotatonal enthalpy
        if (hr == NULL) {
            if (h != NULL)
                hR(Tr, h, PlusEqDiv(Th));
        } else {
            LOOP(hr[i] = 0.0);
            hR(Tr, hr, EqDiv(Th));
            if (h != NULL)
                LOOP_MOLECULES(h[j] += hr[j]);
        }
        
        // Vibrational enthalpy
        if (hv == NULL) {
            if (h != NULL)
                hV(Tv, h, PlusEqDiv(Th));
        } else {
            LOOP(hv[i] = 0.0);
            hV(Tv, hv, EqDiv(Th));
            if (h != NULL)
                LOOP_MOLECULES(h[j] += hv[j]);
        }
            
        // Electronic enthalpy
        if (hel == NULL) {
            if (h != NULL)
                hE(Tel, h, PlusEqDiv(Th));
        } else {
            LOOP(hel[i] = 0.0);
            hE(Tel, hel, EqDiv(Th));
            if (h != NULL)
                LOOP(h[i] += hel[i]);
        }
        
        // Formation enthalpy
        if (hf == NULL) {
            if (h != NULL)
                hF(h, PlusEqDiv(Th));
        } else {
            hF(hf, EqDiv(Th));
            if (h != NULL)
                LOOP(h[i] += hf[i]);
        }
    }
    
//    void hv(double Tv, double* const p_hv)
//    {
//        int ilevel = 0;
//        double sum, fac1;
//        LOOP_MOLECULES(
//            sum = 0.0;
//            for (int k = 0; k < mp_nvib[i]; ++k, ilevel++) {
//                fac1 = mp_vib_temps[ilevel] / Tv;
//                sum += fac1 / (std::exp(fac1) - 1.0);
//            }
//            p_hv[j] = sum;
//        )
//    }
//
//    void hel(double Tel, double* const p_hel)
//    {
//        HelFunctor()(Tel, p_hel, m_elec_data, Eq());
//    }

    /**
     * Computes the unitless species energy \f$ e_i/R_U T_h\f$ of each
     * species in thermal nonequilibrium, which is non-dimensionalized by the 
     * heavy particle translational temperature. 
     *
     * @param Th  - heavy particle translational temperature
     * @param Te  - free electron temperature
     * @param Tr  - mixture rotational temperature
     * @param Tv  - mixture vibrational temperature
     * @param Tel - mixture electronic temperature
     * @param e   - on return, the array of species non-dimensional energies
     * @param et  - if not NULL, the array of species translational energies
     * @param er  - if not NULL, the array of species rotational energies
     * @param ev  - if not NULL, the array of species vibrational energies
     * @param eel - if not NULL, the array of species electronic energies
     * @param ef  - if not NULL, the array of the species formation energies
     */
    void energy(
        double Th, double Te, double Tr, double Tv, double Tel, double* const e,
        double* const et, double* const er, double* const ev, double* const eel,
        double* const ef)
    {
        
        // Special case where we only want the total energy
        if (et == NULL && er == NULL && ev == NULL && eel == NULL && 
            ef == NULL && e != NULL) 
        {
            eT(Th, Te, e, Eq());
            eR(Tr, e, PlusEq());
            eV(Tv, e, PlusEq());
            eE(Tel, e, PlusEq());
            eF(e, PlusEq());
            LOOP(e[i] /= Th);
            return;
        }
        
        // Otherwise selectively choose what we want
        // Translational energy
        if (et == NULL) {
            if (e != NULL)
                eT(Th, Te, e, EqDiv(Th));
        } else {
            eT(Th, Te, et, EqDiv(Th));
            if (e != NULL)
                LOOP(e[i] = et[i]);
        }
        
        // Rotatonal energy
        if (er == NULL) {
            if (e != NULL)
                eR(Tr, e, PlusEqDiv(Th));
        } else {
            LOOP(er[i] = 0.0);
            eR(Tr, er, EqDiv(Th));
            if (e != NULL)
                LOOP_MOLECULES(e[j] += er[j]);
        }
        
        // Vibrational energy
        if (ev == NULL) {
            if (e != NULL)
                eV(Tv, e, PlusEqDiv(Th));
        } else {
            LOOP(ev[i] = 0.0);
            eV(Tv, ev, EqDiv(Th));
            if (e != NULL)
                LOOP_MOLECULES(e[j] += ev[j]);
        }
            
        // Electronic energy
        if (eel == NULL) {
            if (e != NULL)
                eE(Tel, e, PlusEqDiv(Th));
        } else {
            LOOP(eel[i] = 0.0);
            eE(Tel, eel, EqDiv(Th));
            if (e != NULL)
                LOOP(e[i] += eel[i]);
        }
        
        // Formation energy
        if (ef == NULL) {
            if (e != NULL)
                eF(e, PlusEqDiv(Th));
        } else {
            eF(ef, EqDiv(Th));
            if (e != NULL)
                LOOP(e[i] += ef[i]);
        }
    }

//    void ev(double Tv, double* const p_ev)
//    {
//        int ilevel = 0;
//        double sum, fac1;
//        LOOP_MOLECULES(
//            sum = 0.0;
//            for (int k = 0; k < mp_nvib[i]; ++k, ilevel++) {
//                fac1 = mp_vib_temps[ilevel] / Tv;
//                sum += fac1 / (std::exp(fac1) - 1.0);
//            }
//            p_ev[j] = sum;
//        )
//    }
//
//    void eel(double Tel, double* const p_eel)
//    {
//        EelFunctor()(Tel, p_eel, m_elec_data, Eq());
//    }

    /**
     * Computes the unitless species entropy \f$s_i/R_u\f$ allowing for thermal 
     * nonequilibrium.
     *
     * @param Th  - heavy particle translational temperature
     * @param Te  - free electron temperature
     * @param Tr  - mixture rotational temperature
     * @param Tv  - mixture vibrational temperature
     * @param Tel - mixture electronic temperature
     * @param s   - on return, the array of species non-dimensional entropies
     * @param st  - if not NULL, the array of species translational entropies
     * @param sr  - if not NULL, the array of species rotational entropies
     * @param sv  - if not NULL, the array of species vibrational entropies
     * @param sel - if not NULL, the array of species electronic entropies
     */
    void entropy(
        double Th, double Te, double Tr, double Tv, double Tel, double P,
        double* const s, double* const st, double* const sr, double* const sv,
        double* const sel)
    {
        // Special case where we only want the total entropy
        if (st == NULL && sr == NULL && sv == NULL && sel == NULL) {
            sT(Th, Te, P, s, Eq());
            sR(Tr, s, PlusEq());
            sV(Tv, s, PlusEq());
            sE(Tel, s, PlusEq());
            
            // Include spin contribution for free electron entropy
            if (m_has_electron)
                s[0] += std::log(2.0);
            
            return;
        }
        
        // Otherwise collect individual components
        // Translational entropy
        if (st == NULL)
            sT(Th, Te, P, s, Eq());
        else {
            sT(Th, Te, P, st, Eq());
            LOOP(s[i] = st[i]);
        }
        
        // Rotational entropy
        if (sr == NULL)
            sR(Tr, s, PlusEq());
        else {
            LOOP(sr[i] = 0.0);
            sR(Tr, sr, Eq());
            LOOP_MOLECULES(s[j] += sr[j]);
        }
        
        // Vibrational entropy
        if (sv == NULL)
            sV(Tv, s, PlusEq());
        else {
            LOOP(sv[i] = 0.0);
            sV(Tv, sv, Eq());
            LOOP_MOLECULES(s[j] += sv[j]);
        }
        
        // Electronic entropy
        if (sel == NULL)
            sE(Tel, s, PlusEq());
        else {
            LOOP(sel[i] = 0.0);
            sE(Tel, sel, Eq());
            LOOP(s[i] += sel[i]);
        }
        
        // Include spin contribution for free electron entropy
        if (m_has_electron)
            s[0] += std::log(2.0);
    }
    
    /**
     * Computes the unitless Gibbs free energy of each species i, 
     * \f$G_i / R_u T_h\f$ where \f$G_i\f$ is non-dimensionalized by the heavy
     * particle translational temperature.
     *
     * @todo Compute the individual components of the Gibbs function directly
     * instead of H - TS.
     */
    void gibbs(
        double Th, double Te, double Tr, double Tv, double Tel, double P,
        double* const g, double* const gt, double* const gr, double* const gv,
        double* const gel)
    {        
        // First compute the non-dimensional enthalpy
        enthalpy(Th, Te, Tr, Tv, Tel, g, NULL, NULL, NULL, NULL, NULL);

        // Subtract the entropies
        sT(Th, Te, P, g, MinusEq());
        sR(Tr, g, MinusEq());
        sV(Tv, g, MinusEq());
        sE(Tel, g, MinusEq());

        // Account for spin of free electrons
        if (m_has_electron)
            g[0] -= std::log(2.0);
    }

private:

    typedef Equals<double> Eq;
    typedef EqualsYDivAlpha<double> EqDiv;    
    typedef PlusEqualsYDivAlpha<double> PlusEqDiv;
    typedef PlusEquals<double> PlusEq;
    typedef MinusEquals<double> MinusEq;

    class ElecBFacsFunctor
    {
    public:
        typedef ElectronicData DataProvider;

        void operator() (double T, double* p_f, const DataProvider& data) const
        {
            (*this)(T, p_f, data, Equals<double>());
        }

        template <typename OP>
        void operator () (
            double T, double* p_f, const DataProvider& data, const OP& op) const
        {
            double fac;
            int ilevel = 0;

            for (unsigned int i = 0; i < data.nheavy; ++i) {
                p_f[3*i+0] = 0.0;
                p_f[3*i+1] = 0.0;
                p_f[3*i+2] = 0.0;

                for (int k = 0; k < data.p_nelec[i]; ++k, ilevel++) {
                    fac = data.p_levels[ilevel].g *
                        std::exp(-data.p_levels[ilevel].theta / T);
                    p_f[3*i+0] += fac;
                    p_f[3*i+1] += fac * data.p_levels[ilevel].theta;
                    p_f[3*i+2] += fac * data.p_levels[ilevel].theta *
                        data.p_levels[ilevel].theta;
                }
            }
        }
    };
    
protected:

    /**
     * Loads all of the species from the RRaHO database.
     */
    virtual void loadAvailableSpecies(std::list<Species>& species)
    {
        IO::XmlDocument species_doc(databaseFileName("species.xml", "thermo"));
        IO::XmlElement::const_iterator species_iter = species_doc.root().begin();
        
        for ( ; species_iter != species_doc.root().end(); ++species_iter) {
            // Add the species to the list
            species.push_back(*species_iter);
            
            // We can also add all of the excited states as implicitly defined
            // species
            IO::XmlElement::const_iterator rraho_iter =
                species_iter->findTagWithAttribute(
                    "thermodynamics", "type", "RRaHO");
            
            if (rraho_iter == species_iter->end())
                continue;
            
            Species& ground_state = species.back();
            ParticleRRaHO rraho(*rraho_iter);

	    // Add electronic levels as pseudo-species
            //for (size_t i = 0; i < rraho.nElectronicLevels(); ++i)
            //    species.push_back(Species(ground_state, i));

	    // Add vibrational levels as pseudo-species
	    // Maybe just consider ground electronic state ...
	    for (size_t e = 0; e < rraho.nElectronicLevels(); ++e) {
                for (size_t v = 0; v < rraho.nVibrationalLevels(e); ++v) {
                    species.push_back(Species(ground_state, e, v));
                }
            }
        }
    }
    
    /**
     * Load thermodynamic data from the species list.
     */
    virtual void loadThermodynamicData()
    {
        m_ns = species().size();
        m_has_electron = (species()[0].type() == ELECTRON);
        
        // Load the RRaHO models for each of the needed species
        IO::XmlDocument species_doc(databaseFileName("species.xml", "thermo"));
        
        vector<ParticleRRaHO> rrahos;
        map<std::string, const ParticleRRaHO*> to_expand;
        
        for (int i = 0; i < m_ns; ++i) {
            if (species()[i].name() == species()[i].groundStateName()) {
                rrahos.push_back(*(species_doc.root().findTagWithAttribute(
                    "species", "name", species()[i].groundStateName())->
                        findTagWithAttribute("thermodynamics", "type", "RRaHO")));
            }
            else {
                const ParticleRRaHO* p_rraho = to_expand[species()[i].groundStateName()];
                if (p_rraho == NULL) {
                    p_rraho = new ParticleRRaHO(
                        *(species_doc.root().findTagWithAttribute(
                            "species", "name", species()[i].groundStateName())->
                                findTagWithAttribute("thermodynamics", "type", "RRaHO")));
                    to_expand[species()[i].groundStateName()] = p_rraho;
                }
                
                rrahos.push_back(ParticleRRaHO(*p_rraho, species()[i].level()));
            }
        }
        
        map<std::string, const ParticleRRaHO*>::iterator iter =
            to_expand.begin();
        while (iter != to_expand.end()) {
            delete iter->second;
            iter++;
        }

        // Determine the number and indices of the atoms and molecules
        vector<int> atom_indices;
        vector<int> molecule_indices;
        
        LOOP(
            switch(species()[i].type()) {
                case ATOM:
                    atom_indices.push_back(i);
                    break;
                case MOLECULE:
                    molecule_indices.push_back(i);
                    break;
                default:
                    break;
            }
        )
        
        m_na = atom_indices.size();
        m_nm = molecule_indices.size();
        
        // Order the atoms first followed by the molecules
        mp_indices = new int [m_na + m_nm];
        copy(atom_indices.begin(), atom_indices.end(), mp_indices);
        copy(molecule_indices.begin(), molecule_indices.end(), mp_indices+m_na);
        
        // Store the species constants found in the translational entropy term
        mp_lnqtmw = new double [m_ns];
        double qt = 2.5*std::log(KB) + 1.5*std::log(TWOPI/(NA*HP*HP));
        LOOP(mp_lnqtmw[i] = qt + 1.5*std::log(species()[i].molecularWeight()))
        
        // Store the species formation enthalpies in K
        mp_hform = new double [m_ns];
        LOOP(mp_hform[i] = rrahos[i].formationEnthalpy() / RU)

        // Store the species formation enrergies in K
        mp_eform = new double [m_ns];
        LOOP(mp_eform[i] = rrahos[i].formationEnergy() / RU)

        // Store the molecule's rotational energy parameters
        mp_rot_data = new RotData [m_nm];
        LOOP_MOLECULES(
            const ParticleRRaHO& rraho = rrahos[j];
            int linear = rraho.linearity();
            mp_rot_data[i].linearity  = linear / 2.0;
            mp_rot_data[i].ln_omega_t = 
                std::log(rraho.rotationalTemperature()) + 2.0 / linear *
                std::log(rraho.stericFactor());
        )
        
        // Store the vibrational temperatures of all the molecules in a compact
        // form
        mp_nvib = new int [m_nm];
        int nvib = 0;
        LOOP_MOLECULES(
            mp_nvib[i] = rrahos[j].nVibrationalLevels();
            nvib += mp_nvib[i];
        )
        
        mp_vib_temps = new double [nvib];
        int itemp = 0;
        LOOP_MOLECULES(
            const ParticleRRaHO& rraho = rrahos[j];
            for (int k = 0; k < mp_nvib[i]; ++k, itemp++)
                mp_vib_temps[itemp] = rraho.vibrationalEnergy(k);
        )
        
        // Finally store the electronic energy levels in a compact form like the
        // vibrational energy levels
        m_elec_data.p_nelec = new int [m_na + m_nm];
        m_elec_data.nlevels = 0;
        LOOP_HEAVY(
            m_elec_data.p_nelec[i] = rrahos[j].nElectronicLevels();
            m_elec_data.nlevels += m_elec_data.p_nelec[i];
        )
        
        m_elec_data.p_levels = new ElecLevel [m_elec_data.nlevels];
        int ilevel = 0;
        LOOP_HEAVY(
            const ParticleRRaHO& rraho = rrahos[j];
            for (int k = 0; k < m_elec_data.p_nelec[i]; ++k, ilevel++) {
                m_elec_data.p_levels[ilevel].g =
                    rraho.electronicEnergy(k).first;
                m_elec_data.p_levels[ilevel].theta = 
                    rraho.electronicEnergy(k).second;
            }
        )        
        
        m_elec_data.offset = (m_has_electron ? 1 : 0);
        m_elec_data.nheavy = m_na + m_nm;
        
        if (m_use_tables) {
            mp_el_bfac_table = new Mutation::Utilities::LookupTable
                <double, double, ElecBFacsFunctor>(
                50.0, 50000.0, 3*(m_na+m_nm), m_elec_data, 0.005);
        }
        
        mp_el_bfacs = new double [3*(m_na+m_nm)];

        // Compute the contribution of the partition functions at the standard
        // state temperature to the species enthalpies
        mp_part_sst = new double [m_ns];
        double Tss = standardTemperature();
        hT(Tss, Tss, mp_part_sst, Eq());
        hR(Tss, mp_part_sst, PlusEq());
        hV(Tss, mp_part_sst, PlusEq());
        hE(Tss, mp_part_sst, PlusEq());

        eT(Tss, Tss, mp_part_sst, Eq());
        eR(Tss, mp_part_sst, PlusEq());
        eV(Tss, mp_part_sst, PlusEq());
        eE(Tss, mp_part_sst, PlusEq());
    }

private:

    void updateElecBoltzmannFactors(double T)
    {
        if (std::abs(1.0 - m_last_bfacs_T / T) < 1.0e-16)
            return;

        if (m_use_tables)
            mp_el_bfac_table->lookup(T, mp_el_bfacs);
        else
            ElecBFacsFunctor()(T, mp_el_bfacs, m_elec_data);

        m_last_bfacs_T = T;
    }
    
    /**
     * Computes the translational Cp/Ru for each species.
     */
    template <typename OP>
    void cpT(double* const cp, const OP& op) {
        LOOP(op(cp[i], 2.5));
    }

    /**
     * Computes the rotational Cp/Ru for each species.
     */
    template <typename OP>
    void cpR(double* const cp, const OP& op) {
        op(cp[0], 0.0);
        LOOP_ATOMS(op(cp[j], 0.0));
        LOOP_MOLECULES(op(cp[j], mp_rot_data[i].linearity));
    }

    /**
     * Computes the vibratinoal Cp/Ru for each species.
     */
    template <typename OP>
    void cpV(double Tv, double* const cp, const OP& op) {
        int ilevel = 0;
        double sum, fac1, fac2;
        op(cp[0], 0.0);
        LOOP_ATOMS(op(cp[j], 0.0));
        LOOP_MOLECULES(
            sum = 0.0;
            for (int k = 0; k < mp_nvib[i]; ++k, ilevel++) {
                fac1 = mp_vib_temps[ilevel] / Tv;
                fac2 = std::exp(fac1);
                fac1 *= fac1*fac2;
                fac2 -= 1.0;
                sum += fac1/(fac2*fac2);
            }
            op(cp[j], sum);
        )
    }

    /**
     * Computes the electronic specific heat of each species and applies the
     * value to the array using the given operation.
     */
    template <typename OP>
    void cpE(double T, double* const p_cp, const OP& op)
    {
        updateElecBoltzmannFactors(T);
        op(p_cp[0], 0.0);

        double* facs = mp_el_bfacs;
        for (unsigned int i = 0; i < m_elec_data.nheavy; ++i, facs += 3) {
            if (m_elec_data.p_nelec[i] > 1)
                op(p_cp[i+m_elec_data.offset],
                    (facs[2]*facs[0]-facs[1]*facs[1])/(T*T*facs[0]*facs[0]));
            else
                op(p_cp[i+m_elec_data.offset], 0.0);
        }
    }

    /** TODO
     * Computes the translational Cv/Ru for each species.
     */
    template <typename OP>
    void cvT(double* const cv, const OP& op) 
    {
    //    LOOP(op(cv[i], 2.5));
    }

    /** TODO
     * Computes the rotational Cv/Ru for each species.
     */
    template <typename OP>
    void cvR(double* const cv, const OP& op) 
    {
    //    op(cv[0], 0.0);
    //    LOOP_ATOMS(op(cv[j], 0.0));
    //    LOOP_MOLECULES(op(cv[j], mp_rot_data[i].linearity));
    }

    /** TODO
     * Computes the vibratinoal Cv/Ru for each species.
     */
    template <typename OP>
    void cvV(double Tv, double* const cv, const OP& op) 
    {
    //    int ilevel = 0;
    //    double sum, fac1, fac2;
    //    op(cv[0], 0.0);
    //    LOOP_ATOMS(op(cv[j], 0.0));
    //    LOOP_MOLECULES(
    //        sum = 0.0;
    //        for (int k = 0; k < mp_nvib[i]; ++k, ilevel++) {
    //            fac1 = mp_vib_temps[ilevel] / Tv;
    //            fac2 = std::exp(fac1);
    //            fac1 *= fac1*fac2;
    //            fac2 -= 1.0;
    //            sum += fac1/(fac2*fac2);
    //        }
    //        op(cv[j], sum);
    //    )
    }

    /** TODO
     * Computes the electronic specific heat of each species and applies the
     * value to the array using the given operation.
     */
    template <typename OP>
    void cvE(double T, double* const p_cv, const OP& op)
    {
    //    updateElecBoltzmannFactors(T);
    //    op(p_cv[0], 0.0);

    //    double* facs = mp_el_bfacs;
    //    for (unsigned int i = 0; i < m_elec_data.nheavy; ++i, facs += 3) {
    //        if (m_elec_data.p_nelec[i] > 1)
    //            op(p_cv[i+m_elec_data.offset],
    //                (facs[2]*facs[0]-facs[1]*facs[1])/(T*T*facs[0]*facs[0]));
    //        else
    //            op(p_cv[i+m_elec_data.offset], 0.0);
    //    }
    }

    /**
     * Computes the translational enthalpy of each species in K.
     */
    template <typename OP>
    void hT(double T, double Te, double* const h, const OP& op) {
        if (m_has_electron)
            op(h[0], 2.5 * Te);
        LOOP_HEAVY(op(h[j], 2.5 * T))
    }
    
    /**
     * Computes the rotational enthalpy of each species in K.
     */
    template <typename OP>
    void hR(double T, double* const h, const OP& op) {
        LOOP_MOLECULES(op(h[j], mp_rot_data[i].linearity * T))
    }
    
    /**
     * Computes the vibrational enthalpy of each species in K.
     */
    template <typename OP>
    void hV(double T, double* const h, const OP& op) {
        if (T < 10.0) {
            LOOP_MOLECULES(op(h[j], 0.0));
        } else {
            int ilevel = 0;
            double sum;
            LOOP_MOLECULES(
                sum = 0.0;
                for (int k = 0; k < mp_nvib[i]; ++k, ilevel++)
                    sum += mp_vib_temps[ilevel] /
                        (std::exp(mp_vib_temps[ilevel] / T) - 1.0);
                op(h[j], sum);
            )
        }
    }
    
    /**
     * Computes the electronic enthalpy of each species in K and applies the
     * value to the enthalpy array using the given operation.
     */
    template <typename OP>
    void hE(double T, double* const p_h, const OP& op)
    {
        updateElecBoltzmannFactors(T);
        op(p_h[0], 0.0);

        double* facs = mp_el_bfacs;
        for (int i = 0; i < m_elec_data.nheavy; ++i, facs += 3) {
            if (facs[0] > 0)
                op(p_h[i+m_elec_data.offset], facs[1]/facs[0]);
            else
                op(p_h[i+m_elec_data.offset], 0.0);
        }
    }
    
    /**
     * Computes the formation enthalpy of each species in K.
     */
    template <typename OP>
    void hF(double* const h, const OP& op) {
        LOOP(op(h[i], mp_hform[i] - mp_part_sst[i]))
    }
    
    /** TODO
     * Computes the translational energy of each species in K.
     */
    template <typename OP>
    void eT(double T, double Te, double* const e, const OP& op) {
    //    if (m_has_electron)
    //        op(h[0], 2.5 * Te);
    //    LOOP_HEAVY(op(h[j], 2.5 * T))
    }
    
    /** TODO
     * Computes the rotational energy of each species in K.
     */
    template <typename OP>
    void eR(double T, double* const e, const OP& op) {
    //    LOOP_MOLECULES(op(h[j], mp_rot_data[i].linearity * T))
    }
    
    /** TODO
     * Computes the vibrational energy of each species in K.
     */
    template <typename OP>
    void eV(double T, double* const e, const OP& op) {
    //    if (T < 10.0) {
    //        LOOP_MOLECULES(op(h[j], 0.0));
    //    } else {
    //        int ilevel = 0;
    //        double sum;
    //        LOOP_MOLECULES(
    //            sum = 0.0;
    //            for (int k = 0; k < mp_nvib[i]; ++k, ilevel++)
    //                sum += mp_vib_temps[ilevel] /
    //                    (std::exp(mp_vib_temps[ilevel] / T) - 1.0);
    //            op(h[j], sum);
    //        )
    //    }
    }
    
    /** TODO
     * Computes the electronic energy of each species in K and applies the
     * value to the energy array using the given operation.
     */
    template <typename OP>
    void eE(double T, double* const p_e, const OP& op)
    {
    //    updateElecBoltzmannFactors(T);
    //    op(p_h[0], 0.0);

    //    double* facs = mp_el_bfacs;
    //    for (int i = 0; i < m_elec_data.nheavy; ++i, facs += 3) {
    //        if (facs[0] > 0)
    //            op(p_h[i+m_elec_data.offset], facs[1]/facs[0]);
    //        else
    //            op(p_h[i+m_elec_data.offset], 0.0);
    //    }
    }
    
    /** TODO
     * Computes the formation energy of each species in K.
     */
    template <typename OP>
    void eF(double* const e, const OP& op) {
        LOOP(op(e[i], mp_eform[i] - mp_part_sst[i]))
    }

    /**
     * Computes the unitless translational entropy of each species.
     */
    template <typename OP>
    void sT(double Th, double Te, double P, double* const s, const OP& op) {
        double fac = 2.5 * (1.0 + std::log(Th)) - std::log(P);
        if (m_has_electron)
            op(s[0], 2.5 * std::log(Te / Th) + fac + mp_lnqtmw[0]);        
        for (int i = (m_has_electron ? 1 : 0); i < m_ns; ++i)
            op(s[i], fac + mp_lnqtmw[i]);
    }
    
    /**
     * Computes the unitless rotational entropy of each species.
     */
    template <typename OP>
    void sR(double T, double* const s, const OP& op) {
        const double onelnT = 1.0 + std::log(T);
        LOOP_MOLECULES(
            op(s[j], mp_rot_data[i].linearity * (onelnT - 
                mp_rot_data[i].ln_omega_t));
        )
    }
    
    /**
     * Computes the unitless vibrational entropy of each species.
     */
    template <typename OP>
    void sV(double T, double* const s, const OP& op) {
        int ilevel = 0;
        double fac, sum1, sum2;
        LOOP_MOLECULES(
            sum1 = sum2 = 0.0;
            for (int k = 0; k < mp_nvib[i]; ++k, ilevel++) {
                fac  =  std::exp(mp_vib_temps[ilevel] / T);
                sum1 += mp_vib_temps[ilevel] / (fac - 1.0);
                sum2 += std::log(1.0 - 1.0 / fac);
            }
            op(s[j], (sum1 / T - sum2));
        )
    }
    
    /**
     * Computes the unitless electronic entropy of each species.
     */
    template <typename OP>
    void sE(double T, double* const p_s, const OP& op) {
        updateElecBoltzmannFactors(T);
        op(p_s[0], 0.0);

        double* facs = mp_el_bfacs;
        for (int i = 0; i < m_elec_data.nheavy; ++i, facs += 3) {
            if (facs[0] > 0)
                op(p_s[i+m_elec_data.offset],
                    (facs[1]/(facs[0]*T) + std::log(facs[0])));
            else
                op(p_s[i+m_elec_data.offset], 0.0);
        }
    }

private:
    
    int m_ns;
    int m_na;
    int m_nm;
    
    bool m_has_electron;
    bool m_use_tables;
    
    double* mp_lnqtmw;
    double* mp_hform;
    double* mp_eform;
    double* mp_part_sst;
    
    int*       mp_indices;
    RotData*   mp_rot_data;
    
    int*       mp_nvib;
    double*    mp_vib_temps;
    
    ElectronicData m_elec_data;
    Mutation::Utilities::LookupTable<double, double, ElecBFacsFunctor>* mp_el_bfac_table;
    double* mp_el_bfacs;
    double m_last_bfacs_T;

    //Mutation::Utilities::LookupTable<double, double, HelFunctor>* mp_hel_table;
    //Mutation::Utilities::LookupTable<double, double, EelFunctor>* mp_eel_table;
    //Mutation::Utilities::LookupTable<double, double, SelFunctor>* mp_sel_table;
    //Mutation::Utilities::LookupTable<double, double, CpelFunctor>* mp_cpel_table;
    //Mutation::Utilities::LookupTable<double, double, CvelFunctor>* mp_cvel_table;

}; // class RrahoDB

#undef LOOP
#undef LOOP_HEAVY
#undef LOOP_MOLECULES

// Register the RRaHO model with the other thermodynamic databases
Utilities::Config::ObjectProvider<RrahoDB, ThermoDB> rrahoDB("RRaHO");

    } // namespace Thermodynamics
} // namespace Mutation
