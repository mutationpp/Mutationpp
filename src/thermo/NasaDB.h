/**
 * @file NasaDB.h
 *
 * @brief Defines base class for NASA thermodynamic databases.
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

#include <list>
#include <fstream>
#include <vector>
#include <cassert>

#include "ThermoDB.h"
#include "Species.h"
#include "Utilities.h"

namespace Mutation {
    namespace Thermodynamics {

/**
 * Template for a thermodynamic database which uses a set of NASA polynomials.
 * This class abstracts away the implementation of the polynomial type from the
 * database container.
 */
template <typename PolynomialType>
class NasaDB : public ThermoDB
{
public:
    NasaDB()
        : ThermoDB(298.15, 101325.0), m_ns(0)
    { }
    
    virtual ~NasaDB() {};
    
    virtual bool speciesThermoValidAtT(const size_t i, const double T) const {
        assert(i < m_ns);
        return (T > m_polynomials[i].minT() && T <= m_polynomials[i].maxT());
    }
    
    void cp(
        double Th, double Te, double Tr, double Tv, double Tel, 
        double* const cp, double* const cpt, double* const cpr, 
        double* const cpv, double* const cpel);
    
    void enthalpy(
        double Th, double Te, double Tr, double Tv, double Tel, 
        double* const h, double* const ht, double* const hr, double* const hv, 
        double* const hel, double* const hf);
    
    void entropy(
        double Th, double Te, double Tr, double Tv, double Tel, double P,
        double* const s, double* const st, double* const sr, double* const sv, 
        double* const sel);
    
    void gibbs(
        double Th, double Te, double Tr, double Tv, double Tel, double P,
        double* const g, double* const gt, double* const gr, double* const gv, 
        double* const gel);

protected:
    
    void loadAvailableSpecies(std::list<Species>& species_list);
    
    void loadThermodynamicData();
    
    virtual std::string filename() const = 0;
    
    virtual void skipHeader(std::ifstream& is) const = 0;
    
    virtual Species loadSpecies(std::ifstream& is) const = 0;
    
    virtual void loadPolynomials(
        std::ifstream& is, std::vector<PolynomialType>& polynomials) = 0;

private:
    
    size_t m_ns;
    std::vector<PolynomialType> m_polynomials;
    double mp_params[8];
};

template <typename PolynomialType>
void NasaDB<PolynomialType>::cp(
    double Th, double Te, double Tr, double Tv, double Tel,
    double* const cp, double* const cpt, double* const cpr, 
    double* const cpv, double* const cpel)
{
    
    if (cp != NULL) {
        PolynomialType::computeParams(Th, mp_params, PolynomialType::CP);
        for (size_t i = 0; i < m_ns; ++i)
            m_polynomials[i].cp(mp_params, cp[i]);
    }
    
    if (cpt != NULL) std::fill(cpt, cpt+m_ns, 0.0);
    if (cpr != NULL) std::fill(cpr, cpr+m_ns, 0.0);
    if (cpv != NULL) std::fill(cpv, cpv+m_ns, 0.0);
    if (cpel != NULL) std::fill(cpel, cpel+m_ns, 0.0);
}
    
template <typename PolynomialType>
void NasaDB<PolynomialType>::enthalpy(
    double Th, double Te, double Tr, double Tv, double Tel,
    double* const h, double* const ht, double* const hr, double* const hv, 
    double* const hel, double* const hf)
{
    if (h != NULL) {
        PolynomialType::computeParams(Th, mp_params, PolynomialType::ENTHALPY);
        for (size_t i = 0; i < m_ns; ++i)
            m_polynomials[i].enthalpy(mp_params, h[i]);
    }
    
    if (ht != NULL) std::fill(ht, ht+m_ns, 0.0);
    if (hr != NULL) std::fill(hr, hr+m_ns, 0.0);
    if (hv != NULL) std::fill(hv, hv+m_ns, 0.0);
    if (hel != NULL) std::fill(hel, hel+m_ns, 0.0);
    if (hf != NULL) std::fill(hf, hf+m_ns, 0.0);      
}

template <typename PolynomialType>
void NasaDB<PolynomialType>::entropy(
    double Th, double Te, double Tr, double Tv, double Tel, double P,
    double* const s, double* const st, double* const sr, double* const sv, 
    double* const sel)
{
    if (s != NULL) {
        PolynomialType::computeParams(Th, mp_params, PolynomialType::ENTROPY);
        for (size_t i = 0; i < m_ns; ++i)
            m_polynomials[i].entropy(mp_params, s[i]);
    }
    
    if (st != NULL) std::fill(st, st+m_ns, 0.0);
    if (sr != NULL) std::fill(sr, sr+m_ns, 0.0);
    if (sv != NULL) std::fill(sv, sv+m_ns, 0.0);
    if (sel != NULL) std::fill(sel, sel+m_ns, 0.0);
}

template <typename PolynomialType>
void NasaDB<PolynomialType>::gibbs(
    double Th, double Te, double Tr, double Tv, double Tel, double P,
    double* const g, double* const gt, double* const gr, double* const gv, 
    double* const gel)
{
    if (g != NULL) {
        PolynomialType::computeParams(Th, mp_params, PolynomialType::GIBBS);
        for (size_t i = 0; i < m_ns; ++i)
            m_polynomials[i].gibbs(mp_params, g[i]);
    }
    
    if (gt != NULL) std::fill(gt, gt+m_ns, 0.0);
    if (gr != NULL) std::fill(gr, gr+m_ns, 0.0);
    if (gv != NULL) std::fill(gv, gv+m_ns, 0.0);
    if (gel != NULL) std::fill(gel, gel+m_ns, 0.0);
}

template <typename PolynomialType>
void NasaDB<PolynomialType>::loadAvailableSpecies(
    std::list<Species>& species_list)
{
    // Open the database file
    std::string db_path =
        Utilities::databaseFileName(filename(), "thermo", ".dat");
    std::ifstream file(db_path.c_str());

    if (!file.is_open()) {
        throw FileNotFoundError(db_path)
            << "Could not find thermodynamic database.";
    }
    
    // Move to the beginning of the first species
    skipHeader(file);
    
    // Load all of the available species until we reach the end of the file
    species_list.push_back(loadSpecies(file));
    while (species_list.back().name() != "")
        species_list.push_back(loadSpecies(file));
    species_list.pop_back();
    
    // Close the database file
    file.close();
}

template <typename PolynomialType>
void NasaDB<PolynomialType>::loadThermodynamicData()
{
    // Open the database file
    std::string db_path =
        Utilities::databaseFileName(filename(), "thermo", ".dat");
    std::ifstream file(db_path.c_str());
    
    if (!file.is_open()) {
        throw FileNotFoundError(db_path)
            << "Could not find thermodynamic database.";
    }
    
    // Move to the beginning of the first species
    skipHeader(file);
    
    // Ask the concrete class to fill the polynomials
    m_ns = species().size();
    m_polynomials.resize(m_ns);
    loadPolynomials(file, m_polynomials);
    
    // Close the database file
    file.close();
}

    } // namespace Thermodynamics
} // namespace Mutation


    
