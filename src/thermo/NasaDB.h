
//#include <iostream>
//#include <fstream>
//#include <iomanip>
//#include <string>
//#include <map>
//#include <vector>
//
//#include "AutoRegistration.h"
//#include "Nasa7Polynomial.h"
//#include "Nasa9Polynomial.h"
//#include "Utilities.h"

#include <list>
#include <fstream>
#include <vector>

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
        : ThermoDB(298.15, 1.0E5), m_ns(0)
    { }
    
    virtual ~NasaDB() {};
    
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
    
    void loadAvailableSpecies(
        std::list<Species>& species_list, const std::vector<Element>& elements);
    
    void loadThermodynamicData();
    
    virtual std::string filename() const = 0;
    
    virtual void skipHeader(std::ifstream& is) const = 0;
    
    virtual Species loadSpecies(
        std::ifstream& is, const std::vector<Element>& elements) const = 0;
    
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
    PolynomialType::computeParams(Th, mp_params, PolynomialType::CP);
    for (size_t i = 0; i < m_ns; ++i)
        m_polynomials[i].cp(mp_params, cp[i]);
    
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
    PolynomialType::computeParams(Th, mp_params, PolynomialType::ENTHALPY);
    for (size_t i = 0; i < m_ns; ++i)
        m_polynomials[i].enthalpy(mp_params, h[i]);
    
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
    PolynomialType::computeParams(Th, mp_params, PolynomialType::ENTROPY);
    for (size_t i = 0; i < m_ns; ++i)
        m_polynomials[i].entropy(mp_params, s[i]);
    
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
    PolynomialType::computeParams(Th, mp_params, PolynomialType::GIBBS);
    for (size_t i = 0; i < m_ns; ++i)
        m_polynomials[i].gibbs(mp_params, g[i]);
    
    if (gt != NULL) std::fill(gt, gt+m_ns, 0.0);
    if (gr != NULL) std::fill(gr, gr+m_ns, 0.0);
    if (gv != NULL) std::fill(gv, gv+m_ns, 0.0);
    if (gel != NULL) std::fill(gel, gel+m_ns, 0.0);
}

template <typename PolynomialType>
void NasaDB<PolynomialType>::loadAvailableSpecies(
    std::list<Species>& species_list, const std::vector<Element>& elements)
{
    // Open the database file
    std::string db_path =
        Utilities::getEnvironmentVariable("MPP_DATA_DIRECTORY") + "/thermo/" +
        filename();
    std::ifstream file(db_path.c_str());
    
    if (!file.is_open()) {
        std::cerr << "Could not open " << db_path << "!\n" << std::endl;
        exit(1);
    }
    
    // Move to the beginning of the first species
    skipHeader(file);
    
    // Load all of the available species until we reach the end of the file
    species_list.push_back(loadSpecies(file, elements));
    while (species_list.back().name() != "")
        species_list.push_back(loadSpecies(file, elements));
    species_list.pop_back();
    
    // Close the database file
    file.close();
}

template <typename PolynomialType>
void NasaDB<PolynomialType>::loadThermodynamicData()
{
    // Open the database file
    std::string db_path =
        Utilities::getEnvironmentVariable("MPP_DATA_DIRECTORY") + "/thermo/" +
        filename();
    std::ifstream file(db_path.c_str());
    
    if (!file.is_open()) {
        std::cerr << "Could not open " << db_path << "!\n" << std::endl;
        exit(1);
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


    