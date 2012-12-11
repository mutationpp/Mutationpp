
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <map>
#include <vector>

#include "AutoRegistration.h"
#include "Nasa7Polynomial.h"
#include "Nasa9Polynomial.h"
#include "ThermoDB.h"
#include "Species.h"
#include "utilities.h"

/**
 * Template for a thermodynamic database which uses a set of NASA polynomials.
 * This class abstracts away the implementation of the polynomial type from the
 * database container.
 */
template<typename POLYNOMIAL, char* DBNAME>
class NasaDB : public ThermoDB
{
public:

    using ThermoDB::m_ns;

    NasaDB(const std::vector<Species>& species)
        : ThermoDB(species)
    {
         // Open the database file
        std::string thermo_directory = 
            utils::getEnvironmentVariable("MPP_DATA_DIRECTORY") + "/thermo";
        std::string db_path = thermo_directory + "/" + std::string(DBNAME);
        std::ifstream file(db_path.c_str());
        
        if (!file.is_open()) {
            std::cerr << "Could not open " << db_path << "!\n" << std::endl;
            exit(1);
        }
        
        // Make a map of the species names corresponding to the species list
        // that are linked to their position in the list
        std::map<std::string,int> species_names;
        std::map<std::string,int>::iterator iter;
        for (int i = 0; i < species.size(); ++i)
            species_names[species[i].nasaName(POLYNOMIAL::nCoefficients())] = i;
            
        std::string line, name;
        
        // Read everything in the file or load all species
        while (std::getline(file, line)) {
            // Get the first "word" in the line which could be a species name
            name = line.substr(0, line.find(" "));
            
            // "END" signifies the end of the thermodynamic data
            if (name == "END")
                break;
            
            if ((iter = species_names.find(name)) != species_names.end()) {
                // Store the polynomial in the vector and the species position
                // in the species list
                m_polynomials.push_back(std::pair<int, POLYNOMIAL>(
                        iter->second, POLYNOMIAL()));
                
                // Load the polynomial
                file.seekg(
                    -static_cast<int>(line.length()+1), std::ios_base::cur) >> 
                    m_polynomials.back().second;
                
                // Remove the species from the list of species left to load
                species_names.erase(iter);
                
                // If the list is now empty then we don't need to look anymore
                if (species_names.empty())
                    break;
            }
        }
        
        // Be sure to close the file afterward
        file.close();
        
        // Check and see if we got all the species we needed
        if (species_names.empty())
            return;
        
        // If we get here then not all required species were loaded
        std::cout << "Could not find all required species in the NASA polynomial "
                  << "database (" << DBNAME << ")!" << std::endl;
        std::cout << "Missing species:" << std::endl;
        
        std::map<std::string, int>::const_iterator missing_iter = 
            species_names.begin();
        
        for ( ; missing_iter != species_names.end(); ++missing_iter)
            std::cout << std::setw(10) << missing_iter->first << std::endl;
        
        std::exit(1);
    }
    
    virtual ~NasaDB() {};
    
    /**
     * Returns the standard state temperature used for the reference state in
     * this database (298.15 K).
     */
    double standardTemperature() const {
        return 298.15;
    }
    
    /**
     * Returns the standard state pressure used for the reference state in this
     * databse (1 bar = 100,000 Pa).
     */
    double standardPressure() const {
        return 100000.0;
    }
    
    void cp(
        double Th, double Te, double Tr, double Tv, double Tel, 
        double* const cp, double* const cpt, double* const cpr, 
        double* const cpv, double* const cpel)
    {
        POLYNOMIAL::computeParams(Th, mp_params, POLYNOMIAL::CP);
        typename PolynomialList::const_iterator iter = m_polynomials.begin();
        for (; iter != m_polynomials.end(); ++iter)
            iter->second.cp(mp_params, cp[iter->first]);
            
        if (cpt != NULL) std::fill(cpt, cpt+m_ns, 0.0);
        if (cpr != NULL) std::fill(cpr, cpr+m_ns, 0.0);
        if (cpv != NULL) std::fill(cpv, cpv+m_ns, 0.0);
        if (cpel != NULL) std::fill(cpel, cpel+m_ns, 0.0);
    }
    
    void enthalpy(
        double Th, double Te, double Tr, double Tv, double Tel, 
        double* const h, double* const ht, double* const hr, double* const hv, 
        double* const hel, double* const hf)
    {
        POLYNOMIAL::computeParams(Th, mp_params, POLYNOMIAL::ENTHALPY);
        typename PolynomialList::const_iterator iter = m_polynomials.begin();
        for (int i = 0; iter != m_polynomials.end(); ++iter, ++i)
            iter->second.enthalpy(mp_params, h[iter->first]);
        
        if (ht != NULL) std::fill(ht, ht+m_ns, 0.0);
        if (hr != NULL) std::fill(hr, hr+m_ns, 0.0);
        if (hv != NULL) std::fill(hv, hv+m_ns, 0.0);
        if (hel != NULL) std::fill(hel, hel+m_ns, 0.0);
        if (hf != NULL) std::fill(hf, hf+m_ns, 0.0);      
    }

    void entropy(
        double Th, double Te, double Tr, double Tv, double Tel, double P,
        double* const s, double* const st, double* const sr, double* const sv, 
        double* const sel)
    {
        POLYNOMIAL::computeParams(Th, mp_params, POLYNOMIAL::ENTROPY);
        typename PolynomialList::const_iterator iter = m_polynomials.begin();
        for (int i = 0; iter != m_polynomials.end(); ++iter, ++i)
            iter->second.entropy(mp_params, s[iter->first]);
        
        if (st != NULL) std::fill(st, st+m_ns, 0.0);
        if (sr != NULL) std::fill(sr, sr+m_ns, 0.0);
        if (sv != NULL) std::fill(sv, sv+m_ns, 0.0);
        if (sel != NULL) std::fill(sel, sel+m_ns, 0.0);
    }

    void gibbs(
        double Th, double Te, double Tr, double Tv, double Tel, double P,
        double* const g, double* const gt, double* const gr, double* const gv, 
        double* const gel)
    {
        POLYNOMIAL::computeParams(Th, mp_params, POLYNOMIAL::GIBBS);
        typename PolynomialList::const_iterator iter = m_polynomials.begin();
        for (int i = 0; iter != m_polynomials.end(); ++iter, ++i)
            iter->second.gibbs(mp_params, g[iter->first]);
        
        if (gt != NULL) std::fill(gt, gt+m_ns, 0.0);
        if (gr != NULL) std::fill(gr, gr+m_ns, 0.0);
        if (gv != NULL) std::fill(gv, gv+m_ns, 0.0);
        if (gel != NULL) std::fill(gel, gel+m_ns, 0.0);
    }

private:
   
    typedef std::vector<std::pair<int, POLYNOMIAL> > PolynomialList;
    
    PolynomialList m_polynomials;
    double mp_params[8];
};

// File names for the nasa polynomial databases
char nasa7[] = "nasa7.dat";
char nasa9[] = "nasa9.dat";

// Register the NASA-7 and NASA-9 thermodynamic databases as ThermoDB types
Utilities::ObjectProvider<NasaDB<Nasa7Polynomial, nasa7>, ThermoDB> nasa7DB("NASA-7");
Utilities::ObjectProvider<NasaDB<Nasa9Polynomial, nasa9>, ThermoDB> nasa9DB("NASA-9");


