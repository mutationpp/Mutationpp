
#include "NasaDB.h"
#include "Nasa9Polynomial.h"
#include "Utilities.h"

#include <cstdlib>
#include <string>
#include <iostream>

namespace Mutation {
    namespace Thermodynamics {
    
using namespace std;
using namespace Utilities;

class Nasa9DB : public NasaDB<Nasa9Polynomial>
{
public:
    Nasa9DB(int arg)
    { }

protected:
    /**
     * Returns the name of the file where this database resides.
     */
    std::string filename() const {
        return std::string("nasa9.dat");
    }

    /**
     * Positions the stream at the beginning of the species data in the 
     * database file.
     */
    void skipHeader(std::ifstream& is) const
    {
        // Start by skipping all of the comments
        std::string line;
        int length = 0;
        while (std::getline(is, line)) {
            length = line.length();
            line = String::trim(line);
            if (line != "" && line[0] != '!') break;
        }
        
        // Skip the "thermo" header if it exists
        if (length > 6 && String::toLowerCase(line.substr(0,6)) == "thermo")
            std::getline(is, line);
        // Else put the curser back to the beginning of the line
        else
            is.seekg(-static_cast<int>(length+1), std::ios_base::cur);
    }
    
    /**
     * Loads the species that is at the current location in the file stream.
     */
    Species loadSpecies(
        std::ifstream& is, const std::vector<Element>& elements) const
    {
        // Skip all comments and blank lines first
        std::string line;
        while (std::getline(is, line))
            if (String::trim(line).length() > 1 && line[0] != '!') break; 
        
        if (is.eof() || String::toLowerCase(line.substr(0,3)) == "end")
            return Species();
        
        // Species name
        std::string name = Utilities::String::trim(line.substr(0,24));

        // Phase
        getline(is, line);
        PhaseType phase;
        switch (line[51]) {
            case '0':
                phase = GAS;
                break;
            default:
                if (name.substr(name.size()-3,3) == "(L)")
                    phase = LIQUID;
                else
                    phase = SOLID;
        }
        
        // Stoichiometry
        std::vector< std::pair<std::string, int> > stoich; 
        for (size_t i = 10; i < 50; i += 8) {
            // Element name
            std::string el = Utilities::String::trim(line.substr(i,2));
            if (el == "") break;
            
            if (el == "E" || el == "E-" || el == "e")
                el = "e-";
            else {
                el[0] = std::toupper(el[0]);
                el[1] = std::tolower(el[1]);
            }
            
            // Number of atoms
            int n = (int)std::atof(line.substr(i+2,6).c_str());
            stoich.push_back(make_pair(el, n));
        }
        
        // Skip over the polynomial data
        int ranges = std::atoi(line.substr(0,2).c_str());
        
        if (ranges == 0)
            std::getline(is, line);
        else {
            for (int i = 0; i < ranges; ++i) {
                std::getline(is, line);
                std::getline(is, line);
                std::getline(is, line);
            }
        }
        
        return Species(name, phase, stoich, elements);
    }
    
    /**
     * Loads the thermodynamic data corresponding to the required species.  In
     * this case, the data are the 7 coefficient polynomials.
     */
    void loadPolynomials(
        std::ifstream& is, std::vector<Nasa9Polynomial>& polynomials)
    {
        // First build a map of the species names with their indices in order
        // to speed up the searching for needed names
        std::map<std::string, size_t> species_names;
        std::map<std::string, size_t>::iterator iter;
        for (size_t i = 0; i < species().size(); ++i)
            species_names.insert(make_pair(species()[i].name(), i));
        
        // Keep reading species from file until we have found all the ones we
        // need
        std::string line, name;
        while (species_names.size() > 0) {
            // Read each line until encountering a species that we want
            std::getline(is, line);
            
            // Is this the first line of a species that we need?
            iter = species_names.find(
                Utilities::String::trim(line.substr(0,24)));                
            if (iter != species_names.end()) {
                is.seekg(
                    -static_cast<int>(line.length()+1), std::ios_base::cur) >>
                    polynomials[iter->second];
                species_names.erase(iter);
            }
        }
    }
};

// Register this database type
Utilities::Config::ObjectProvider<Nasa9DB, ThermoDB> nasa9DB("NASA-9");


    } // namespace Thermodynamics
} // namespace Mutation
