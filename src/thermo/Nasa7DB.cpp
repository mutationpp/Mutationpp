/**
 * @file Nasa7DB.cpp
 *
 * @brief Provides Nasa-7 thermodynamic database.
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


#include "NasaDB.h"
#include "Nasa7Polynomial.h"
#include "Utilities.h"

#include <cstdlib>
#include <string>
#include <map>

namespace Mutation {
    namespace Thermodynamics {

/**
 * @brief Adds support for the new NASA-7 polynomial formatted database.
 */
class Nasa7DB : public NasaDB<Nasa7Polynomial>
{
public:
    Nasa7DB(int arg)
    { }

protected:
    /**
     * Returns the name of the file where this database resides.
     */
    std::string filename() const {
        return std::string("nasa7.dat");
    }

    /**
     * Positions the stream at the beginning of the species data in the 
     * database file.
     */
    void skipHeader(std::ifstream& is) const
    {
        std::string line;
        while (std::getline(is, line))
            if (line.substr(0,6) == "THERMO") break;
        std::getline(is, line);
    }
    
    /**
     * Loads the species that is at the current location in the file stream.
     */
    Species loadSpecies(std::ifstream& is) const
    {
        std::string line;
        if (!std::getline(is, line) || line.substr(0,3) == "END")
            return Species();
        
        // Species name
        std::string name = Utilities::String::trim(line.substr(0,18));
        
        // Phase
        PhaseType phase;
        switch (line[44]) {
            case 'G':
                phase = GAS;
                break;
            case 'S':
                phase = SOLID;
                break;
            case 'L':
                phase = LIQUID;
                break;
            default:
                if (name.substr(name.size()-3,3) == "(L)")
                    phase = LIQUID;
                else
                    phase = SOLID;
        }
        
        // Stoichiometry
        Species::StoichList stoich;
        for (size_t i = 24; i < 40; i += 5) {
            // Element name
            std::string el = Utilities::String::trim(line.substr(i,2));
            if (el == "") break;
            
            if (el == "E")
                el = "e-";
            else {
                el[0] = std::toupper(el[0]);
                el[1] = std::tolower(el[1]);
            }
            
            // Number of atoms
            int n = (int)std::atof(line.substr(i+2,3).c_str());
            stoich(el, n);
        }
        
        // Skip the next 3 lines
        std::getline(is, line);
        std::getline(is, line);
        std::getline(is, line);
        
        return Species(name, phase, stoich);
    }
    
    /**
     * Loads the thermodynamic data corresponding to the required species.  In
     * this case, the data are the 7 coefficient polynomials.
     */
    void loadPolynomials(
        std::ifstream& is, std::vector<Nasa7Polynomial>& polynomials)
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
            std::getline(is, line);
            
            // Do we need this species?
            iter = species_names.find(
                Utilities::String::trim(line.substr(0,18)));
            if (iter == species_names.end()) {
                // No, skip remaining lines
                std::getline(is, line);
                std::getline(is, line);
                std::getline(is, line);
            } else {
                // Yes, load the polynomial and remove species from list
                is.seekg(
                    -static_cast<int>(line.length()+1), std::ios_base::cur) >>
                    polynomials[iter->second];
                species_names.erase(iter);
            }
        }
    }
};

// Register this database type
Utilities::Config::ObjectProvider<Nasa7DB, ThermoDB> nasa7DB("NASA-7");


    } // namespace Thermodynamics
} // namespace Mutation
