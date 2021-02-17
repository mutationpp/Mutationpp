/**
 * @file Nasa9DB.cpp
 *
 * @brief Provides Nasa-9 thermodynamic database.
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
#include "Nasa9Polynomial.h"
#include "Utilities.h"

#include <cstdlib>
#include <string>
#include <iostream>

namespace Mutation {
    namespace Thermodynamics {
    
using namespace std;
using namespace Utilities;

/**
 * @brief Adds support for the new NASA-9 polynomial formatted database.
 */
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
        if (length > 5 && String::toLowerCase(line.substr(0,6)) == "thermo")
            std::getline(is, line);
        // Else put the curser back to the beginning of the line
        else
            is.seekg(-static_cast<int>(length+1), std::ios_base::cur);
    }
    
    /**
     * Loads the species that is at the current location in the file stream.
     */
    Species loadSpecies(std::ifstream& is) const
    {
        // Skip all comments and blank lines first
        std::string line;
        while (std::getline(is, line))
            if (String::trim(line).length() > 1 && line[0] != '!') break; 
        
        if (is.eof() || String::toLowerCase(line.substr(0,3)) == "end")
            return Species();
        
        // Species name
        std::vector<std::string> tokens;
        std::string name =
            Utilities::String::tokenize(line.substr(0,24), tokens, " ")[0];

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
        Species::StoichList stoich;
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
            stoich(el, n);
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
        
        return Species(name, phase, stoich);
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
        std::vector<std::string> tokens;
        while (species_names.size() > 0) {
            // Read each line until encountering a species that we want
            std::getline(is, line);
            
            // Is this the first line of a species that we need?
            tokens.clear();
            name =
                Utilities::String::tokenize(line.substr(0,24), tokens, " ")[0];
            iter = species_names.find(name);
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


/**
 * @brief Adds support for the new NASA-9 polynomial formatted database.
 * 
 * The database is detailed in Scoggins et al. Aerospace Science and Technology 
 * 66:177-192, 2017.  Since it is still being vetted, loading it will display a 
 * warning to the user.
 * 
 * @todo Finalize the new NASA-9 database and provide a better name.
 */
class Nasa9NewDB : public Nasa9DB
{
public:
    Nasa9NewDB(int arg) : Nasa9DB(arg)
    { 
        std::cout 
            << "Warning: the NASA-9-New thermodynamic database is assembled "
            << "from a collection of sources and is still being vetted. See "
            << "Scoggins et al. Aerospace Science and Technology 66:177-192, "
            << "2017. for more details.  The name of this database may change "
            << "in the future." << std::endl;
    }

protected:
    /**
     * Returns the name of the file where this database resides.
     */
    std::string filename() const {
        return std::string("nasa9_new.dat");
    }
};

// Register this database type
Utilities::Config::ObjectProvider<Nasa9NewDB, ThermoDB> nasa9newDB("NASA-9-New");



    } // namespace Thermodynamics
} // namespace Mutation
