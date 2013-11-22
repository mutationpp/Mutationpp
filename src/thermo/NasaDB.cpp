
#include "NasaDB.h"
#include "Utilities.h"

#include <string>

namespace Mutation {
    namespace Thermodynamics {


//    NasaDB(int 0)
//        : ThermoDB(298.0, 1.0E5)
//    {
//        // Open the database file
//        std::string thermo_directory = 
//            getEnvironmentVariable("MPP_DATA_DIRECTORY") + "/thermo";
//        std::string db_path = thermo_directory + "/" + std::string(DBNAME);
//        std::ifstream file(db_path.c_str());
//        
//        if (!file.is_open()) {
//            std::cerr << "Could not open " << db_path << "!\n" << std::endl;
//            exit(1);
//        }
//        
//        // Make a map of the species names corresponding to the species list
//        // that are linked to their position in the list
//        std::map<std::string,int> species_names;
//        std::map<std::string,int>::iterator iter;
//        for (int i = 0; i < species.size(); ++i)
//            species_names[species[i].nasaName(POLYNOMIAL::nCoefficients())] = i;
//            
//        std::string line, name;
//        
//        // Read everything in the file or load all species
//        while (std::getline(file, line)) {
//            // Get the first "word" in the line which could be a species name
//            name = line.substr(0, line.find(" "));
//            
//            // "END" signifies the end of the thermodynamic data
//            if (name == "END")
//                break;
//            
//            if ((iter = species_names.find(name)) != species_names.end()) {
//                // Store the polynomial in the vector and the species position
//                // in the species list
//                m_polynomials.push_back(std::pair<int, POLYNOMIAL>(
//                        iter->second, POLYNOMIAL()));
//                
//                // Load the polynomial
//                file.seekg(
//                    -static_cast<int>(line.length()+1), std::ios_base::cur) >> 
//                    m_polynomials.back().second;
//                
//                // Remove the species from the list of species left to load
//                species_names.erase(iter);
//                
//                // If the list is now empty then we don't need to look anymore
//                if (species_names.empty())
//                    break;
//            }
//        }
//        
//        // Be sure to close the file afterward
//        file.close();
//        
//        // Check and see if we got all the species we needed
//        if (species_names.empty())
//            return;
//        
//        // If we get here then not all required species were loaded
//        std::cout << "Could not find all required species in the NASA polynomial "
//                  << "database (" << DBNAME << ")!" << std::endl;
//        std::cout << "Missing species:" << std::endl;
//        
//        std::map<std::string, int>::const_iterator missing_iter = 
//            species_names.begin();
//        
//        for ( ; missing_iter != species_names.end(); ++missing_iter)
//            std::cout << std::setw(10) << missing_iter->first << std::endl;
//        
//        std::exit(1);
//    }

    

    
    

//private:
   
    //typedef std::vector<std::pair<int, POLYNOMIAL> > PolynomialList;
    
    //PolynomialList m_polynomials;
    //double mp_params[8];


// File names for the nasa polynomial databases
//char nasa7[] = "nasa7.dat";
//char nasa9[] = "nasa9.dat";

// Register the NASA-7 and NASA-9 thermodynamic databases as ThermoDB types
//Config::ObjectProvider<NasaDB<Nasa7Polynomial, nasa7>, ThermoDB> nasa7DB("NASA-7");
//Config::ObjectProvider<NasaDB<Nasa9Polynomial, nasa9>, ThermoDB> nasa9DB("NASA-9");

    } // namespace Thermodynamics
} // namespace Mutation


