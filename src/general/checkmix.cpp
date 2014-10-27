/**
 * @file checkmix.cpp
 *
 * @brief Loads a mixture and prints information about the mixture to the
 * console. @see @ref checkmix
 */

/*
 * Copyright 2014 von Karman Institute for Fluid Dynamics (VKI)
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

#include <iostream>
#include <iomanip>
#include <typeinfo>

#include "mutation++.h"

using std::cout;
using std::endl;
using std::string;
using std::setw;
using std::vector;
using std::pair;

using namespace Mutation;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Kinetics;

/**
 * @page checkmix
 * __Usage 1__:
 *
 *    checkmix mixture
 *
 * __Usage 2__:
 *
 *    checkmix [NASA-7 | NASA-9 | RRHO] species-descriptor
 *
 * This program will load a mixture and print out information about the various
 * elements, species, and reactions in the mixture.  A mixture is loaded just
 * as it would be in any other application, so this tool is useful to check for
 * any syntax errors or missing data in a given mixture before using it
 * elsewhere.  It can also be used to see exactly which order species and
 * reactions are stored internally in Mutation++ for a given mixture.
 */

int main(int argc, char** argv)
{
    MixtureOptions opts;
    
    if (argc < 2 || argc > 3) {
        cout << "- checkmix mixture-name" << endl;
        cout << "           or          " << endl;
        cout << "- checkmix database(NASA-7, NASA-9, RRHO) species-descriptor" << endl;
        exit(1);
    } else if (argc == 2) {
        opts.loadFromFile(argv[1]);
    } else {
        opts.setThermodynamicDatabase(argv[1]);
        opts.setSpeciesDescriptor(argv[2]);
    }
    
    Mixture mixture(opts);
        
    const int ne = mixture.nElements();
    const int ns = mixture.nSpecies();
    const int nr = mixture.nReactions();
    const int ng = mixture.nGas();
    const int nc = mixture.nCondensed();
 
    cout << "location: " << opts.getSource() << endl;
    cout << ns << " species containing " << ne << " elements" << endl;
    cout << nr << " reactions" << endl;
    cout << endl;
    
    int width = 0;
	for (int i = 0; i < ns; ++i)
		width = std::max(width, (int) mixture.speciesName(i).size());
	width++;

    cout << "Species info:" << endl;
    cout << "-------------" << endl;  
    cout << setw(width) << " ";
    for (int i = 0; i < ne; ++i)
        cout << setw(4) << mixture.elementName(i);
    
    cout << setw(12) << "Mw (g/mol)" << setw(10) << "Charge";
    cout << setw(12) << "Phase" << endl;
    cout << "Gas Species (" << ng << "):" << endl;
    for (int i = 0; i < ng; ++i) {
        cout.setf(std::ios::left, std::ios::adjustfield);
        cout << setw(width) << mixture.speciesName(i);

        cout.setf(std::ios::right, std::ios::adjustfield);
        for (int j = 0; j < ne; ++j)
            cout << setw(4) << mixture.elementMatrix()(i,j);

        cout << setw(12) << mixture.speciesMw(i) * 1000.0
             << setw(10) << mixture.species()[i].charge();

        PhaseType phase = mixture.species()[i].phase();
        cout << setw(12) << (phase == GAS ? "gas" :
            (phase == LIQUID ? "liquid" : "solid")) << endl;
    }

    if (nc > 0) {
        cout << "Condensed Species (" << nc << "):" << endl;
        for (int i = ng; i < ns; ++i) {
            cout.setf(std::ios::left, std::ios::adjustfield);
            cout << setw(width) << mixture.speciesName(i);

            cout.setf(std::ios::right, std::ios::adjustfield);
            for (int j = 0; j < ne; ++j)
                cout << setw(4) << mixture.elementMatrix()(i,j);

            cout << setw(12) << mixture.speciesMw(i) * 1000.0
                 << setw(10) << mixture.species()[i].charge();

            PhaseType phase = mixture.species()[i].phase();
            cout << setw(12) << (phase == GAS ? "gas" :
                (phase == LIQUID ? "liquid" : "solid")) << endl;
        }
    }
    
    cout << endl;
    
    cout << "Default elemental composition:" << endl;
    cout << "------------------------------" << endl;
    for (int i = 0; i < ne; ++i) {
        cout.setf(std::ios::left, std::ios::adjustfield);
        cout << "   " << setw(3) << mixture.elementName(i) << ": ";
        cout.setf(std::ios::right, std::ios::adjustfield);
        cout << setw(5) << mixture.getDefaultComposition(i) << endl;
    }
    
    cout << endl;
    
    if (nr == 0) return 0;
    
    cout << "Reaction info:" << endl;
    cout << "--------------" << endl;
    
    // Make a key for the reaction type
    std::set<ReactionType> types;
    for (int i = 0; i < nr; ++i)
        types.insert(mixture.reactions()[i].type());
    
    cout << "Type ID Key" << endl;
    std::set<ReactionType>::const_iterator iter = types.begin();
    for ( ; iter != types.end(); ++iter)
        cout << setw(4) << *iter << ": " << reactionTypeString(*iter) << endl;
    
    // Reaction table header
    cout << endl;
    cout << "Reactions" << endl;
    cout << setw(4) << "#";
    cout.setf(std::ios::left, std::ios::adjustfield);
    cout << setw(22) << "  Formula";
    cout << setw(6)  << "Type";
    cout << setw(12) << "Rate Law";
    cout.setf(std::ios::right, std::ios::adjustfield);
    cout << setw(12) << "A (m,s,mol)";
    cout << setw(7) << "n";
    cout << setw(10) << "Ta (K)" << endl;
    for (int i = 0; i < nr; ++i) {
        const Reaction& r = mixture.reactions()[i];
        
        // Print reaction formula and number
        cout.setf(std::ios::right, std::ios::adjustfield);
        cout << setw(4) << i+1 << ": ";        
        cout.setf(std::ios::left, std::ios::adjustfield);
        cout << setw(20) << r.formula();
        cout << setw(6)  << r.type();
        
        // Print out rate constants
        if (typeid(*(r.rateLaw())) == typeid(Arrhenius)) {
            const Arrhenius& rate = 
                dynamic_cast<const Arrhenius&>(*(r.rateLaw()));
            cout << setw(12) << "Arrhenius: ";
            
            cout.setf(std::ios::right, std::ios::adjustfield);
            cout.setf(std::ios::scientific, std::ios::floatfield);
            cout.precision(3);
            cout << setw(12) << rate.A();
            cout.setf(std::ios::fixed, std::ios::floatfield);
            cout.precision(2);
            cout << setw(7)  << rate.n();
            cout.precision(1);
            cout << setw(10) << rate.T();
        }
        
        cout << endl;
        
        // If this is a thirdbody reaction, print out thirdbody efficiency
        // factors
        if (r.isThirdbody() && r.efficiencies().size() > 0) {
            
            vector<pair<int, double> >::const_iterator iter =
                r.efficiencies().begin();
            
            cout.precision(2);
            cout << setw(6) << "" << mixture.speciesName(iter->first) << ": "
                 << iter->second;
                
            for (iter++; iter != r.efficiencies().end(); ++iter) {
                cout << ", " << mixture.speciesName(iter->first) << ": "
                     << iter->second;
            }        
        
            cout << endl;
        }
    }
    
    cout << endl;
}
