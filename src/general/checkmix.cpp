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
 * @page checkmix checkmix
 *
 * @tableofcontents
 *
 * This program will load a mixture and print out information about the various
 * elements, species, and reactions in the [mixture](@ref mixtures).  A mixture
 * is loaded just as it would be in any other application, so this tool is
 * useful to check for any syntax errors or missing data in a given mixture
 * before using it elsewhere.  It can also be used to see
 * [exactly which order species](@ref species_order) and reactions are stored
 * internally in Mutation++ for a given mixture.
 *
 * @section checkmix_usage Usage
 * @subsection checkmix_usage_1 Using a Mixture File
 *
 * `checkmix mixture`
 *
 * @subsection checkmix_usage_2 Using a Species Descriptor
 *
 * `checkmix [NASA-7 | NASA-9 | RRHO] "species-descriptor"`
 *
 * `species-descriptor` should follow the rules given
 * [here](@ref species_list).
 *
 * @section checkmix_example Example
 * Assuming you have a [mixture file](@ref mixtures) called `air5.xml` and a
 * [mechanism file](@ref reaction_mechanisms) called `air5-chem.xml` (in either
 * your local directory or the `data/mixtures` and `data/mechanisms`
 * directories) which are as follows
 *
 * @code{.xml}
 * <mixture name="air5" mechanism="air5-chem">
 *     <species>
 *         N O NO N2 O2
 *     </species>
 * </mixture>
 * @endcode
 *
 * and
 *
 * @code{.xml}
 * <mechanism name="air5">
 *     <arrhenius_units A="mol,cm,s,K" E="kcal,mol,K" />
 *     <!-- 1 -->
 *     <reaction formula="N2+M=2N+M">
 *         <arrhenius A="3.0E+22" n="-1.6" T="113200.0" />
 *        <M>N2:0.2333, NO:0.2333, O2:0.2333</M>
 *     </reaction>
 *     <!-- 2 -->
 *     <reaction formula="O2+M=2O+M">
 *         <arrhenius A="1.0E+22" n="-1.5" T="59360.0" />
 *         <M>N2:0.5, NO:0.5, O2:0.5</M>
 *     </reaction>
 *     <!-- 3 -->
 *     <reaction formula="NO+M=N+O+M">
 *         <arrhenius A="5.0E15" n="+0.0" T="75500.0" />
 *         <M>NO:20.0, N:20.0, O:20.0</M>
 *     </reaction>
 *     <!-- 4 -->
 *     <reaction formula="N2+O=NO+N">
 *         <arrhenius A="5.69E+12" n="+0.42" T="42938.0" />
 *     </reaction>
 *     <!-- 5 -->
 *     <reaction formula="O2+N=NO+O">
 *         <arrhenius A="2.49E+09" n="+1.18" T="4005.5" />
 *     </reaction>
 * </mechanism>
 * @endcode
 *
 * then the command " `checkmix air5` " will produce the following output:
 *
 * @code{.unparsed}
 * 5 species containing 2 elements
 * 5 reactions
 *
 * Species info:
 * -------------
 *     N   O  Mw (g/mol)    Charge       Phase
 * Gas Species (5):
 * N     1   0     14.0067         0         gas
 * O     0   1     15.9994         0         gas
 * NO    1   1     30.0061         0         gas
 * N2    2   0     28.0134         0         gas
 * O2    0   2     31.9988         0         gas
 *
 * Default elemental composition:
 * ------------------------------
 *    N  :   0.5
 *    O  :   0.5
 *
 * Reaction info:
 * --------------
 * Type ID Key
 *    3: heavy particle impact dissociation
 *    7: exchange
 *
 * Reactions
 *    #  Formula             Type  Rate Law     A (m,s,mol)      n    Ta (K)
 *    1: N2+M=2N+M           3     Arrhenius:     3.000e+16  -1.60  113200.0
 *       N2: 0.23, NO: 0.23, O2: 0.23
 *    2: O2+M=2O+M           3     Arrhenius:     1.000e+16  -1.50   59360.0
 *       N2: 0.50, NO: 0.50, O2: 0.50
 *    3: NO+M=N+O+M          3     Arrhenius:     5.000e+09   0.00   75500.0
 *       NO: 20.00, N: 20.00, O: 20.00
 *    4: N2+O=NO+N           7     Arrhenius:     5.690e+06   0.42   42938.0
 *    5: O2+N=NO+O           7     Arrhenius:     2.490e+03   1.18    4005.5
 * @endcode
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
