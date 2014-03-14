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
    for (int i = 0; i < ns; ++i) {
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
