#include <iostream>
#include <fstream>
#include <iomanip>
#include <typeinfo>
//#include <random> c++x11
#include <cstdlib>
#include <ctime>
#include <vector>
#include <set>
#include <cmath>
#include <algorithm>

#include "mutation++.h"

using namespace Mutation;

#define NCOMPOSITIONS 1000

#define T_MIN     500.0
#define T_MAX    5000.0
#define T_POINTS   200

#define P_MIN     (10.0)
#define P_MAX     (10.0*ONEATM)
#define P_POINTS   20

/**
 * Simple class to represent a reduced species list at a given tolerance.
 */
class ReducedSet
{
public:
    
    ReducedSet(double tol)
        : m_tol(tol)
    { }
    
    double tolerance() const {
        return m_tol;
    }
    
    const std::set<int>& species() const {
        return m_species;
    }
    
    void updateFromComposition(const std::vector<double>& xs)
    {
        for (int i = 0; i < xs.size(); ++i)
            if (xs[i] >= m_tol) m_species.insert(i);
    }
    
    std::ostream& writeReducedSet(std::ostream& os, const Mixture& mix) const
    {    
        std::set<int>::const_iterator iter = m_species.begin();
        for ( ; iter != m_species.end(); ++iter)
            os << mix.speciesName(*iter) << " ";
        os << std::endl;
        
        return os;
    }

private:

    double m_tol;
    std::set<int> m_species;
};

/**
 * Computes a random composition using a uniform distribution and ensuring mole
 * fractions sum to 1.
 */
void randomComposition(std::vector<double>& xe, const Mixture& mix)
{
    // Compute a random element fraction distribution
    double sum = 0.0;
    for (int k = 0; k < mix.nElements(); ++k) {
        xe[k] = (mix.elementName(k) == "e-" ? 0.0 :
            double(rand()) / double(RAND_MAX));// * (1.0-sum));
        sum += xe[k];
    }
    //xe[mix.nElements()-1] = std::max(1.0-sum, 0.0);
    for (int k = 0; k < mix.nElements(); ++k)
        xe[k] /= sum;
}

/**
 * Performs the reduction for a given temperature and pressure.
 */
void reduceTP(
    double T, double P, Mixture& mix, std::vector<ReducedSet>& reductions)
{
    const int ne = mix.nElements();
    const int ns = mix.nSpecies();
    const int mod = (NCOMPOSITIONS > 20 ? NCOMPOSITIONS/20 : 1);
    
    std::vector<double> xe(ne);
    std::vector<double> xs(ns);
    
    for (int i = 0; i < NCOMPOSITIONS; ++i) {
        if (i % mod == 0)
            cout << i << " "; cout.flush();
        
        randomComposition(xe, mix);
        
        //for (int k = 0; k < ne; ++k)
        //    cout << mix.elementName(k) << ": " << xe[k] << endl;
        
        // Compute equilibrium mixture
        mix.equilibriumComposition(T, P, &xe[0], &xs[0]);
        
        // Update reductions
        for (int k = 0; k < reductions.size(); ++k)
            reductions[k].updateFromComposition(xs);
    }
    cout << endl;
}

void writeLatexTable(
    std::ostream& os, const std::set<int>& species, const Mixture& mix)
{
    std::vector<std::string> tokens;

    std::set<int>::const_iterator iter = species.begin();
    for (int i = 0 ; iter != species.end(); ++iter, ++i) {
        Utilities::String::tokenize(mix.speciesName(*iter), tokens, ",_");
        os << "\\ce{" << tokens[0] << "}";
        if (tokens.size() == 2)
            os << " (" << tokens[1] << ")";
        if (i == species.size()-1 || ((i + 1) % 5 == 0))
            os << "\\\\" << std::endl;
        else
            os << " & ";
    }
    os << std::endl;
}

void computeErrorTable(
    const ReducedSet& reduction, Mixture& full, MixtureOptions opts)
{
    // Create a new mixture based on the reduced set of species
    opts.setMechanism("none");
    
    std::string species_names("");
    std::set<int>::const_iterator iter = reduction.species().begin();
    for ( ; iter != reduction.species().end(); ++iter)
        species_names += " " + full.speciesName(*iter);
    
    opts.setSpeciesDescriptor(species_names);
    Mixture mix(opts);
    
    for (int k = 0; k < mix.nElements(); ++k)
        assert(mix.elementName(k) == full.elementName(k));
    
    // Now loop over a range of compositions, temperatures and pressures to
    // compute the error in Cp, H, S, G
    // Loop over temperature and pressure values
    double TP [2];
    double cpr, cpf, cp_error = 0.0;
    double hr, hf, h_error = 0.0;
    double sr, sf, s_error = 0.0;
    double gr, gf, g_error = 0.0;
    
    std::vector<double> xe(mix.nElements());
    for (int t = 0; t < T_POINTS; ++t) {
        TP[0] = double(t)/double(T_POINTS-1)*(T_MAX-T_MIN)+T_MIN;
        
        for (int p = 0; p < P_POINTS; ++p) {
            TP[1] = std::exp(double(p)/double(P_POINTS-1)*std::log(P_MAX/P_MIN)+
                std::log(P_MIN));
            
            for (int i = 0; i < std::max(NCOMPOSITIONS / 10, 100); ++i) {
                randomComposition(xe, mix);
                
                mix.setState(TP, &xe[0]);
                cpr = mix.mixtureEquilibriumCpMass();
                hr  = mix.mixtureHMass();
                sr  = mix.mixtureSMass();
                //gr  = mix.mixtureGMass();
                
                full.setState(TP, &xe[0]);
                cpf = full.mixtureEquilibriumCpMass();
                hf  = full.mixtureHMass();
                sf  = full.mixtureSMass();
                //gf  = full.mixtureGMass();
                
                cp_error = std::max(cp_error, std::abs((cpf-cpr)/cpf));
                h_error  = std::max(h_error, std::abs((hf-hr)/hf));
                s_error  = std::max(s_error, std::abs((sf-sr)/sf));
                //g_error  = std::max(g_error, std::abs((gf-gr)/gf));
            }
        }
    }
    
    // Now print the results
    cout << setw(10) << reduction.tolerance();
    cout << setw(14) << cp_error*100;
    cout << setw(14) << h_error*100;
    cout << setw(14) << s_error*100;
    //cout << setw(14) << g_error;
    cout << endl;

}

void reduceRandomComposition(Mixture& mix, std::vector<ReducedSet>& reductions, std::vector<double>& xe)
{
    const int ne = mix.nElements();
    const int ns = mix.nSpecies();
    const int ncases = T_POINTS*P_POINTS;

    std::vector<double> xs(ns);

    // Get a random composition
    randomComposition(xe, mix);

    // Loop over temperature and pressure values
    double T, P;
    int tpcase = 0;
    for (int t = 0; t < T_POINTS; ++t) {
        T = double(t)/double(T_POINTS-1)*(T_MAX-T_MIN)+T_MIN;

        for (int p = 0; p < P_POINTS; ++p, ++tpcase) {
            P = std::exp(double(p)/double(P_POINTS-1)*std::log(P_MAX/P_MIN)+
                std::log(P_MIN));

            // Compute the composition
            mix.equilibriumComposition(T, P, &xe[0], &xs[0], Thermodynamics::GLOBAL);

            // Update reductions
            for (int k = 0; k < reductions.size(); ++k)
                reductions[k].updateFromComposition(xs);

            if (tpcase % (ncases/20) == 0) {
                cout << ".";
                cout.flush();
            }

        }
    }

    cout << ">";
}

int main(int argc, char** argv)
{
    MixtureOptions* p_opts;
    
    if (argc < 2 || argc > 3) {
        cout << "- reduce mixture-name" << endl;
        cout << "           or          " << endl;
        cout << "- reduce database(NASA-7, NASA-9, RRHO) species-descriptor" << endl;
        exit(1);
    } else if (argc == 2) {
        p_opts = new MixtureOptions(argv[1]);
    } else {
        p_opts = new MixtureOptions();
        p_opts->setThermodynamicDatabase(argv[1]);
        p_opts->setSpeciesDescriptor(argv[2]);
    }
    
    p_opts->setStateModel("EquilTP");
    Mixture mix(*p_opts);
    delete p_opts;
    
    srand(time(NULL));
    
    // Define the different tolerance levels
    std::vector<ReducedSet> reductions;
    reductions.push_back(1.0E-2);
    reductions.push_back(1.0E-3);
    reductions.push_back(1.0E-4);
    reductions.push_back(1.0E-5);
    reductions.push_back(1.0E-7);
    reductions.push_back(1.0E-10);
    
    // Save the history of the results
    std::ofstream file("hist.dat");
    // Header
    file << "#    ";
    for (int i = 0; i < mix.nElements(); ++i)
        file << setw(15) << mix.elementName(i);
    for (int i = 0; i < reductions.size(); ++i)
        file << setw(10) << reductions[i].tolerance();
    file << endl;


    // Loop over compositions
    std::vector<double> xe(mix.nElements());
    for (int k = 0; k < NCOMPOSITIONS; ++k) {
        cout << "(" << setw(5) << k+1 << "/" << setw(5) << NCOMPOSITIONS << ") ";
        reduceRandomComposition(mix, reductions, xe);
            
        // Write out the current reductions based on the previous
        // equilibrate calls
        for (int i = 0; i < reductions.size(); ++i) {
            cout << setw(7) << reductions[i].tolerance() << " : "
                 << setw(3) << reductions[i].species().size() << "   ";
        }
        cout << endl;

        // Save the history of the results
        file << setw(5) << k;
        for (int i = 0; i < mix.nElements(); ++i)
            file << setw(15) << xe[i];
        for (int i = 0; i < reductions.size(); ++i)
            file << setw(10) << reductions[i].species().size();
        file << endl;

    }
    file.close();
    
    
    // Write out the latex tables starting with the first reduction and then
    // adding species that are not in previous reductions
    file.open("table.tex");
    writeLatexTable(file, reductions[0].species(), mix); file << endl;
    std::set<int> old_species(
        reductions[0].species().begin(), reductions[0].species().end());
    std::set<int> new_species;
    
    for (int i = 1; i < reductions.size(); ++i) {
        std::set_difference(
            reductions[i].species().begin(), reductions[i].species().end(),
            old_species.begin(), old_species.end(),
            std::inserter(new_species, new_species.end()));
        old_species.insert(
            reductions[i].species().begin(), reductions[i].species().end());
        
        writeLatexTable(file, new_species, mix);
        new_species.clear();
        file << endl;
    }
    file.close();
    
    char filename [20];
    for (int i = 0; i < reductions.size(); i++) {
        sprintf(filename, "list-%02d.txt", (int)-log10(reductions[i].tolerance()));
        file.open(filename);
        reductions[i].writeReducedSet(file, mix);
        file.close();
    }
    
    // Compute a random sampling of composition points which can be used to show
    // good coverage of the composition space
    /*std::ofstream file("compositions.dat");
    for (int i = 0; i < mix.nElements(); ++i)
        file << setw(14) << mix.elementName(i);
    file << endl;
    std::vector<double> xe(mix.nElements());
    for (int i = 0; i < NCOMPOSITIONS; ++i) {
        randomComposition(xe, mix);
        for (int k = 0; k < mix.nElements(); ++k)
            file << setw(14) << xe[k];
        file << endl;
    }
    file.close();*/
    
    // Compute the mixture Cp, H, S, and G and determine the maximum error in
    // each for each mixture
//    cout << "Maximum errors thermodynamic properties..." << endl;
//    cout << setw(10) << "Tol";
//    cout << setw(14) << "Cp_mix [%]";
//    cout << setw(14) << "H_mix [%]";
//    cout << setw(14) << "S_mix [%]" << endl;
//    for (int i = 0; i < reductions.size(); ++i) {
//        computeErrorTable(reductions[i], mix, opts);
//    }
}
