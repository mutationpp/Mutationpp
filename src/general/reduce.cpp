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

#define XE_POINTS 1000

#define T_MIN     300.0
#define T_MAX    6000.0
#define T_POINTS   101

#define P_MIN     (0.001*ONEATM)
#define P_MAX     (10.0*ONEATM)
#define P_POINTS   10

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
        os << "Tolerance: " << m_tol << " (" << m_species.size()
           << " species)" << std::endl;
    
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
    const int mod = XE_POINTS/20;
    
    std::vector<double> xe(ne);
    std::vector<double> xs(ns);
    
    for (int i = 0; i < XE_POINTS; ++i) {
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
    opts.clearSpeciesNames();
    opts.setMechanism("none");
    std::set<int>::const_iterator iter = reduction.species().begin();
    for ( ; iter != reduction.species().end(); ++iter)
        opts.addSpeciesName(full.speciesName(*iter));
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
            
            for (int i = 0; i < std::max(XE_POINTS / 10, 100); ++i) {
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

int main(int argc, char** argv)
{
    if (argc <= 1) {
        cout << "Provide a reaction mechanism." << endl;
        exit(1);
    }
    
    srand(time(NULL));
    
    // Load the mixture to be reduced
    MixtureOptions opts(argv[1]);
    opts.setStateModel("EquilTPXe");
    Mixture mix(opts);
    
    // Define the different tolerance levels
    std::vector<ReducedSet> reductions;
    reductions.push_back(1.0E-3);
    reductions.push_back(1.0E-4);
    reductions.push_back(1.0E-5);
    reductions.push_back(1.0E-7);
    
    // Loop over temperature and pressure values
    double T, P;
    for (int t = 0; t < T_POINTS; ++t) {
        T = double(t)/double(T_POINTS-1)*(T_MAX-T_MIN)+T_MIN;
        
        for (int p = 0; p < P_POINTS; ++p) {
            P = std::exp(double(p)/double(P_POINTS-1)*std::log(P_MAX/P_MIN)+
                std::log(P_MIN));
            cout << "T = " << T << ", P = " << P << endl;
            
            // Perform the reduction at the given T and P
            reduceTP(T, P, mix, reductions);
            cout << endl;
            
            // Write out the current reductions based on the previous
            // equilibrate calls
            for (int i = 0; i < reductions.size(); ++i)
                reductions[i].writeReducedSet(cout, mix) << endl;
            cout << endl;
        }
    }
    
    
    // Write out the latex tables starting with the first reduction and then
    // adding species that are not in previous reductions
    writeLatexTable(cout, reductions[0].species(), mix); cout << endl;
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
        
        writeLatexTable(cout, new_species, mix);
        new_species.clear();
        cout << endl;
    }
    cout << endl;
    
    // Compute a random sampling of composition points which can be used to show
    // good coverage of the composition space
    std::ofstream file("compositions.dat");
    for (int i = 0; i < mix.nElements(); ++i)
        file << setw(14) << mix.elementName(i);
    file << endl;
    std::vector<double> xe(mix.nElements());
    for (int i = 0; i < XE_POINTS; ++i) {
        randomComposition(xe, mix);
        for (int k = 0; k < mix.nElements(); ++k)
            file << setw(14) << xe[k];
        file << endl;
    }
    file.close();
    
    // Compute the mixture Cp, H, S, and G and determine the maximum error in
    // each for each mixture
    cout << "Maximum errors thermodynamic properties..." << endl;
    cout << setw(10) << "Tol";
    cout << setw(14) << "Cp_mix [%]";
    cout << setw(14) << "H_mix [%]";
    cout << setw(14) << "S_mix [%]" << endl;
    for (int i = 0; i < reductions.size(); ++i) {
        computeErrorTable(reductions[i], mix, opts);
    }
}
