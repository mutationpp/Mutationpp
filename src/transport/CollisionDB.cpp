#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <iomanip>

#include "Constants.h"
#include "CollisionDB.h"
#include "LookupTable.h"

using namespace std;
using namespace Numerics;

// Charge-charge collision functions
const CollisionFunc5 CollisionDB::sm_Q11_att(
    -7.9270465e-01,5.8867723e-01,-8.7607125e-02,9.3018130e-03,-4.1315208e-04);

const CollisionFunc5 CollisionDB::sm_Q11_rep(
    -1.3980752e+00,8.0482070e-01,-9.4801647e-02,5.2812176e-03,-8.2652059e-05);

const CollisionFunc5 CollisionDB::sm_Q22_att(
    -8.1145738e-01,7.0419264e-01,-1.2219724e-01,1.3234707e-02,-5.6994085e-04);

const CollisionFunc5 CollisionDB::sm_Q22_rep(
    -1.1089170e+00,7.7460857e-01,-1.0168153e-01,6.7600878e-03,-1.5622212e-04);

const CollisionFunc5 CollisionDB::sm_Bst_att(
    2.8616335e-01,-4.8858107e-02,1.3372770e-03,5.6686310e-04,-4.4964175e-05);

const CollisionFunc5 CollisionDB::sm_Bst_rep(
    3.2606973e-01,-2.3742462e-02,-9.6987749e-03,1.7799200e-03,-8.3938754e-05);

const CollisionFunc5 CollisionDB::sm_Cst_att(
    -6.4672764e-01,-1.0737242e-01,1.7962157e-02,-1.8653043e-03,8.0134495e-05);

const CollisionFunc5 CollisionDB::sm_Cst_rep(
    -5.0077127e-01,-1.0639836e-01,-1.1478149e-03,1.8355926e-03,-1.1830337e-04);


CollisionDB::CollisionDB(const Thermodynamics& thermo)
    : m_ns(thermo.nSpecies()), m_ncollisions((m_ns*(m_ns + 1))/2), 
      m_mass(m_ns), m_mass_sum(m_ns), m_mass_prod(m_ns), m_red_mass(m_ns),
      m_em_index(thermo.speciesIndex("e-")), m_Q11(m_ns), m_Q22(m_ns), 
      m_Ast(m_ns), m_Bst(m_ns), m_Cst(m_ns), m_eta(m_ns), m_Dij(m_ns)
{
    // Load collision integrals
    loadCollisionIntegrals(thermo.species());
    
    //LookupTable<double, double> table(
    //    200.0, 20000.0, CollisionDB::computeQ11, this, m_ncollisions);
    
    // Next, compute the mass quantities
    for (int i = 0; i < m_ns; ++i)
        m_mass(i) = thermo.speciesMw(i) / NA;
    
    for (int i = 0; i < m_ns; ++i) {
        for (int j = i; j < m_ns; ++j) {
            m_mass_sum(i,j)  = m_mass(i) + m_mass(j);
            m_mass_prod(i,j) = m_mass(i) * m_mass(j);
            m_red_mass(i,j) = m_mass_prod(i,j) / m_mass_sum(i,j);
        }
    }
    
    for (int i = 0; i < DATA_SIZE; ++i)
        mp_last_T[i] = -1.0;
}

void CollisionDB::loadCollisionIntegrals(const vector<Species>& species)
{
    // First step is to determine all of the collision pairs that are needed and
    // what index they belong to in the collision function lists
    map<CollisionPair, int> collision_map;
    int index = 0;
    
    for (int i = 0; i < m_ns; ++i) {
        for (int j = i; j < m_ns; ++j, ++index)
             collision_map[
                CollisionPair(species[i], species[j])
             ] = index;
    }
    
    // With the collision map generated, look through the database and load all
    // collisions that are found in the map
    string transport_dir = 
        utils::getEnvironmentVariable("MPP_DATA_DIRECTORY") + "/transport";
    string heavy_path = transport_dir + "/heavy.dat";
    ifstream file(heavy_path.c_str(), ios::in);
    
    if (!file.is_open()) {
        cerr << "Could not open file " + heavy_path + "!\n" << endl;
        exit(1);
    }
    
    string str;    
    file >> str;
    
    int ineutral = 0;
    CollisionFunc4 func;
    std::map<CollisionPair, int>::iterator iter;
    
    while (str != "STOP" && !collision_map.empty()) {
        // Check if we have landed on a collision identifier
        if (isValidCollisionString(str)) {
            // If so, is this a collision we need to load?
            if ((iter = collision_map.find(str)) != collision_map.end()) {
                // Load Q11, Q22, and B* functions
                file >> func;
                m_Q11_funcs.push_back(func);
                
                file >> func;
                m_Q22_funcs.push_back(func);
                
                file >> func;
                m_Bst_funcs.push_back(func);
                
                // Add the collision to the list of found indices
                m_neutral_indices.push_back(iter->second);
                
                // Remove collision from collision map to speed up search of
                // remaining collisions
                collision_map.erase(iter);
            }
        }
        
        // Read next record
        file >> str;
    }
    
    file.close();
    
    // Charge-charge collisions that were not given explicit functions in the 
    // database should use integrals computed with a Coulomb potential screened
    // with the Debye length.  These are split up into attractive and repulsive
    // collisions.
    if (!collision_map.empty()) {
        iter = collision_map.begin();
        while (iter != collision_map.end()) {
            if (iter->first.charge() == 0)
                iter++;
            else {
                if (iter->first.charge() < 0)
                    m_attract_indices.push_back(iter->second);
                else
                    m_repulse_indices.push_back(iter->second);
                collision_map.erase(iter++);
            }
        }
    }
    
    // If there are still collisions left over at this point then let the user 
    // know that they will be represented as zeros
    if (!collision_map.empty()) {
        cout << endl;
        cout << "The following collision pairs were not found!" << endl;
        //func = CollisionFunc4();
        iter = collision_map.begin();
        for ( ; iter != collision_map.end(); ++iter) {
            cout << "\t" << iter->first.name() << endl;
            m_Q11_funcs.push_back(m_Q11_funcs[0]); // they are not evaluated as
            m_Q22_funcs.push_back(m_Q22_funcs[0]); // zero anymore...
            m_Bst_funcs.push_back(m_Bst_funcs[0]);
            m_neutral_indices.push_back(iter->second);
        }
        cout << "They will be evaluated as the first integral..." << endl;
        cout << endl;
    }
    
    // Find collision pairs amongst those that were loaded which best resemble
    // those that were not loaded
    /*int size_found = m_neutral_indices.size();
    iter = collision_map.begin();
    for ( ; iter != collision_map.end(); ++iter) {
        for (int i = 0; i < size_found; ++i) {
            // Compute the similarity of this collision pair with the one that
            // wasn't found
            
        }
    }*/
}

bool CollisionDB::isValidCollisionString(const std::string &str)
{
    // Must be at least 5 characters long
    if (str.length() < 5)
        return false;
    
    // Must have periods bounding the string but nowhere else
    if (str[0] != '.' || str.find(".", 1) != str.length() - 1)
        return false;
        
    // Must have dash somewhere in the middle
    if (str.find("-", 2) > str.length() - 3)
        return false;
    
    // No spaces
    if (str.find(" ") != string::npos)
        return false;
    
    // Assume that above constraints are good enough for now...
    // Note we could use regex functions for a more accurate analysis
    // later when it is added to the standard library
    return true;
}

void CollisionDB::updateCollisionData(
    const double Th, const double Te, const double nd, const double *const p_x, 
    const CollisionData data)
{
    // Return if we already computed this stuff for the same temperature
    if (std::abs(Th - mp_last_T[data]) < std::sqrt(RealConsts::eps))
        return;
    //if (mp_last_T[data] > 0.0) return;
    
    // Update the needed values
    const double lnT = log(Th);
    
    const size_t nn = m_neutral_indices.size();
    const size_t na = m_attract_indices.size();
    const size_t nr = m_repulse_indices.size();
    
    // Average closest impact parameters
    const double bfac = QE * QE / (8.0 * PI * EPS0 * KB);
    const double be   = bfac / Te;   // electron-electron
    const double bh   = bfac / Th;   // ion-ion
    
    // Compute quantities needed for the Debye-Huckel potential integrals
    // Debye length (set to zero if there are no electrons)
    const double lambdaD = 
        (m_em_index >= 0 ?
        std::min(
            std::sqrt(EPS0 * KB * Te / 
                (2.0 * nd * std::max(p_x[m_em_index], 1.0e-99) * QE * QE)),
            10000.0 * (be + bh)) :
        0.0); 
    
    // Reduced temperatures
    const double Tste = std::max(0.5 * lambdaD / be, 0.1);
    const double Tsth = std::max(0.5 * lambdaD / bh, 0.1);
    
    // Factors used in the curve-fitting of the reduced collision integrals
    const double efac   = PI * lambdaD * lambdaD / (Tste * Tste);
    const double hfac   = PI * lambdaD * lambdaD / (Tsth * Tsth);
    const double lnTste = log(Tste);
    const double lnTsth = log(Tsth);
    
    switch (data) {
        
        case Q11IJ: {
            // Neutral collisions
            for (int i = 0; i < nn; ++i)
                m_Q11(m_neutral_indices[i]) = 1.0E-20 * m_Q11_funcs[i](lnT);
                
            // Charged collisions
            const double Q11_rep = hfac * sm_Q11_rep(lnTsth);
            for (int i = 0; i < nr; ++i)
                m_Q11(m_repulse_indices[i]) = Q11_rep;
            
            const double Q11_att = hfac * sm_Q11_att(lnTsth);
            for (int i = 0; i < na; ++i)
                m_Q11(m_attract_indices[i]) = Q11_att;
            } break;
        
        case Q22IJ: {
            for (int i = 0 ; i < nn; ++i)
                m_Q22(m_neutral_indices[i]) = 1.0E-20 * m_Q22_funcs[i](lnT);
            
            // Charged collisions
            const double Q22_rep = hfac * sm_Q22_rep(lnTsth);
            for (int i = 0; i < nr; ++i)
                m_Q22(m_repulse_indices[i]) = Q22_rep;
            
            const double Q22_att = hfac * sm_Q22_att(lnTsth);
            for (int i = 0; i < na; ++i)
                m_Q22(m_attract_indices[i]) = Q22_att;
            } break;
            
        case ASTAR:
            updateCollisionData(Th, Te, nd, p_x, Q11IJ);
            updateCollisionData(Th, Te, nd, p_x, Q22IJ);
            for (int i = 0; i < m_ncollisions; ++i)
                m_Ast(i) = m_Q22(i) / m_Q11(i);
            break;
                
        case BSTAR: {
            // Neutral collisions from the database
            for (int i = 0; i < nn; ++i)
                m_Bst(m_neutral_indices[i]) = m_Bst_funcs[i](lnT);
            
            // Repulsive collisions
            const double Bst_rep = hfac * sm_Bst_rep(lnTsth);
            for (int i = 0; i < nr; ++i)
                m_Bst(m_repulse_indices[i]) = Bst_rep;
            
            // Attractive collisions
            const double Bst_att = hfac * sm_Bst_att(lnTsth);
            for (int i = 0; i < na; ++i)
                m_Bst(m_attract_indices[i]) = Bst_att;
            } break;
            
        case CSTAR: {        
            // Neutral-electron collisions are treated as 1 for now
            m_Cst = 1.0;
            
            // Ion-electron (repulsive) interactions
            int index = 0;
            int j = m_repulse_indices[index];
            const double Cst_rep = hfac * sm_Cst_rep(lnTsth);
            while (j < m_ns) {
                m_Cst(j) = Cst_rep;
                j = m_repulse_indices[++index];
            }
            
            // Ion-electron (attractive) interactions
            index = 0;
            j = m_attract_indices[index];
            const double Cst_att = hfac * sm_Cst_rep(lnTsth);
            while (j < m_ns) {
                m_Cst(j) = Cst_att;
                j = m_attract_indices[++index];
            }
            } break;
        
        case ETAI: {
            updateCollisionData(Th, Te, nd, p_x, Q22IJ);
            
            // electron does not contribute to viscosity
            m_eta(0) = 0.0;
            
            // heavy species
            for (int i = 1; i < m_ns; ++i)
                m_eta(i) = 
                    5.0 / 16.0 * sqrt(PI * KB * Th * m_mass(i)) / m_Q22(i,i);
            } break;
            
        case NDIJ: {
            updateCollisionData(Th, Te, nd, p_x, Q11IJ);
            
            // electron-electron
            m_Dij(0) = 3.0 / 8.0 * sqrt(PI * KB * Te / m_mass(0)) / m_Q11(0);
            
            // electron-heavy interactions
            for (int i = 1; i < m_ns; ++i)
                m_Dij(i) = 
                    3.0 / 16.0 * sqrt(TWOPI * KB * Te / m_mass(0)) / m_Q11(i);
            
            // heavy-heavy interactions
            for (int i = m_ns; i < m_ncollisions; ++i)
                m_Dij(i) = 3.0 / 16.0 * sqrt(TWOPI * KB * Th / 
                    m_red_mass(i)) / m_Q11(i);
            } break;
    }
    
    mp_last_T[data] = Th;
}


