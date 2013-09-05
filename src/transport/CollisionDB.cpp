#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <iomanip>

#include "Constants.h"
#include "CollisionDB.h"

//#define VERBOSE

using namespace std;
using namespace Mutation::Numerics;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities;

namespace Mutation {
    namespace Transport {

using Mutation::Thermodynamics::Thermodynamics;

// Charge-charge collision functions
const CollisionFunc5 CollisionDB::sm_Q11_att(
    -7.9270465e-01,5.8867723e-01,-8.7607125e-02,9.3018130e-03,-4.1315208e-04);

const CollisionFunc5 CollisionDB::sm_Q11_rep(
    -1.3980752e+00,8.0482070e-01,-9.4801647e-02,5.2812176e-03,-8.2652059e-05);

const CollisionFunc5 CollisionDB::sm_Q14_att(
    -2.3064366e+00,3.9187956e-01,-5.0541034e-02,4.9996677e-03,-2.1443795e-04);

const CollisionFunc5 CollisionDB::sm_Q14_rep(
    -2.6061028e+00,5.7602188e-01,-8.9813364e-02,8.4176238e-03,-3.1611961e-04);

const CollisionFunc5 CollisionDB::sm_Q15_att(
    -2.6207402e+00,3.6889858e-01,-4.5983762e-02,4.4531740e-03,-1.8883719e-04);

const CollisionFunc5 CollisionDB::sm_Q15_rep(
    -2.8756074e+00,5.3793608e-01,-8.6195947e-02,8.5216824e-03,-3.3719465e-04);

const CollisionFunc5 CollisionDB::sm_Q22_att(
    -8.1145738e-01,7.0419264e-01,-1.2219724e-01,1.3234707e-02,-5.6994085e-04);

const CollisionFunc5 CollisionDB::sm_Q22_rep(
    -1.1089170e+00,7.7460857e-01,-1.0168153e-01,6.7600878e-03,-1.5622212e-04);

const CollisionFunc5 CollisionDB::sm_Q24_rep(
    -1.7664116e+00,6.4675076e-01,-9.8071543e-02,8.6280547e-03,-3.0256691e-04);

const CollisionFunc5 CollisionDB::sm_Bst_att(
    2.8616335e-01,-4.8858107e-02,1.3372770e-03,5.6686310e-04,-4.4964175e-05);

const CollisionFunc5 CollisionDB::sm_Bst_rep(
    3.2606973e-01,-2.3742462e-02,-9.6987749e-03,1.7799200e-03,-8.3938754e-05);

const CollisionFunc5 CollisionDB::sm_Cst_att(
    -6.4672764e-01,-1.0737242e-01,1.7962157e-02,-1.8653043e-03,8.0134495e-05);

const CollisionFunc5 CollisionDB::sm_Cst_rep(
    -5.0077127e-01,-1.0639836e-01,-1.1478149e-03,1.8355926e-03,-1.1830337e-04);

const CollisionFunc5 CollisionDB::sm_Est_att(
    -3.8922959e-01,-8.5334423e-02,8.7430000e-03,1.3038899e-04,-4.3583770e-05);
    
const CollisionFunc5 CollisionDB::sm_Est_rep(
    -3.6035379e-01,-7.3192205e-02,1.3517297e-03,1.1249000e-03,-8.3717249e-05);


CollisionDB::CollisionDB(const Thermodynamics& thermo)
    : m_ns(thermo.nSpecies()),
      m_nh(thermo.nHeavy()),
      m_ncollisions((m_ns*(m_ns + 1))/2),
      m_em_index(thermo.speciesIndex("e-")),
      m_mass(m_ns),
      m_mass_sum(m_ns),
      m_mass_prod(m_ns),
      m_red_mass(m_ns),
      m_Q11(m_ns),
      m_Q12ei(m_ns),
      m_Q13ei(m_ns),
      m_Q14ei(m_ns),
      m_Q15ei(m_ns),
      m_Q22(m_ns),
      m_Ast(m_ns),
      m_Bst(m_ns),
      m_Cst(m_ns),
      m_eta(m_ns),
      m_Dij(m_ns)
{
    // Load collision integrals
    loadCollisionIntegrals(thermo.species());
    
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
#ifdef USE_COLLISION_INTEGRAL_TABLES
    std::vector<CollisionFunc4> Q11_funcs;
    std::vector<CollisionFunc4> Q22_funcs;
    std::vector<CollisionFunc4> Bst_funcs;
#else
    std::vector<CollisionFunc4>& Q11_funcs = m_Q11_funcs;
    std::vector<CollisionFunc4>& Q22_funcs = m_Q22_funcs;
    std::vector<CollisionFunc4>& Bst_funcs = m_Bst_funcs;
#endif
    
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
        getEnvironmentVariable("MPP_DATA_DIRECTORY") + "/transport";
    string heavy_path = transport_dir + "/heavy.dat";
    ifstream file(heavy_path.c_str(), ios::in);
    
    if (!file.is_open()) {
        cerr << "Could not open file " + heavy_path + "!\n" << endl;
        exit(1);
    }
    
    string str;    
    file >> str;
    
    CollisionFunc4 func;
    std::map<CollisionPair, int>::iterator iter;
    
    while (str != "STOP" && !collision_map.empty()) {
        // Check if we have landed on a collision identifier
        if (isValidCollisionString(str)) {
            // If so, is this a collision we need to load?
            if ((iter = collision_map.find(str)) != collision_map.end()) {
                // Load Q11, Q22, and B* functions
                file >> func;
                Q11_funcs.push_back(func);
                
                file >> func;
                Q22_funcs.push_back(func);
                
                file >> func;
                Bst_funcs.push_back(func);
                
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
#ifdef VERBOSE
        cout << endl;
        cout << "The following collision pairs were not found!" << endl;
#endif
        //func = CollisionFunc4();
        iter = collision_map.begin();
        for ( ; iter != collision_map.end(); ++iter) {
#ifdef VERBOSE
            cout << "\t" << iter->first.name() << endl;
#endif
            Q11_funcs.push_back(Q11_funcs[0]); // they are not evaluated as
            Q22_funcs.push_back(Q22_funcs[0]); // zero anymore...
            Bst_funcs.push_back(Bst_funcs[0]);
            m_neutral_indices.push_back(iter->second);
        }
#ifdef VERBOSE
        cout << "They will be evaluated as the first integral..." << endl;
        cout << endl;
#endif
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
    
#ifdef USE_COLLISION_INTEGRAL_TABLES
    // Tabulate the collision integral functions using 10 data points
    mp_Q11_table = new LookupTable<double, double, QijTableFunction>(
        std::log(200.0), std::log(20100.0), 10, Q11_funcs.size(), Q11_funcs);//, 0.000001, LINEAR);
    mp_Q22_table = new LookupTable<double, double, QijTableFunction>(
        std::log(200.0), std::log(20100.0), 10, Q22_funcs.size(), Q22_funcs);//, 0.000001, LINEAR);
    mp_Bst_table = new LookupTable<double, double, QRatioTableFunction>(
        std::log(200.0), std::log(20100.0), 10, Bst_funcs.size(), Bst_funcs);//, 0.000001, LINEAR);
    mp_work = new double [Q11_funcs.size()];
    //mp_Q11_table->save("Q11Table.dat");
#endif
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
    const double lnT = std::log(Th);
    
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
#ifdef USE_COLLISION_INTEGRAL_TABLES
            mp_Q11_table->lookup(lnT, mp_work, LINEAR);
            for (int i = 0; i < nn; ++i)
                m_Q11(m_neutral_indices[i]) = mp_work[i];
#else
            for (int i = 0; i < nn; ++i)
                m_Q11(m_neutral_indices[i]) = 1.0E-20 * m_Q11_funcs[i](lnT);
#endif
                
            // Charged collisions
            const double Q11_rep = hfac * sm_Q11_rep(lnTsth);
            for (int i = 0; i < nr; ++i)
                m_Q11(m_repulse_indices[i]) = Q11_rep;
            
            const double Q11_att = hfac * sm_Q11_att(lnTsth);
            for (int i = 0; i < na; ++i)
                m_Q11(m_attract_indices[i]) = Q11_att;
            } break;
        
        case Q12EI: {
            updateCollisionData(Th, Te, nd, p_x, Q11IJ);
            updateCollisionData(Th, Te, nd, p_x, CSTAR);
            for (int i = 0; i < m_ns; ++i)
                m_Q12ei(i) = m_Q11(i)*m_Cst(i);
            } break;
            
        case Q13EI: {
            updateCollisionData(Th, Te, nd, p_x, Q12EI);
            updateCollisionData(Th, Te, nd, p_x, BSTAR);
            for (int i = 0; i < m_ns; ++i)
                m_Q13ei(i) = 1.25*m_Q12ei(i) - 0.25*m_Q11(i)*m_Bst(i);
            } break;
        
        case Q14EI: {
            // Neutral-electron collisions take value of Q(1,3)
            updateCollisionData(Th, Te, nd, p_x, Q13EI);
            m_Q14ei = m_Q13ei;
            
            // Ion-electron (repulsive) interactions
            int index = 0;
            int j = m_repulse_indices[index];
            const double rep = efac * sm_Q14_rep(lnTste);
            while (j < m_ns && index < nr) {
                m_Q14ei(j) = rep;
                j = m_repulse_indices[++index];
            }
            
            // Ion-electron (attractive) interactions
            index = 0;
            j = m_attract_indices[index];
            const double att = efac * sm_Q14_att(lnTste);
            while (j < m_ns && index < na) {
                m_Q14ei(j) = att;
                j = m_attract_indices[++index];
            }
            } break;
            
        case Q15EI: {
            // Neutral-electron collisions take value of Q(1,3)
            updateCollisionData(Th, Te, nd, p_x, Q13EI);
            m_Q15ei = m_Q13ei;
            
            // Ion-electron (repulsive) interactions
            int index = 0;
            int j = m_repulse_indices[index];
            const double rep = efac * sm_Q15_rep(lnTste);
            while (j < m_ns && index < nr) {
                m_Q15ei(j) = rep;
                j = m_repulse_indices[++index];
            }
            
            // Ion-electron (attractive) interactions
            index = 0;
            j = m_attract_indices[index];
            const double att = efac * sm_Q15_att(lnTste);
            while (j < m_ns && index < na) {
                m_Q15ei(j) = att;
                j = m_attract_indices[++index];
            }
            } break;
        
        case Q22IJ: {
#ifdef USE_COLLISION_INTEGRAL_TABLES
            mp_Q22_table->lookup(lnT, mp_work, LINEAR);
            for (int i = 0; i < nn; ++i)
                m_Q22(m_neutral_indices[i]) = mp_work[i];
#else
            for (int i = 0; i < nn; ++i)
                m_Q22(m_neutral_indices[i]) = 1.0E-20 * m_Q22_funcs[i](lnT);
#endif
            
            // Charged collisions
            const double Q22_rep = hfac * sm_Q22_rep(lnTsth);
            for (int i = 0; i < nr; ++i)
                m_Q22(m_repulse_indices[i]) = Q22_rep;
            
            const double Q22_att = hfac * sm_Q22_att(lnTsth);
            for (int i = 0; i < na; ++i)
                m_Q22(m_attract_indices[i]) = Q22_att;
            } break;
        
        case Q23EE: {
            updateCollisionData(Th, Te, nd, p_x, Q22IJ);
            m_Q23ee = m_Q22(0) * sm_Est_rep(lnTste);
            } break;
        
        case Q24EE: {
            m_Q24ee = efac * sm_Q24_rep(lnTste);
            } break;
            
        case ASTAR:
            updateCollisionData(Th, Te, nd, p_x, Q11IJ);
            updateCollisionData(Th, Te, nd, p_x, Q22IJ);
            for (int i = 0; i < m_ncollisions; ++i)
                m_Ast(i) = m_Q22(i) / m_Q11(i);
            break;
                
        case BSTAR: {
#ifdef USE_COLLISION_INTEGRAL_TABLES
            mp_Bst_table->lookup(lnT, mp_work, LINEAR);
            for (int i = 0; i < nn; ++i)
                m_Bst(m_neutral_indices[i]) = mp_work[i];
#else
            // Neutral collisions from the database
            for (int i = 0; i < nn; ++i)
                m_Bst(m_neutral_indices[i]) = m_Bst_funcs[i](lnT);
#endif
            
            // Repulsive collisions
            const double Bst_rep = sm_Bst_rep(lnTsth);
            for (int i = 0; i < nr; ++i)
                m_Bst(m_repulse_indices[i]) = Bst_rep;
            
            // Attractive collisions
            const double Bst_att = sm_Bst_att(lnTsth);
            for (int i = 0; i < na; ++i)
                m_Bst(m_attract_indices[i]) = Bst_att;
            } break;
            
        case CSTAR: {        
            // Neutral-electron collisions are treated as 1 for now
            m_Cst = 1.0;
            
            // Ion-electron (repulsive) interactions
            int index = 0;
            int j = m_repulse_indices[index];
            const double Cst_rep = sm_Cst_rep(lnTste);
            while (j < m_ns && index < nr) {
                m_Cst(j) = Cst_rep;
                j = m_repulse_indices[++index];
            }
            
            // Ion-electron (attractive) interactions
            index = 0;
            j = m_attract_indices[index];
            const double Cst_att = sm_Cst_att(lnTste);
            while (j < m_ns && index < na) {
                m_Cst(j) = Cst_att;
                j = m_attract_indices[++index];
            }
            } break;
        
        case ETAI: {
            updateCollisionData(Th, Te, nd, p_x, Q22IJ);
            
            // Electron does not contribute to viscosity
            int k = nSpecies() - nHeavy();
            m_eta(0) = 0.0;
            
            // Heavy species
            for (int i = k; i < m_ns; ++i)
                m_eta(i) = 
                    5.0 / 16.0 * sqrt(PI * KB * Th * m_mass(i)) / m_Q22(i,i);
            } break;
            
        case NDIJ: {
            updateCollisionData(Th, Te, nd, p_x, Q11IJ);
            
            int k = nSpecies() - nHeavy();
            
            if (k == 1) {
                // Electron-Electron
                m_Dij(0) = 3.0 / 8.0 * sqrt(PI * KB * Te / m_mass(0)) / m_Q11(0);
                
                // Electron-heavy interactions
                for (int i = 1; i < m_ns; ++i)
                    m_Dij(i) = 
                        3.0 / 16.0 * sqrt(TWOPI * KB * Te / m_mass(0)) / m_Q11(i);
            }
            
            // Heavy-heavy interactions
            for (int i = k * m_ns; i < m_ncollisions; ++i)
                m_Dij(i) = 3.0 / 16.0 * sqrt(TWOPI * KB * Th / 
                    m_red_mass(i)) / m_Q11(i);
            } break;
        
        default:
            break;  // Will never get here
    }
    
    mp_last_T[data] = Th;
}

    } // namespace Transport
} // namespace Mutation


