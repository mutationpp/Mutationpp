
#include "MultiPhaseEquilSolver.h"
#include <set>

using std::cout;
using std::endl;

using namespace Mutation::Numerics;

namespace Mutation {
    namespace Thermodynamics {

//==============================================================================

MultiPhaseEquilSolver::MultiPhaseEquilSolver(const Thermodynamics& thermo)
    : m_thermo(thermo), m_phase(m_thermo.nSpecies())
{ 
    m_ns = m_thermo.nSpecies();
    m_ne = m_thermo.nElements();
    m_B  = m_thermo.elementMatrix();
    
    initPhases();
    
    cout << "# species:  " << m_ns << endl;
    cout << "# elements: " << m_ne << endl;
    cout << "# phases:   " << m_np << endl; 
    cout << "phases" << endl;
    
    for (int i = 0; i < m_ns; ++i)
        cout << m_phase(i) << endl;
        
    cout << "B" << endl;
    cout << m_B << endl;
}

//==============================================================================

void MultiPhaseEquilSolver::initPhases()
{
    std::set<int> phases;
    
    // Assign an integer value to each species corresponding to the phase it
    // belongs to (also keep track of unique phases)
    for (int i = 0; i < m_ns; ++i) {
        m_phase(i) = (int) m_thermo.species(i).phase();
        phases.insert(m_phase(i));
    }
    
    // Number of phases = number of unique phase numbers
    m_np = phases.size();
    
    // Make sure all the phase numbers are in the set {0,...,np-1}
    bool empty;
    for (int p = 0; p < m_np; ++p) {
        empty = true;
        for (int i = 0; i < m_ns; ++i) {
            if (m_phase(i) == p) {
                empty = false;
                break;
            }
        }
        
        if (empty) {
            for (int i = 0; i < m_ns; ++i)
                if (m_phase(i) > p) m_phase(i)--;
            p--;
        }
    }
}

//==============================================================================

std::pair<int,int> MultiPhaseEquilSolver::equilibrate(
    double T, double P, const double *const p_ev, double *const p_sv)
{
    static RealVector gtp(m_ns);
    
    // Compute the unitless Gibbs function for each species
    m_thermo.speciesGOverRT(T, P, p_sv);
    gtp = asVector(p_sv, m_ns);
    
    cout << "G/RT(T, P)" << endl;
    cout << gtp << endl;
}

//==============================================================================
/*
void MultiPhaseEquilSolver::maxmin(const RealVector& c, RealVector& nmm) const
{
    // Setup the tableau for input to lp()
    RealMatrix tableau(m_nrc, m_nsu+1);
    
    // B'
    tableau(0,m_nrc,0,m_nsu) = m_B_rt;    
    
    // sum(B)'
    for (size_t i = 0; i < m_nrc; ++i)
        for (size_t j = 0; j < m_nsu; ++j)
            tableau(i,m_nsu) += tableau(i,j);
    
    // Objective function
    RealVector f(m_nsu+1);
    f(m_nsu) = 1.0;
    
    // Temporary solution vector (to make room for zmin)
    RealVector x(m_nsu+1);
    
    // Now solve the LP problem
    double zmin;
    LpResult lp_result = lp(f, MAXIMIZE, tableau, c, 0, 0, x, zmin);
    
    if (lp_result != SOLUTION_FOUND) {
        cout << "Error finding max-min solution!" << endl;
        exit(1);
    }
    
    // Transfer solution
    nmm = x(0,m_nsu) + zmin;
}
*/
//==============================================================================
/*
void MultiPhaseEquilSolver::ming(
    const double& T, const double& P, const RealVector& cr, RealVector& gu,
    RealVector& ng) const
{
    double g[m_ns];
    
    // Compute the standard state normalized molar specific Gibbs energy for
    // each species in the mechanism
    m_thermo.speciesGOverRT(T, m_thermo.standardStateP(), g);
    
    // Store the standard state g's for the undetermined species and update to
    // account for pressure
    double logP = std::log(P / m_thermo.standardStateP());
    for (int i = 0; i < m_nsu; ++i)
        gu(i) = g[m_species_order[m_nsd + i]] + logP;
    
    // Solve the linear programming problem
    LpResult lp_result =
        lp(gu, MINIMIZE, m_B_rt, cr, 0, 0, ng, g[0]);
    
    if (lp_result != SOLUTION_FOUND) {
        cout << "Error finding min-g solution!" << endl;
        exit(1);
    }
}
*/
//==============================================================================

    } // namespace Thermodynamics
} // namespace Mutation
