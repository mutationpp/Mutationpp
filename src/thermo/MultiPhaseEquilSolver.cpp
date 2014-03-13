
#include "MultiPhaseEquilSolver.h"
#include "minres.h"
#include "gmres.h"

#include <set>
#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>

//#define VERBOSE

using std::cout;
using std::endl;
using std::setw;

using namespace Mutation::Numerics;

namespace Mutation {
    namespace Thermodynamics {
    

//==============================================================================

MultiPhaseEquilSolver::MultiPhaseEquilSolver(
    const Thermodynamics& thermo, const bool pure_condensed) :
        m_thermo(thermo),
        m_pure_condensed(pure_condensed), m_solution(thermo)
{
    // Sizing information
    m_ns  = m_thermo.nSpecies();
    m_ne  = m_thermo.nElements();
    m_nc  = m_ne;
    
    // System element matrix
    m_B  = m_thermo.elementMatrix();
    
    m_tableau_capacity = (m_nc+2)*(m_ns+2);
    mp_tableau = new double [m_tableau_capacity];
    mp_ming    = new double [m_ns];
    mp_maxmin  = new double [m_ns];
    mp_g       = new double [m_ns];
    mp_g0      = new double [m_ns];
    mp_c       = new double [m_nc];
    
    // Compute phase information
    mp_phase = new int [m_ns];
    initPhases();
    
    // Initialize storage for the solution object
    m_solution.initialize(m_np, m_nc, m_ns);
}

//==============================================================================

MultiPhaseEquilSolver::~MultiPhaseEquilSolver()
{
    delete mp_phase;
    delete mp_tableau;
    delete mp_ming;
    delete mp_maxmin;
    delete mp_g;
    delete mp_g0;
    delete mp_c;
};

//==============================================================================

void MultiPhaseEquilSolver::addConstraint(const double *const p_A)
{
//    // Save constraints
//    m_constraints.push_back(asVector(p_A, m_ns));
//    
//    // Update the B matrix
//    m_B = zeros<double>(m_ns, ++m_nc);
//    m_B(0,m_ns,0,m_ne) = m_thermo.elementMatrix();
//    
//    for (int i = 0; i < m_nc-m_ne; ++i)
//        m_B.col(m_ne+i) = m_constraints[i];
}

//==============================================================================
        
void MultiPhaseEquilSolver::clearConstraints()
{
//    m_constraints.clear();
//    m_B = m_thermo.elementMatrix();
//    m_nc = m_ne;
}

//==============================================================================
        
void MultiPhaseEquilSolver::dXdT(double *const p_dxdt) const
{
//    // Compute the dg/dT = [-H/RT]/T values
//    RealVector dgdt(m_nsr);
//    m_thermo.speciesHOverRT(m_T, p_dxdt);
//    for (int j = 0; j < m_nsr; ++j)
//        dgdt(j) = -p_dxdt[mp_sjr[j]] / m_T;
//    
//    partialOfX(dgdt, p_dxdt);
}

//==============================================================================

void MultiPhaseEquilSolver::dXdP(double *const p_dxdp) const
{
//    // Compute the dg/dP = 1/P for gas species
//    RealVector dgdp(m_nsr, 1.0 / m_P);
//    
//    // Zero for all other species
//    for (int j = 0; j < m_nsr; ++j)
//        if (mp_phase[mp_sjr[j]] > 0) dgdp(j) = 0.0;
//    
//    partialOfX(dgdp, p_dxdp);
}

//==============================================================================

void MultiPhaseEquilSolver::partialOfX(
    const RealVector& dg, double* const p_dx) const
{
//    // RHS
//    double temp;
//    RealVector rhs(m_ncr+m_np, 0.0);
//    RealVector sol = rhs;
//    
//    for (int k = 0; k < m_nsr; k++) {
//        temp = m_N(k)*dg(k);
//        for (int i = 0; i < m_ncr; ++i)
//            rhs(i) += temp * m_Br(k,i);
//        rhs(m_ncr+mp_phase[mp_sjr[k]]) += temp;
//    }
//    
//    // Now compute the system solution
//    RealSymMat A(m_ncr+m_np);
//    formSystemMatrix(A);
//    
//    // Solve the system for dlambda/dalpha
//    LeastSquares<double> ls(A);
//    ls.solve(sol, rhs);
//    
//    // Compute dXj/dalpha = Xj*[-dgj/dalpha + sum_i Bji dlambdai/dalpha]
//    RealVector Nbar(m_np, 0.0);
//    for (int j = 0; j < m_nsr; ++j)
//        Nbar(mp_phase[mp_sjr[j]]) += m_N(j);
//    
//    std::fill(p_dx, p_dx+m_ns, 0.0);
//    
//    for (int j = 0; j < m_nsr; ++j) {
//        p_dx[mp_sjr[j]] = -dg(j);
//        for (int i = 0; i < m_ncr; ++i)
//            p_dx[mp_sjr[j]] += m_Br(j,i) * sol(i);
//        p_dx[mp_sjr[j]] *= m_N(j) / Nbar(mp_phase[mp_sjr[j]]);
//    }
}


//==============================================================================

void MultiPhaseEquilSolver::initPhases()
{
    std::set<int> phases;
    
    // Assign an integer value to each species corresponding to the phase it
    // belongs to (also keep track of unique phases)
    int npure = 0;
    for (int i = 0; i < m_ns; ++i) {
        mp_phase[i] = (int) m_thermo.species(i).phase();
        if (mp_phase[i] > 0 && m_pure_condensed)
            mp_phase[i] = ++npure;
        phases.insert(mp_phase[i]);
    }
    
    // Number of phases = number of unique phase numbers
    m_np = phases.size();
    
    // Make sure all the phase numbers are in the set {0,...,np-1}
    bool empty;
    for (int p = 0; p < m_np; ++p) {
        empty = true;
        for (int i = 0; i < m_ns; ++i) {
            if (mp_phase[i] == p) {
                empty = false;
                break;
            }
        }
        
        if (empty) {
            for (int i = 0; i < m_ns; ++i)
                if (mp_phase[i] > p) mp_phase[i]--;
            p--;
        }
    }
}

//==============================================================================

std::pair<int, int> MultiPhaseEquilSolver::equilibrate(
    double T, double P, const double* const p_cv, double* const p_sv)
{
    // Save T, P, c
    m_T  = T;
    m_P  = P;
    std::copy(p_cv, p_cv+m_nc, mp_c);
    #ifdef VERBOSE
    std::cout << "equilibrate(" << T << " K, " << P << " Pa,";
    for (int i = 0; i < m_ne; ++i)
        std::cout << m_thermo.elementName(i) << " " << p_cv[i] << (i == m_ne-1 ? ")" : ",");
    std::cout << std::endl;
    #endif
    
    // Compute species gibbs free energies at the given temperature and pressure
    m_thermo.speciesGOverRT(T, P, mp_g);
    
    // Setup species ordering and remove "determined" species from consideration
    // (reduced problem)
    checkForDeterminedSpecies();
    
    // Compute the initial conditions lambda(0), Nbar(0), N(0), and g(0)
    initialConditions();
    
    #ifdef VERBOSE
//    cout << "Gj/RT = " << endl;
//    for (int j = 0; j < m_ns; ++j)
//        cout << mp_g[j] << endl;
//    m_solution.printOrder();
    cout << "Initial conditions: " << endl;
    m_solution.printSolution();
//    m_solution.printG();
    #endif

    double s  = 0.0;
    double ds = 1.0;
    Solution last_solution(m_thermo);
    
    // Integrate the equilibrium solution from s = 0 to s = 1
    m_niters = 0;
    m_nnewts = 0;
    while (s < 1.0) {
        // Save current solution incase we need to take a smaller step
        last_solution = m_solution;
        
        // Compute the rates for lambda and Nbar
        RealVector dx(m_solution.ncr()+m_solution.npr());
        rates(dx);
        #ifdef VERBOSE
        cout << "iter = " << m_niters << ", s = " << s << ", ds = " << ds << endl;
        cout << "dx = " << endl;
        cout << dx << endl;
        #endif

        // Keep trying to obtain a solution at progressively smaller steps
        // until a solution is obtained
        bool solution_obtained = false;
        while (!solution_obtained) {
            // Get trial solution
            m_solution.setG(mp_g0, mp_g, s+ds);
            m_solution.update(dx*ds, m_B);
            //m_solution.printSolution();

            // Use newton to reduce the residual
            solution_obtained = newton();
            //while (solution_obtained && phaseRedistribution())
            //    solution_obtained = newton();
            
            // If newton fails to converge then reduce the step size and try
            // again
            if (!solution_obtained) {
                #ifdef VERBOSE
                cout << "Could not converge to solution, reducing step size" << endl;
                #endif
                ds *= 0.25;
                m_solution = last_solution;
            }
            
            if (m_niters + m_nnewts > 1000) {
                cout << "Touble converging..." << endl;
                exit(1);
            }
        }
        
        // We found the solution at s, increase step size and move forward
        #ifdef VERBOSE
        cout << "Step succeded, increasing step size" << endl;
        #endif
        s += ds;
        ds = std::min(ds*4.0, 1.0-s);
        m_niters++;
    }
    
    // Finally, unwrap the solution for the user and return convergence stats
    for (int j = 0; j < m_ns; ++j)
        p_sv[j] = 0.0;
    double sum = 0.0;
    for (int j = 0; j < m_solution.nsr(); ++j) {
        p_sv[m_solution.sjr()[j]] = m_solution.y()[j]*m_solution.y()[j];
        sum += p_sv[m_solution.sjr()[j]];
    }
    for (int j = 0; j < m_ns; ++j)
        p_sv[j] /= sum;
    
    return std::make_pair(m_niters, m_nnewts);
}


//std::pair<int,int> MultiPhaseEquilSolver::equilibrate(
//    double T, double P, const double *const p_cv, double *const p_sv)
//{
//    m_T = T;
//    m_P = P;
//
//#ifdef VERBOSE
//    std::cout << "equilibrate(" << T << " K, " << P << " Pa,";
//    for (int i = 0; i < m_ne; ++i)
//        std::cout << m_thermo.elementName(i) << " " << p_cv[i] << (i == m_ne-1 ? ")" : ",");
//    std::cout << std::endl;
//#endif
//
//    static RealVector g0;
//    
//    // Begin by reducing zero species
//    m_c = asVector(p_cv, m_nc);
//    initialize();
//
//    m_niters = 0;
//    m_nnewts = 0;
//    
//    bool still_solving = true;
//    while (still_solving) {
//    
//
//    m_Nbar = ones<double>(m_npr);
//    RealVector Nbar_trial;
//    
//    //RealVector lambda(m_ncr);
//    m_lambda = ones<double>(m_ncr);
//    RealVector lambda_trial;
//    
//    // Resize vectors if need be
//    m_dg = ones<double>(m_nsr);
//    
//    // Compute the unitless Gibbs function for each nonzero species
//    // (temporarily store in dg)
//    updateMinGSolution();
//    for (int i = 0; i < m_nsr; ++i)
//        m_dg(i) = mp_g[mp_sjr[i]];
//
//    // Compute the initial conditions for the integration
//    initialConditions(m_dg, m_lambda, m_Nbar, g0);
//    m_dg -= g0;
//    m_g  =  g0;
//    //updatePreConditioner();
//    
//    RealSymMat A(m_ncr+m_npr);
//    RealVector r(m_ncr+m_npr);
//    RealVector dx(m_ncr+m_npr);
//    
//    double s  = 0.0;
//    double ds = 1.0;
//    double error_tol = 1.0e-6;
//    
//    bool update_dx = true;
//    
//    //updateSpeciesMoles(m_lambda, m_Nbar, g0);
//    
//    // Integrate quantities from s = 0 to s = 1
//    while (s < 1.0) {
//        m_niters++;
//        
//#ifdef VERBOSE
//        std::cout << "Step: iter = " << m_niters << ", s = " << s << ", ds = "
//                  << ds << ", newt = " << m_nnewts << endl;
//#endif
//        
//        if (update_dx) {
//            m_g = g0 + s*m_dg;
//            rates(dx);
//            //cout << "dx from continuation step" << endl;
//            //cout << dx << endl;
//        }
//
//        
//        // Update trial values
//        lambda_trial = m_lambda + dx(0,m_ncr)*ds;
//        Nbar_trial   = m_Nbar * exp(dx(m_ncr,m_ncr+m_npr)*ds);
//        m_g          = g0 + m_dg*(s+ds);
//
//        // Compute the 
//        if (!updateSpeciesMoles(lambda_trial, Nbar_trial, m_g)) {
//            ds *= 0.5;
//            update_dx = false;
//            continue;
//        }
//        
//        computeResidual(Nbar_trial, r);
//        
//        std::pair<int,double> res = 
//            newton(lambda_trial, Nbar_trial, m_g, r, A, 1.0e-6, 5);
//        m_nnewts += res.first;
//        
//        if ((m_niters + m_nnewts) > 1000) {
//            cout << "Trouble converging...." << endl;
//            cout << "equilibrate(" << T << " K, " << P << " Pa,";
//            for (int i = 0; i < m_ne; ++i)
//                cout << m_thermo.elementName(i) << " " << p_cv[i]
//                     << (i == m_ne-1 ? ")" : ",");
//            cout << endl;
//
//            cout << "Step: iter = " << m_niters << ", s = " << s << ", ds = " << ds
//                 << ", newt = " << m_nnewts << endl;
//            exit(1);
//        }
//        
//        update_dx = false;
//        if (m_nsr != g0.size()) {
//            dx = RealVector(m_ncr+m_npr);
//            g0 = m_g - m_dg*(s+ds);
//            update_dx = true;
//        }
//        
//        // If newton did not converge, reduce the step size and start over unless
//        // the stepsize was already very small, in that case just take the step
//        if (res.second > error_tol) {
//#ifdef VERBOSE
//            std::cout << "Newton did not converge, reducing step size" << std::endl;
//#endif
//            ds *= 0.25;
//            continue;
//        }
//        
//        error_tol = res.second + std::max(1.0e-6, 0.05*res.second);
//        
//        // Update the system variables
//        m_lambda = lambda_trial;
//        m_Nbar   = Nbar_trial;
//        s += ds;
//        //ds *= 2.0;
//        
//        ds = std::min(ds*4.0, 1.0-s);
//        
//        if (ds == 0.0 && s < 1.0) {
//            std::cerr << "Continuation step size dropped to zero in equil solver!" << endl;
//            std::exit(1);
//        }
//        
//        update_dx = true;
//    }
//    
//#ifdef VERBOSE
//    cout << "s = 1 achieved, reducing residual..." << endl;
//#endif
//    
//    // Use Newton iterations to improve residual for final answer
//    std::pair<int,double> res = newton(m_lambda, m_Nbar, m_g, r, A, 1.0e-8, 20);
//    m_nnewts += res.first;
//    
//    still_solving = phaseRedistribution();
//    
//    } // still solving
//    
//    // Compute the species mole fractions
//    for (int j = 0; j < m_nsr; ++j)
//        p_sv[mp_sjr[j]] = m_N(j);
//    for (int j = m_nsr; j < m_ns; ++j)
//        p_sv[mp_sjr[j]] = 0.0;
//    
//    m_Nbar = RealVector(m_np, 0.0);
//    for (int j = 0; j < m_ns; ++j)
//        m_Nbar(mp_phase[j]) += p_sv[j];
//
//    double total_moles = m_Nbar.sum();
//    for (int j = 0; j < m_ns; ++j) {
//        p_sv[j] /= total_moles;
//    }
//    
////    double pmole;
////    for (int j = 0; j < m_ns; ++j) {
////        pmole = m_Nbar(mp_phase[j]);
////        if (pmole > 0.0)
////            p_sv[j] /= pmole;
////    }
//    
//    //exit(1);
//#ifdef VERBOSE
//    cout << "Solution obtained :) - iters = " << m_niters << ", newts = "
//         << m_nnewts << endl;
//#endif
//    return std::make_pair(m_niters, m_nnewts);
//}

//void MultiPhaseEquilSolver::removePhase(const int phase)
//{
//#ifdef VERBOSE
//    cout << "Removing phase " << phase << " ( ";
//    for (int j = mp_sizes[phase]; j < mp_sizes[phase+1]; ++j)
//        cout << m_thermo.speciesName(mp_sjr[j]) << " ";
//    cout << ")" << endl;
//#endif
//    
//    // Check that the phase number is feasible
//    assert(phase >= 0);
//    assert(phase < m_npr);
//    assert(m_npr > 1);
//    
//    // If the phase is not the last non-empty phase, we need to shift it to the
//    // right to make the non-empty species list contiguous
//    const int size = mp_sizes[phase+1] - mp_sizes[phase];
//    if (phase != m_npr-1) {
//        int temp [size];
//        
//        // First copy the phase to be removed into temporary array
//        int* p = temp;
//        for (int i = mp_sizes[phase]; i < mp_sizes[phase+1]; ++i)
//            *p++ = mp_sjr[i];
//        
//        // Shift all remaining phases to fill the space
//        p = mp_sjr + mp_sizes[phase];
//        for (int i = mp_sizes[phase+1]; i < mp_sizes[m_np]; ++i)
//            *p++ = mp_sjr[i];
//        
//        // Place the removed phase at the end of the list (before determined
//        // species)
//        for (int i = phase+1; i < m_np; ++i)
//            mp_sizes[i] = mp_sizes[i+1] - size;
//        p = temp;
//        for (int i = mp_sizes[m_np-1]; i < mp_sizes[m_np]; ++i)
//            mp_sjr[i] = *p++;
//    }
//    
//    // Update reduced species phase
//    m_npr--;
//    m_nsr -= size;
//#ifdef VERBOSE
//    cout << "New phase ordering:" << endl;
//    cout << "(Included)" << endl;
//    for (int m = 0; m < m_npr; ++m) {
//        cout << m << ": ";
//        for (int j = mp_sizes[m]; j < mp_sizes[m+1]; ++j)
//            cout << m_thermo.speciesName(mp_sjr[j]) << " ";
//        cout << endl;
//    }
//    cout << "(Excluded)" << endl;
//    for (int m = m_npr; m < m_np; ++m) {
//        cout << m << ": ";
//        for (int j = mp_sizes[m]; j < mp_sizes[m+1]; ++j)
//            cout << m_thermo.speciesName(mp_sjr[j]) << " ";
//        cout << endl;
//    }
//#endif
//}

//void MultiPhaseEquilSolver::addPhase(const int phase)
//{
//#ifdef VERBOSE
//    cout << "Adding phase " << phase << " ( ";
//    for (int j = mp_sizes[phase]; j < mp_sizes[phase+1]; ++j)
//        cout << m_thermo.speciesName(mp_sjr[j]) << " ";
//    cout << ")" << endl;
//#endif
//
//    assert(phase >= m_npr);
//    assert(phase < m_np);
//    
//    int size = mp_sizes[phase+1] - mp_sizes[phase];
//    
//    // If the phase is not the first empty phase, then we need to shift it to
//    // the front
//    if (phase > m_npr) {
//        int temp [size];
//        
//        // First copy the phase to be added into temporary array
//        int* p = temp;
//        for (int i = mp_sizes[phase]; i < mp_sizes[phase+1]; ++i)
//            *p++ = mp_sjr[i];
//        
//        // Shift phases that come before to the back
//        p = mp_sjr + mp_sizes[phase+1]-1;
//        for (int i = mp_sizes[phase]-1; i >= m_nsr; --i)
//            *p-- = mp_sjr[i];
//        
//        // Place the added phase at the end of the included list
//        for (int i = phase+1; i > m_npr; --i)
//            mp_sizes[i] = mp_sizes[i-1] + size;
//        p = temp;
//        for (int i = mp_sizes[m_npr]; i < mp_sizes[m_npr+1]; ++i)
//            mp_sjr[i] = *p++;
//    }
//    
//    m_npr++;
//    m_nsr += size;
//#ifdef VERBOSE
//    cout << "New phase ordering:" << endl;
//    cout << "(Included)" << endl;
//    for (int m = 0; m < m_npr; ++m) {
//        cout << m << ": ";
//        for (int j = mp_sizes[m]; j < mp_sizes[m+1]; ++j)
//            cout << m_thermo.speciesName(mp_sjr[j]) << " ";
//        cout << endl;
//    }
//    cout << "(Excluded)" << endl;
//    for (int m = m_npr; m < m_np; ++m) {
//        cout << m << ": ";
//        for (int j = mp_sizes[m]; j < mp_sizes[m+1]; ++j)
//            cout << m_thermo.speciesName(mp_sjr[j]) << " ";
//        cout << endl;
//    }
//#endif
//}

bool MultiPhaseEquilSolver::phaseRedistribution()
{
    #ifdef VERBOSE
    cout << "Phase redistribution:" << endl;
    #endif
    
    const double phase_tol = std::log(1.0e-20);

    const int* const p_sizes = m_solution.sizes();
    const int* const p_sjr   = m_solution.sjr();
    const int* const p_cir   = m_solution.cir();
    
    // Remove any phase which has negligable moles
    int m;
//    int m = 1;
//    while (m < m_solution.npr()) {
//        if (m_solution.lnNbar()[m] < phase_tol)
//            m_solution.removePhase(m);
//        else
//            m++;
//    }
    
    // Check to see if there are any condensed phases that could reduce the
    // residual further as the phase rule allows
    int j = p_sizes[m_solution.npr()], jk, lowest_phase = m_solution.npr()-1;
    double test, lowest_value = 1.0;
    
    for (m = m_solution.npr(); m < m_np; ++m) {
        for ( ; j < p_sizes[m+1]; ++j) {
            jk = p_sjr[j];
            test = m_solution.g()[jk];
            for (int i = 0; i < m_solution.ncr(); ++i)
                test -= m_solution.lambda()[i]*m_B(jk,p_cir[i]);
            
            if (test < lowest_value) {
                lowest_value = test;
                lowest_phase = m;
            }
        }
    }
    
    bool phase_redist = (lowest_value < 0.0);
    
    // Add the phase with the most negative vapor pressure test value
    if (phase_redist) {
        int phase = m_solution.addPhase(lowest_phase);
        double lambda [m_solution.ncr()];
        double lnNbar [m_solution.npr()];
        
        std::copy(m_solution.lambda(), m_solution.lambda()+m_solution.ncr(), lambda);
        std::copy(m_solution.lnNbar(), m_solution.lnNbar()+m_solution.npr(), lnNbar);
        lnNbar[phase] = std::log(1.0e-4);
        
        for (int i = 0; i < m_solution.npr(); ++i)
            lnNbar[i] = std::exp(lnNbar[i]);
        m_solution.setSolution(lambda, lnNbar, m_B);

        // Adhere to the phase rule (npr <= ncr)
        if (m_solution.npr() > m_solution.ncr()) {
            lowest_phase = 1;
            lowest_value = m_solution.lnNbar()[1];
            
            for (int m = 2; m < m_solution.npr()-1; ++m) {
                if (m_solution.lnNbar()[m] < lowest_value) {
                    lowest_value = m_solution.lnNbar()[m];
                    lowest_phase = m;
                }
            }
            
            m_solution.removePhase(lowest_phase);
        }
    }
    
    return phase_redist;
}



void MultiPhaseEquilSolver::rates(RealVector& dx)
{
    using std::sqrt;
    using std::exp;
    
    // Get some of the solution variables
    const int npr = m_solution.npr();
    const int ncr = m_solution.ncr();
    const int nsr = m_solution.nsr();
    
    const int* const p_sjr   = m_solution.sjr();
    const int* const p_cir   = m_solution.cir();
    const int* const p_sizes = m_solution.sizes();
    
    const double* const p_y = m_solution.y();

    // Compute a least squares factorization of H
    HMatrix H = m_solution.hMatrix(m_B);
    LeastSquares<double> ls(H);
    
    // Use tableau for temporary storage
    double* p_rhs   = mp_tableau;
    double* p_dlamg = p_rhs + nsr;
    double* p_dlamy = p_dlamg + ncr;

    // Compute dlamg
    int jk;
    for (int j = 0; j < nsr; ++j) {
        jk = p_sjr[j];
        p_rhs[j] = p_y[j] * (mp_g[jk] - mp_g0[jk]);
    }
    ls.solve(p_dlamg, p_rhs);
    
    // Compute dlambda_m for each phase m
    int j;
    for (int m = 0; m < npr; ++m) {
        for (j = 0 ; j < p_sizes[m]; ++j)
            p_rhs[j] = 0.0;
        for ( ; j < p_sizes[m+1]; ++j)
            p_rhs[j] = p_y[j];
        for ( ; j < nsr; ++j)
            p_rhs[j] = 0.0;
        ls.solve(p_dlamy+m*ncr, p_rhs);
    }
    
    // Compute the linear system to be solved for the d(lnNbar)/ds variables
    RealMatrix A(npr, npr);
    
    double sum;
    for (int p = 0; p < npr; ++p) {
        for (int m = 0; m < npr; ++m) {
            for (int j = p_sizes[m]; j < p_sizes[m+1]; ++j) {
                sum = 0.0;
                for (int i = 0; i < ncr; ++i)
                    sum += H(j,i)*p_dlamy[p*ncr+i];
                A(m,p) += p_y[j]*sum;
            }
        }
    }
    
    RealVector b(npr);
    for (int m = 0; m < npr; ++m) {
        b(m) = 0.0;
        for (int j = p_sizes[m]; j < p_sizes[m+1]; ++j) {
            sum = 0.0;
            for (int i = 0; i < ncr; ++i)
                sum += H(j,i)*p_dlamg[i];
            jk = p_sjr[j];
            b(m) += p_y[j]*(sum - p_y[j]*(mp_g[jk]-mp_g0[jk]));
        }
    }
    
    #ifdef VERBOSE
    cout << "A in rates: " << endl;
    cout << A << endl;
    cout << "b in rates: " << endl;
    cout << b << endl;
    #endif
    
    // Solve the linear system to get d(lnNbar)/ds
    switch (npr) {
        case 1:
            dx(ncr) = b(0)/A(0);
            break;
        case 2: {
            double det  = A(0)*A(3) - A(1)*A(2);
            dx(ncr)   = (A(3)*b(0)-A(1)*b(1))/det;
            dx(ncr+1) = (A(0)*b(1)-A(2)*b(0))/det;
            break;
        } case 3: {
            double c11 = A(4)*A(8)-A(5)*A(7);
            double c12 = A(3)*A(8)-A(5)*A(6);
            double c13 = A(3)*A(7)-A(4)*A(6);
            double det = A(0)*c11 - A(1)*c12 + A(2)*c13;
            dx(ncr) = (
                b(0)*c11-
                A(1)*(b(1)*A(8)-A(5)*b(2))+
                A(2)*(b(1)*A(7)-A(4)*b(2)))/det;
            dx(ncr+1) = (
                A(0)*(b(1)*A(8)-A(5)*b(2))-
                b(0)*c12+
                A(2)*(A(3)*b(2)-b(1)*A(6)))/det;
            dx(ncr+2) = (
                A(0)*(A(4)*b(2)-b(1)*A(5))-
                A(1)*(A(3)*b(2)-b(1)*A(6))+
                b(0)*c13)/det;
            break;
        } default:
            RealVector x(npr);
            //QRP<double> qrp(A);
            //qrp.solve(x, b);
            SVD<double> svd(A);
            svd.solve(x, b);
            dx(ncr,ncr+npr) = x;
            
        }
    
    // Finally compute the d(lambda)/ds vector
    for (int i = 0; i < ncr; ++i)
        dx(i) = p_dlamg[i];
    for (int m = 0; m < npr; ++m)
        for (int i = 0; i < ncr; ++i)
            dx(i) -= dx(ncr+m)*p_dlamy[m*ncr+i];
}

//==============================================================================

bool MultiPhaseEquilSolver::newton()
{
    #ifdef VERBOSE
    cout << "newton:" << endl;
    #endif
    const double tolerance = 1.0e-6;
    const int    max_iters = 5;
    const double phase_tol = std::sqrt(1.0e-6);
    
    int npr = m_solution.npr();
    const int ncr = m_solution.ncr();
    const int nsr = m_solution.nsr();
    const int* const p_sjr = m_solution.sjr();
    const int* const p_cir = m_solution.cir();
    const int* const p_sizes = m_solution.sizes();
    const double* const p_lnNbar = m_solution.lnNbar();
    
    // First compute the residual
    RealSymMat A(ncr+npr);
    RealVector r(ncr+npr);
    RealVector dx(ncr+npr);
    computeResidual(r);
    
    double res = r.norm2();
    
    LDLT<double> ldlt;
    
    if (res > 1.0)
        return false;
    
    int iter = 0;
    while (res > tolerance && iter < max_iters) {
        #ifdef VERBOSE
        cout << "newton, iter = " << iter << ", " << "res = " << res << endl;
        m_solution.printSolution();
        #endif
        // First check if a phase needs to be removed (always include gas phase)
        int m = 1;
        while(m < m_solution.npr()) {
            double sum = 0.0;
            for (int j = p_sizes[m]; j < p_sizes[m+1]; ++j)
                sum += m_solution.y()[j];
            if (sum < phase_tol)
                m_solution.removePhase(m);
            else
                m++;
        }
        
        // Update residual if we need to
        if (npr != m_solution.npr()) {
            npr = m_solution.npr();
            r.resize(ncr+npr);
            computeResidual(r);
            res = r.norm2();
            A = RealSymMat(ncr+npr);
            dx = RealVector(ncr+npr);
        }
        
        // Compute the system jacobian
        formSystemMatrix(A);
        #ifdef VERBOSE
        cout << "jacobian matrix = " << endl;
        cout << A << endl;
        cout << "residual vector = " << endl;
        cout << r << endl;
        #endif
        
        // Solve the linear system
        //ldlt.setMatrix(A);
        //ldlt.solve(dx, -r);
        QRP<double> qrp(A);
        #ifdef VERBOSE
        cout << "R matrix = " << endl;
        cout << qrp.R() << endl;
        #endif
        r = -r;
        qrp.solve(dx, r);
        
        #ifdef VERBOSE
        cout << "dx = " << endl;
        cout << dx << endl;
        #endif
        
        // Update the solution
        m_solution.update(dx, m_B);
        
        // Update residual
        computeResidual(r);
        double new_res = r.norm2();
        #ifdef VERBOSE
        cout << "new res = " << new_res << endl;
        #endif
        iter++;
        if (new_res > res)
            break;
        else
            res = new_res;
    }
    
    m_nnewts += iter;
    return (res <= tolerance);
}

////==============================================================================
//
//std::pair<int,double> MultiPhaseEquilSolver::newton(
//    RealVector& lambda, RealVector& Nbar, RealVector& g,
//    RealVector& r, RealSymMat& A, const double tol, const int max_iters)
//{
//    using std::exp;
//    using std::log;
//    
//    const double alpha = 1.0e-6;
//    const double phase_tol = 1.0e-30;
//    int n = m_ncr+m_npr;
//    
//    double res0, res, lam;
//    int iters = 0;
//    
//    static LDLT<double> ldlt;
//    RealVector xt(n);
//    RealVector dx(n);
//    RealVector x(n);
//    
//    x(0,m_ncr) = lambda;
//    x(m_ncr,n) = log(Nbar);
//    
//    res0 = r.norm2();
//
//#ifdef VERBOSE
//    cout << "Newton: res = " << res0 << endl;
//#endif
//    while (res0 > tol && iters++ < max_iters) {
//        // Check if we need to remove any phases before continuing
//        int m = 0;
//        while (m < m_npr) {
//            if (Nbar(m) < phase_tol) {
//#ifdef VERBOSE
//                cout << "Removing phase " << m << "..." << endl;
//                cout << "Nbar(m) = " << Nbar(m) << endl;
//#endif
//                // Update Nbar
//                RealVector temp(m_npr-1);
//                for (int p = 0; p < m; ++p)
//                    temp(p) = Nbar(p);
//                for (int p = m; p < m_npr-1; ++p)
//                    temp(p) = Nbar(p+1);
//                Nbar = temp;
//                
//                for (int p = 0; p < m; ++p)
//                    temp(p) = m_Nbar(p);
//                for (int p = m; p < m_npr-1; ++p)
//                    temp(p) = m_Nbar(p+1);
//                m_Nbar = temp;
//                
//                // Update g
//                int size = mp_sizes[m+1]-mp_sizes[m];
//                temp = RealVector(m_nsr-size);
//                for (int j = 0; j < mp_sizes[m]; ++j)
//                    temp(j) = g(j);
//                for (int j = mp_sizes[m]; j < m_nsr-size; ++j)
//                    temp(j) = g(j+size);
//                g = temp;
//                
//                // Update dg
//                for (int j = 0; j < mp_sizes[m]; ++j)
//                    temp(j) = m_dg(j);
//                for (int j = mp_sizes[m]; j < m_nsr-size; ++j)
//                    temp(j) = m_dg(j+size);
//                m_dg = temp;
//                
//                // Update N
//                for (int j = 0; j < mp_sizes[m]; ++j)
//                    temp(j) = m_N(j);
//                for (int j = mp_sizes[m]; j < m_nsr-size; ++j)
//                    temp(j) = m_N(j+size);
//                m_N = temp;
//                
//                // Update the species order
//                removePhase(m);
//                
//                n = m_ncr+m_npr;
//                A = RealSymMat(n);
//                xt = RealVector(n);
//                dx = RealVector(n);
//                x  = RealVector(n);
//                r = RealVector(n);
//                
//                computeResidual(Nbar, r);
//                
//                x(0,m_ncr) = lambda;
//                x(m_ncr,n) = log(Nbar);
//            } else
//                m++;
//        }
//    
//        // Compute the jacobian matrix
//        formSystemMatrix(Nbar, A);
//#ifdef VERBOSE
//        cout << "Newton iter = " << iters << endl;
//        cout << "Jacobian = " << endl << A << endl;
//        cout << "r = " << endl << r << endl;
//#endif
////        // Solve the system to obtain dx
////        QRP<double> qrp(A);
////        r = -r;
////        qrp.solve(dx, r);
//        ldlt.setMatrix(A);
//        ldlt.solve(dx, -r);
////        SVD<double> svd(A);
////        r = -r;
////        svd.solve(dx, r);
//#ifdef VERBOSE
//        cout << "dx = " << endl << dx << endl;
//#endif
//        lam = 1.0;
//        xt = x + dx;
//        lambda = xt(0,m_ncr);
//        Nbar = exp(xt(m_ncr,n));
//        updateSpeciesMoles(lambda, Nbar, g);
//        computeResidual(Nbar, r);
//        res = r.norm2();
//#ifdef VERBOSE
//        cout << "x = " << endl << xt << endl;
//        cout << "res = " << res << endl;
//#endif
//        if (res > res0)
//            return std::make_pair(iters, res);
//        
//        int nsearch = 0;
//        while (res > (1.0-alpha*lam)*res0 && nsearch < 10) {
//            nsearch++;
//            lam *= 0.5;
//            xt = x + lam*dx;
//            lambda = xt(0,m_ncr);
//            Nbar = exp(xt(m_ncr,n));
//            if (!updateSpeciesMoles(lambda, Nbar, g))
//                continue;
//            computeResidual(Nbar, r);
//            res = r.norm2();
//#ifdef VERBOSE
//            cout << "line search, lamda = " << lam << ", res = " << res << endl;
//#endif
//        }
//        
//        if (nsearch == 10)
//            return std::make_pair(iters, res);
//
//        
//        x = xt;
//        res0 = res;
//    }
//    
//    return std::make_pair(iters, res0);
//}

//==============================================================================

void MultiPhaseEquilSolver::checkForDeterminedSpecies()
{
    // First determine which species must be zero (either due to constraints or
    // to temperature limits on condensed species)
    int species_group [m_ns];
    bool zero_constraint [m_nc];
    
    // Group each species
    for (int i = 0; i < m_ns; ++i) {
        species_group[i] = mp_phase[i];
        if (species_group[i] > 0 && !m_thermo.speciesThermoValidAtT(i, m_T))
            species_group[i] = m_np;
    }
    
    // Check for constraints which force species to be zero
    for (int i = 0; i < m_nc; ++i) {
        zero_constraint[i] = false;
        if (mp_c[i] == 0.0) {
            // Determine first non-zero sign
            bool pos;
            int j;
            for (j = 0; j < m_ns; ++j) {
                if (m_B(j,i) != 0.0) {
                    pos = m_B(j,i) > 0.0;
                    break;
                }
            }
            
            // Compare all the other signs to the first non-zero sign
            for ( ; j < m_ns; ++j)
                if ((m_B(j,i) != 0.0) && (pos != (m_B(j,i) > 0.0)))
                    break;
            
            if (j == m_ns) {
                // Keep track of the species which will be zeroed
                for (j = 0; j < m_ns; ++j) {
                    if (m_B(j,i) != 0.0)
                        species_group[j] = m_np;
                }
                
                // Keep track of the constraint that is zero
                zero_constraint[i] = true;
            }
        }
    }
    
    m_solution.setupOrdering(species_group, zero_constraint, m_pure_condensed);
}

//==============================================================================

void MultiPhaseEquilSolver::initialConditions()
{
    const double alpha = 1.0e-3;
    const int npr = m_solution.npr();
    const int ncr = m_solution.ncr();
    const int nsr = m_solution.nsr();
    const int* const p_sjr = m_solution.sjr();
    const int* const p_cir = m_solution.cir();
    const int* const p_sizes = m_solution.sizes();
    
    // Compute the min-g and max-min solutions
    updateMinGSolution(mp_g);
    updateMaxMinSolution();
    
    // Use the tableau as temporary storage
    double* p_N = mp_tableau;
    double* p_Nbar = p_N + nsr;
    double* p_lambda = p_Nbar + npr;
    
    // Form estimate of species moles from min-g and max-min vectors
    // (temporarily store in g0 vector)
    int j, jk;
    for (j = 0; j < nsr; ++j)
        p_N[j] = mp_ming[j]*(1.0-alpha) + mp_maxmin[j]*alpha;
    
    // Compute initial phase moles
    for (int m = 0, j = 0; m < npr; ++m) {
        p_Nbar[m] = 0.0;
        for ( ; j < p_sizes[m+1]; ++j)
            p_Nbar[m] += p_N[j];
    }
    
    // Compute RHS of lambda least-squares problem (temporarily store in g0)
    for (int m = 0, j = 0; m < npr; ++m) {
        const double nbar = p_Nbar[m];
        for ( ; j < p_sizes[m+1]; ++j) {
            p_N[j] = std::log(p_N[j] / nbar);
            mp_g0[j] = p_N[j] + mp_g[p_sjr[j]];
        }
    }
    
    // Compute SVD of Br
    SVD<double> svd(m_solution.reducedMatrix(m_B));
    svd.solve(p_lambda, mp_g0);
    
    // Now compute the g0 which satisfies the constraints
    for (int j = 0; j < nsr; ++j) {
        jk = p_sjr[j];
        mp_g0[jk] = -p_N[j];
        for (int i = 0; i < ncr; ++i)
            mp_g0[jk] += m_B(jk,p_cir[i])*p_lambda[i];
    }
    
    m_solution.setG(mp_g0, mp_g, 0.0);
    m_solution.setSolution(p_lambda, p_Nbar, m_B);
    
    // Make sure that we initially satisfy the phase rule
    while (m_solution.npr() > m_solution.ncr()) {
        double min_value = 0.0;
        int min_index = 1;
        for (int j = p_sizes[1]; j < p_sizes[2]; ++j)
            min_value += m_solution.y()[j];
        for (int m = 2; m < m_solution.npr(); ++m) {
            double value = 0.0;
            for (int j = p_sizes[m]; j < p_sizes[m+1]; ++j)
                value += m_solution.y()[j];
            if (value < min_value) {
                min_value = value;
                min_index = m;
            }
        }
        m_solution.removePhase(min_index);
    }
}

//void MultiPhaseEquilSolver::initialConditions(
//    RealVector& gtp, RealVector& lambda, RealVector& Nbar, RealVector& g)
//{    
//    //RealVector Nmm(m_nsr);
//    //RealVector Nmg(m_nsr);
//    
//    // Initial conditions on N (N = a*Nmg + (1-a)*Nmm)
//    //ming(m_cr, gtp, Nmg);
//    
//    //cout << "ming:" << endl;
//    //cout << Nmg << endl;
//    
//    //perturb(Nmm);
//    
//    //cout << "max-min:" << endl;
//    //cout << Nmm << endl;
//    
//    // Compute the min-g solution
//    //updateMinGSolution();
//    
//    // Check if the solid phase is actually present in the min-g solution
////    int m = 1;
////    while (m < m_npr) {
////        bool empty = true;
////        for (int j = mp_sizes[m]; j < mp_sizes[m+1]; ++j) {
////            if (mp_ming[j] > 0.0) {
////                empty = false;
////                break;
////            }
////        }
////        
////        if (!empty) {
////            m++;
////            continue;
////        }
////        
////        int size = mp_sizes[m+1]-mp_sizes[m];
////        for (int j = mp_sizes[m]; j < m_nsr - size; ++j)
////            mp_ming[j] = mp_ming[j+size];
////
////        RealVector temp(m_nsr-size);
////        for (int j = 0; j < mp_sizes[m]; ++j)
////            temp(j) = gtp(j);
////        for (int j = mp_sizes[m]; j < m_nsr - size; ++j)
////            temp(j) = gtp(j+size);
////        gtp = temp;
////        
////        temp = RealVector(m_npr-1);
////        for (int p = 0; p < m; ++p)
////            temp(p) = Nbar(p);
////        for (int p = m; p < m_npr-1; ++p)
////            temp(p) = Nbar(p+1);
////        Nbar = temp;
////        
////        removePhase(m);
////    }
//    
//    updateMaxMinSolution();
//    
//    const double alpha = 1.0E-3;
//    m_N = ones<double>(m_nsr);
//    for (int i = 0; i < m_nsr; ++i)
//        m_N(i) = (1.0-alpha)*mp_ming[i] + alpha*mp_maxmin[i];
//
//    
//    // Initial phase moles
//    for (int m = 0; m < m_npr; ++m) {
//        Nbar(m) = 0.0;
//        for (int j = mp_sizes[m]; j < mp_sizes[m+1]; ++j)
//            Nbar(m) += m_N(j);
//    }
//    
//    #ifdef VERBOSE
//        cout << "Initial phase moles:" << endl;
//        cout << Nbar << endl;
//        cout << "Initial species moles:" << endl;
//        cout << m_N << endl;
//    #endif
//    
//    // Mole fractions
//    for (int m = 0; m < m_npr; ++m) {
//        for (int j = mp_sizes[m]; j < mp_sizes[m+1]; ++j)
//            m_N(j) /= Nbar(m);
//    }
//    
//    #ifdef VERBOSE
//        cout << "Initial species mole fractions:" << endl;
//        cout << m_N << endl;
//    #endif
//    
//    RealMatrix Br(m_nsr,m_ncr);
//    for (int j = 0; j < m_nsr; ++j)
//        for (int i = 0; i < m_ncr; ++i)
//            Br(j,i) = m_B(mp_sjr[j], mp_cir[i]);
//    
//    #ifdef VERBOSE
//        cout << "Reduce B matrix in initialConditions()" << endl;
//        cout << Br << endl;
//    #endif
//    
//    // Initial lambda
//    SVD<double> svd(Br);
//    g = log(m_N) + gtp;
//    
//    #ifdef VERBOSE
//        cout << "RHS of least-squares equation:" << endl;
//        cout << g << endl;
//    #endif
//    lambda = ones<double>(m_ncr);
//    svd.solve(lambda, g);
//    
//    // Initial g
//    g = Br*lambda - log(m_N);
//    
//#ifdef VERBOSE
//    cout << "intial lambda: " << endl;
//    cout << lambda << endl;
//    cout << "initial Nbar: " << endl;
//    cout << Nbar << endl;
//    cout << "initial g:" << endl;
//    cout << g << endl;
//#endif
//}

//void MultiPhaseEquilSolver::perturb(RealVector& nmm)
//{
//    // Possibly perturb the constraints
//    double ne_max = m_cr.maxValue();
//    double ne_low = 1.0e-99 * ne_max;
//    m_cr = max(m_cr, ne_low);
//    
//    // Determine upper bound on moles of undetermined species (ie a species
//    // contains all the moles of a given element and thus cannot be any bigger)
//    RealVector nk_max(m_nsr, 0.0);
//    for (int i = 0; i < m_nsr; ++i) {
//        for (int j = 0; j < m_ncr; ++j)
//            nk_max(i) = std::max(nk_max(i), m_Br(i,j) / m_cr(j));
//        nk_max(i) = 1.0 / nk_max(i);
//    }
//
//    // Compute the max-min composition of the undetermined species using the
//    // perturbed constraints
//    maxmin(m_cr, nmm);
//    
//    // Compute a (possibly) perturbed vector of undetermined species moles
//    // (just storing in nk_max)
//    nk_max = max(nmm, nk_max * 1.0e-99);
//    
//    // Now back out the (possibly) perturbed constraint vector
//    m_cr = nk_max * m_Br;
//}

//==============================================================================

//void MultiPhaseEquilSolver::maxmin(const RealVector& c, RealVector& Nmm)
//{
//    // Setup the tableau for input to lp()
//    RealMatrix tableau(m_ncr, m_nsr+1);
//    
//    // B'
//    tableau(0,m_ncr,0,m_nsr) = m_Br.transpose();    
//    
//    // sum(B)'
//    for (size_t i = 0; i < m_ncr; ++i)
//        for (size_t j = 0; j < m_nsr; ++j)
//            tableau(i,m_nsr) += tableau(i,j);
//    
//    // Objective function
//    RealVector f(m_nsr+1, 0.0);
//    f(m_nsr) = 1.0;
//    
//    // Temporary solution vector (to make room for zmin)
//    RealVector x(m_nsr+1);
//    
//    // Now solve the LP problem
//    double zmin;
//    int* iposv;
//    LpResult lp_result = lp(f, MAXIMIZE, tableau, c, 0, 0, x, zmin, &iposv);
//    
//    if (lp_result != SOLUTION_FOUND && lp_result != LINEARLY_DEPENDENT) {
//        /*if (lp_result == LINEARLY_DEPENDENT) {
//            cout << "The following constraints are linearly dependent on the others:" << endl;
//            for (int i = 0; i < m_nc; ++i) {
//                if (iposv[i] > m_ns) {
//                    int index = iposv[i]-m_ns-1;
//                    if (index < m_ne)
//                        cout << m_thermo.elementName(index) << endl;
//                    else
//                        cout << "Constraint #" << index-m_ne+1 << endl;
//                }
//            }
//        } else {*/
//            cout << "Error finding max-min solution!  Either the solution is "
//                 << "not bounded or there is no feasible solution with the "
//                 << "given constraints. Exiting..." << endl;
//            exit(1);
//        //}
//    }
//    
//    delete [] iposv;
//    
//    // Transfer solution
//    Nmm = x(0,m_nsr) + zmin;
//}

//==============================================================================

void MultiPhaseEquilSolver::updateMinGSolution(const double* const p_g)
{
    // Use the current solution's ordering
    const int nsr = m_solution.nsr();
    const int ncr = m_solution.ncr();
    const int* const p_sjr = m_solution.sjr();
    const int* const p_cir = m_solution.cir();
    
    // Build the tableau
    // 0 -f'
    double* p = mp_tableau;
    *p++ = 0;
    for (int i = 0; i < nsr; ++i)
        *p++ = -p_g[p_sjr[i]];
    
    // c -B'
    int ik;
    for (int i = 0; i < ncr; ++i) {
        ik = p_cir[i];
        *p++ = mp_c[ik];
        for (int j = 0; j < nsr; ++j)
            *p++ = -m_B(p_sjr[j], ik);
    }
    
    // 0 0
    for (int i = 0; i < nsr+1; ++i)
        *p++ = 0.0;
    
    // Use the simplex algorithm to get the min-g solution
    int izrov [nsr];
    int iposv [ncr];
    int ret = simplex(mp_tableau, ncr, nsr, 0, 0, izrov, iposv, 1.0E-9);
    
    // Error check
    if (ret != 0) {
        cout << "Error in computing the min-g solution!" << endl;
        exit(1);
    }
    
    // Unravel the solution
    for (int i = 0; i < nsr; ++i)
        mp_ming[i] = 0.0;
    for (int i = 0; i < ncr; ++i) {
        if (iposv[i] < nsr)
            mp_ming[iposv[i]] = mp_tableau[(nsr+1)*(i+1)];
        else {
            cout << "Linearly dependent in min-g!" << endl;
            exit(1);
        }
    }
    
//    cout << "min-g solution:" << endl;
//    for (int i = 0; i < nsr; ++i)
//        cout << m_thermo.speciesName(p_sjr[i]) << " " << mp_ming[i] << endl;
}

void MultiPhaseEquilSolver::updateMaxMinSolution()
{
    // Use the current solution's ordering
    const int nsr = m_solution.nsr();
    const int ncr = m_solution.ncr();
    const int* const p_sjr = m_solution.sjr();
    const int* const p_cir = m_solution.cir();
    
    // Build the tableau
    // 0 0' 1
    double* p = mp_tableau;
    for (int i = 0; i < nsr+1; ++i)
        *p++ = 0.0;
    *p++ = 1.0;
    
    // c -B' -sum(B')
    double sum;
    int ik;
    for (int i = 0; i < ncr; ++i) {
        ik = p_cir[i];
        *p++ = mp_c[ik];
        sum = 0.0;
        for (int j = 0; j < nsr; ++j) {
            *p = -m_B(p_sjr[j], ik);
            sum += *p++;
        }
        *p++ = sum;
    }
    
    // 0 0 0
    for (int i = 0; i < nsr+2; ++i)
        *p++ = 0.0;
    
    // Use the simplex algorithm to get the max-min solution
    int izrov [nsr];
    int iposv [ncr];
    int ret = simplex(mp_tableau, ncr, nsr+1, 0, 0, izrov, iposv, 1.0E-9);
    
    // Error check
    if (ret != 0) {
        cout << "Error in computing the max-min solution!" << endl;
        exit(1);
    }
    
    // Unravel the solution
    for (int i = 0; i < nsr; ++i)
        mp_maxmin[i] = *mp_tableau;
    for (int i = 0; i < ncr; ++i) {
        if (iposv[i] < nsr)
            mp_maxmin[iposv[i]] += mp_tableau[(nsr+2)*(i+1)];
    }
    
//    cout << "max-min solution:" << endl;
//    for (int i = 0; i < nsr; ++i)
//        cout << m_thermo.speciesName(p_sjr[i]) << " " << mp_maxmin[i] << endl;
}

////void MultiPhaseEquilSolver::ming(
////    const RealVector& c, const RealVector& g, RealVector& Nmg)
////{
////    // Solve the linear programming problem
////    double temp;
////    int* iposv;
////    LpResult lp_result =
////        lp(g, MINIMIZE, m_Br.transpose(), c, 0, 0, Nmg, temp, &iposv);
////    
////    int nld = 0;
////    
////    if (lp_result != SOLUTION_FOUND) {
////        if (lp_result == LINEARLY_DEPENDENT) {
////            
////            bool lindep [m_ncr];
////            std::fill(lindep, lindep+m_ncr, false);
////            
////            for (int i = 0; i < m_ncr; ++i) {
////                if (iposv[i] >= m_nsr) {
////                    lindep[iposv[i]-m_nsr] = true;
////                    nld++;
////                }
////            }
////            
////            RealVector li(m_ncr-nld);
////            for (int i = 0, index = 0; i < m_ncr; ++i)
////                if (!lindep[i]) li(index++) = mp_cir[i];
////            
////            m_ncr -= nld;
////            m_cir = li;
////            
////            m_Br = zeros<double>(m_nsr, m_ncr);
////            for (int j = 0; j < m_nsr; ++j)
////                for (int i = 0; i < m_ncr; ++i)
////                    m_Br(j,i) = m_B(mp_sjr[j], mp_cir[i]);
////            
////            m_cr = ones<double>(m_ncr);
////            for (int i = 0; i < m_ncr; ++i)
////                m_cr(i) = m_c(mp_cir[i]);
////            
////        } else {
////            cout << "Error finding min-g solution!  Either the solution is "
////                 << "not bounded or there is no feasible solution with the "
////                 << "given constraints. " << endl;
////            cout << "Constraints:" << endl;
////            for (int i = 0; i < m_ne; ++i)
////                cout << setw(3) << m_thermo.elementName(i) << ": "
////                     << m_c(i) << endl;
////            for (int i = m_ne; i < m_nc; ++i)
////                cout << " C" << i-m_ne+1 << ": " << m_c(i) << endl;
////            exit(1);
////        }
////    }
////    
////    // Compute the base species
////    m_base = ones<int>(m_ncr);
////    for (int i = 0, index = 0; i < m_ncr+nld; ++i)
////        if (iposv[i] < m_nsr) m_base(index++) = iposv[i];
////    
////    delete [] iposv;
////}
//
////==============================================================================
//
void MultiPhaseEquilSolver::formSystemMatrix(RealSymMat& A) const
{
    const int npr = m_solution.npr();
    const int ncr = m_solution.ncr();
    const int nsr = m_solution.nsr();
    const int* const p_sjr = m_solution.sjr();
    const int* const p_cir = m_solution.cir();
    const int* const p_sizes = m_solution.sizes();
    const double* const p_y = m_solution.y();
    const double* const p_lnNbar = m_solution.lnNbar();
    
    A = 0.0;
    
    double temp, Nj;
    int mk, jk;
    
    for (int m = 0; m < npr; ++m) {
        mk = ncr+m;
        
        for (int j = p_sizes[m]; j < p_sizes[m+1]; ++j) {
            jk = p_sjr[j];
            Nj = p_y[j]*p_y[j];
            for (int i = 0; i < ncr; ++i) {
                temp = Nj * m_B(jk, p_cir[i]);
                for (int k = i; k < ncr; ++k)
                    A(i,k) += temp * m_B(jk, p_cir[k]);
                A(i,mk) += temp;
            }
            A(mk,mk) += Nj;
        }
        
        A(mk,mk) -= std::exp(p_lnNbar[m]);
    }
}
//
//
//
//
////==============================================================================
//
//bool MultiPhaseEquilSolver::updateSpeciesMoles(
//    const RealVector& lambda, const RealVector& Nbar, const RealVector& g)
//{
//    #ifdef VERBOSE
//    cout << "updateSpeciesMoles: " << endl;
//    #endif
//    
//    if (Nbar.sum() > m_c.sum()*10.0) {
//        #ifdef VERBOSE
//        cout << "Nbar to big: " << endl;
//        cout << Nbar << endl;
//        #endif
//        return false;
//    }
//    
//    int jk;
//    for (int j = 0; j < m_nsr; j++) {
//        m_N(j) = -g(j);
//        jk = mp_sjr[j];
//        for (int i = 0; i < m_ncr; ++i)
//            m_N(j) += m_B(jk,mp_cir[i])*lambda(i);
//    }
//    
//    for (int m = 0; m < m_npr; ++m) {
//        for (int j = mp_sizes[m]; j < mp_sizes[m+1]; ++j) {
//            if (m_N(j) > 10) {
//                #ifdef VERBOSE
//                cout << "N to big before exp():" << endl;
//                cout << m_N << endl;
//                #endif
//                return false;
//            }
//            else 
//                m_N(j) = std::exp(m_N(j)) * Nbar(m);
//        }
//    }
//    
//#ifdef VERBOSE
//    cout << "N = " << endl;
//    cout << m_N << endl;
//#endif
//    return true;
//}
//
////==============================================================================
//
void MultiPhaseEquilSolver::computeResidual(RealVector& r) const
{
    const int npr = m_solution.npr();
    const int ncr = m_solution.ncr();
    const int nsr = m_solution.nsr();
    const int* const p_sjr = m_solution.sjr();
    const int* const p_cir = m_solution.cir();
    const int* const p_sizes = m_solution.sizes();
    const double* const p_y = m_solution.y();
    const double* const p_lnNbar = m_solution.lnNbar();
    
    
    for (int i = 0; i < ncr; ++i)
        r(i) = -mp_c[p_cir[i]];
    
    int jk;
    for (int j = 0; j < nsr; ++j) {
        jk = p_sjr[j];
        for (int i = 0; i < ncr; ++i)
            r(i) += m_B(jk, p_cir[i]) * p_y[j] * p_y[j];
    }
    
    for (int m = 0; m < npr; ++m) {
        jk = ncr + m;
        r(jk) = -std::exp(p_lnNbar[m]);
        for (int j = p_sizes[m]; j < p_sizes[m+1]; j++)
            r(jk) += p_y[j] * p_y[j];
    }
}

//==============================================================================

//void MultiPhaseEquilSolver::updatePreConditioner()
//{    
//    const int nb = m_base.size();
//    
//    m_M = Identity<double>(nb+m_np);
//    
//    if (nb == 1) return;
//    
//    RealMatrix A(nb, nb);
//    RealVector b(nb);
//    RealVector x(nb);
//
//    // Loop over each row of the pre-conditioner and solve linear system
//    for (int n = 0; n < nb; ++n) {
//        for (int k = 0; k < nb; ++k) {
//            A(n,k) = m_Br(m_base(n),k);
//        }
//    }
//    
//    QRP<double> qrp(A);
//    
//    for (int n = 0; n < nb; ++n) {
//        b = 0.0; b(n) = 1.0;
//        qrp.solve(x, b);
//        
//        for (int k = 0; k < nb; ++k)
//            m_M(n,k) =
//                (std::abs(x(k)) > NumConst<double>::eps*10.0 ? x(k) : 0.0);
//    }
//    
//    cout << "Base species" << endl;
//    for (int n = 0; n < nb; ++n)
//        cout << m_thermo.speciesName(mp_sjr[m_base(n)]) << endl;
//    
//    cout << "Preconditioner" << endl;
//    cout << m_M;
// 
//}

//==============================================================================

    } // namespace Thermodynamics
} // namespace Mutation
