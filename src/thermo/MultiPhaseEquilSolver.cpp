
#include "MultiPhaseEquilSolver.h"
#include "minres.h"
#include "gmres.h"

#include <set>
#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>

//#define VERBOSE
#define SAVE_PATH
#include "Utilities.h"

using std::cout;
using std::endl;
using std::setw;

using namespace Mutation::Numerics;

namespace Mutation {
    namespace Thermodynamics {
    

// Define static constants for class MultiPhaseEquilSolver
const double MultiPhaseEquilSolver::ms_eps_rel = 0.001;
const double MultiPhaseEquilSolver::ms_eps_abs = 1.0e-8;
const double MultiPhaseEquilSolver::ms_ds_inc  = 4.0;
const double MultiPhaseEquilSolver::ms_ds_dec  = 0.25;
const double MultiPhaseEquilSolver::ms_max_ds  = 0.01;


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
#ifdef SAVE_PATH
    std::ofstream of("path.dat");
#endif

    // Save T, P, c
    m_T  = T;
    m_P  = P;
    std::copy(p_cv, p_cv+m_nc, mp_c);
    /*#ifdef VERBOSE
    std::cout << "equilibrate(" << T << " K, " << P << " Pa,";
    for (int i = 0; i < m_ne; ++i)
        std::cout << m_thermo.elementName(i) << " " << p_cv[i] << (i == m_ne-1 ? ")" : ",");
    std::cout << std::endl;
    #endif*/
    DEBUG("equilibrate(" << T << " K, " << P << " Pa,")
    for (int i = 0; i < m_ne; ++i)
		DEBUG(m_thermo.elementName(i) << " " << p_cv[i] << (i == m_ne-1 ? ")" : ","))
	DEBUG(endl)
    
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
    double ds = ms_max_ds;
    double res, resk = ms_eps_abs;
    Solution last_solution(m_thermo);
    
    // Integrate the equilibrium solution from s = 0 to s = 1
    m_niters = 0;
    m_nnewts = 0;
    while (s < 1.0) {
        // Save current solution in case we need to take a smaller step
        last_solution = m_solution;
        
        // Compute the rates for lambda and Nbar
        RealVector dx(m_solution.ncr()+m_solution.npr());
        rates(dx);
        #ifdef VERBOSE
        cout << "iter = " << m_niters << ", s = " << s << ", ds = " << ds << endl;
        cout << "dx = " << endl;
        cout << dx << endl;
        #endif

#ifdef SAVE_PATH
        // Compute the Gibbs free energy of the system
        of << setw(12) << s;

        double G = 0.0;
        for (int i = 0; i < m_solution.ncr(); ++i)
            G += m_solution.lambda()[i] * mp_c[m_solution.cir()[i]];
        of << setw(12) << G;

        for (int i = 0; i < m_solution.ncr(); ++i)
            of << setw(12) << m_solution.lambda()[i];
        for (int m = 0; m < m_solution.npr(); ++m)
            of << setw(12) << std::exp(m_solution.lnNbar()[m]);

        of << endl;
#endif

        // Keep trying to obtain a solution at progressively smaller steps
        // until a solution is obtained
        bool solution_obtained = false;
        while (!solution_obtained) {
            // Get trial solution
            m_solution.setG(mp_g0, mp_g, s+ds);
            m_solution.update(dx*ds, m_B);

            // Use newton to reduce the residual
            res = newton();
            solution_obtained =
            		(res < std::max((1.0+ms_eps_rel)*resk, ms_eps_abs));

            // If newton fails to converge then reduce the step size and try
            // again, otherwise we can continue
            if (solution_obtained)
            	resk = res;
            else {
                #ifdef VERBOSE
                cout << "Could not converge to solution, reducing step size" << endl;
                #endif
                ds *= ms_ds_dec;
                m_solution = last_solution;
            }
            
            if (m_niters + m_nnewts > 1000) {
                cout << "equil iter: " << m_niters << ", newt: " << m_nnewts << endl;
                cout << "s = " << s << ", ds = " << log(ds) << ", res = " << res << ", resk = " << resk << endl;
                std::cout << "equilibrate(" << T << " K, " << P << " Pa,";
                    for (int i = 0; i < m_ne; ++i)
                        std::cout << m_thermo.elementName(i) << " " << p_cv[i] << (i == m_ne-1 ? ")" : ",");
                    std::cout << std::endl;
                exit(1);
            }
        }
        
        // We found the solution at s, increase step size and move forward
        #ifdef VERBOSE
        cout << "Step succeded, increasing step size" << endl;
        #endif

        s += ds;
        ds = std::min(ms_max_ds, std::min(ds*ms_ds_inc, 1.0-s));
        m_niters++;
    }
    
    // Check one last newton
    resk = newton();

#ifdef SAVE_PATH
    // Compute the Gibbs free energy of the system
    of << setw(12) << s;

    double G = 0.0;
    for (int i = 0; i < m_solution.ncr(); ++i)
        G += m_solution.lambda()[i] * mp_c[m_solution.cir()[i]];
    of << setw(12) << G;

    for (int i = 0; i < m_solution.ncr(); ++i)
        of << setw(12) << m_solution.lambda()[i];
    for (int m = 0; m < m_solution.npr(); ++m)
        of << setw(12) << std::exp(m_solution.lnNbar()[m]);

    of << endl;
#endif

    if (resk > ms_eps_abs) {
    	cout << "Warning: equilibrium solver finished with residual of "
    	     << resk << "!";
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
    
#ifdef SAVE_PATH
    of.close();
#endif

    return std::make_pair(m_niters, m_nnewts);
}

//==============================================================================

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

//==============================================================================

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

double MultiPhaseEquilSolver::newton()
{
    #ifdef VERBOSE
    cout << "newton:" << endl;
    #endif
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
        return res;
    
    int iter = 0;
    while (res > ms_eps_abs && iter < max_iters) {
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
        
        // Solve the linear system (if it is singular then don't bother)
        if (ldlt.setMatrix(A))
        	ldlt.solve(dx, -r);
        else
			return res;
        //QRP<double> qrp(A);
        //#ifdef VERBOSE
        //cout << "R matrix = " << endl;
        //cout << qrp.R() << endl;
        //#endif
        //r = -r;
        //qrp.solve(dx, r);
        
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
    return res;
}

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

//==============================================================================

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

//==============================================================================

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

    } // namespace Thermodynamics
} // namespace Mutation
