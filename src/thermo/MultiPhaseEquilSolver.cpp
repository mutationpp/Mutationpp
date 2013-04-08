
#include "MultiPhaseEquilSolver.h"
#include "minres.h"

#include <set>
#include <utility>

using std::cout;
using std::endl;

using namespace Mutation::Numerics;

namespace Mutation {
    namespace Thermodynamics {

//==============================================================================

MultiPhaseEquilSolver::MultiPhaseEquilSolver(const Thermodynamics& thermo)
    : m_thermo(thermo), m_phase(m_thermo.nSpecies())
{
    // Sizing information
    m_ns = m_thermo.nSpecies();
    m_ne = m_thermo.nElements();
    m_nc = m_ne;
    
    // System element matrix
    m_B  = m_thermo.elementMatrix();
    
    // Compute phase information
    initPhases();
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
    static RealVector dg;
    static RealVector g;
    static RealVector c;
    static RealVector N;
    static RealVector Nbar(m_np, 0.0);
    static RealVector lambda(m_nc);
    static RealSymMat A(m_nc+m_np);
    static RealVector r(m_nc+m_np);
    
    // Compute the unitless Gibbs function for each species
    // (temporarily store in N because we only need it to compute g and dg)
    m_thermo.speciesGOverRT(T, P, p_sv);
    N = asVector(p_sv, m_ns);
    
    // Compute the maxmin composition
    c = asVector(p_ev, m_nc);

    // Compute the initial conditions for the integration
    initialConditions(c, N, lambda, Nbar, g);
    dg = N - g;
    
    double s = 0.0;
    double ds = 0.2;
    int iter = 0;
    int newt = 0;
    
    computeSpeciesMoles(lambda, Nbar, g, N);
    
    // Integrate quantities from s = 0 to s = 1
    while (s < 1.0) {
        iter++;
        
        // Compute the system matrix
        formSystemMatrix(N, A);
        
        // Solve for dlambda/ds and dNbar/ds with given dg/ds
        double temp;
        r = 0.0;
        for (int k = 0; k < m_ns; k++) {
            temp = N(k)*dg(k);
            
            for (int i = 0; i < m_nc; ++i)
                r(i) += temp * m_B(k,i);
            
            r(m_nc+m_phase(k)) += temp;
        }
        
        minres(A, r, r);
        
        // Update values
        lambda += r(0,m_nc)*ds;
        Nbar *= exp(r(m_nc,m_nc+m_np)*ds);
        g += dg*ds;
        s += ds; 
        
        computeSpeciesMoles(lambda, Nbar, g, N);
        computeResidual(c, Nbar, N, r);
        
        // Newton Iteration to reduce residual to acceptable tolerance
        while (r.norm2() > 1.0e-3) {
            newt++;
            formSystemMatrix(N, A);
            minres(A, r, -1.0*r);
            
            lambda += r(0,m_nc);
            Nbar *= exp(r(m_nc,m_nc+m_np));
            
            computeSpeciesMoles(lambda, Nbar, g, N);
            computeResidual(c, Nbar, N, r);
        }
    }
    
    // Use Newton iterations to improve residual for final answer
    while (r.norm2() > 1.0e-8) {
        newt++;
        formSystemMatrix(N, A);
        minres(A, r, -1.0*r);
        
        lambda += r(0,m_nc);
        Nbar *= exp(r(m_nc,m_nc+m_np));
        
        computeSpeciesMoles(lambda, Nbar, g, N);
        computeResidual(c, Nbar, N, r);
    }
    
    // Compute the species mole fractions
    double moles = Nbar.sum();
    for (int i = 0; i < m_ns; ++i)
        p_sv[i] = N(i) / moles;
    
    return std::make_pair(iter, newt);
}

//==============================================================================

void MultiPhaseEquilSolver::initialConditions(
    const RealVector& c, const RealVector& gtp, RealVector& lambda, 
    RealVector& Nbar, RealVector& g) const
{
    static RealVector Nmm(m_ns);
    static RealVector Nmg(m_ns);
    static RealVector N(m_ns);
    
    // Initial conditions on N (N = a*Nmg + (1-a)*Nmm)
    maxmin(c, Nmm);
    ming(c, gtp, Nmg);
    
    double alpha = 0.999;
    N = alpha*Nmg + (1.0-alpha)*Nmm;
    
    // Initial phase moles
    Nbar = 0.0;
    for (int i = 0; i < m_ns; ++i)
        Nbar(m_phase(i)) += N(i);
    
    // Mole fractions
    for (int i = 0; i < m_ns; ++i)
        N(i) /= Nbar(m_phase(i));
    
    // Initial lambda
    static QRP<double> qrB(m_B);
    g = log(N) + gtp;
    qrB.solve(lambda, g);
    
    // Initial g
    g = m_B*lambda - log(N);
}

//==============================================================================

void MultiPhaseEquilSolver::maxmin(const RealVector& c, RealVector& Nmm) const
{
    // Setup the tableau for input to lp()
    RealMatrix tableau(m_nc, m_ns+1);
    
    // B'
    tableau(0,m_nc,0,m_ns) = m_B.transpose();    
    
    // sum(B)'
    for (size_t i = 0; i < m_nc; ++i)
        for (size_t j = 0; j < m_ns; ++j)
            tableau(i,m_ns) += tableau(i,j);
    
    // Objective function
    RealVector f(m_ns+1, 0.0);
    f(m_ns) = 1.0;
    
    // Temporary solution vector (to make room for zmin)
    RealVector x(m_ns+1);
    
    // Now solve the LP problem
    double zmin;
    LpResult lp_result = lp(f, MAXIMIZE, tableau, c, 0, 0, x, zmin);
    
    if (lp_result != SOLUTION_FOUND) {
        cout << "Error finding max-min solution!" << endl;
        exit(1);
    }
    
    // Transfer solution
    Nmm = x(0,m_ns) + zmin;
}

//==============================================================================

void MultiPhaseEquilSolver::ming(
    const RealVector& c, const RealVector& g, RealVector& Nmg) const
{
    // Solve the linear programming problem
    double temp;
    LpResult lp_result =
        lp(g, MINIMIZE, m_B.transpose(), c, 0, 0, Nmg, temp);
    
    if (lp_result != SOLUTION_FOUND) {
        cout << "Error finding min-g solution!" << endl;
        exit(1);
    }
}

//==============================================================================

void MultiPhaseEquilSolver::formSystemMatrix(
    const RealVector& N, RealSymMat& A) const
{
    // nc*nc*ns*3 + nc*ns*2 = nc*ns(2+3*nc)
    // dW/dlambda
    A = 0.0;
    /*for (int i = 0; i < m_nc; ++i)
        for (int j = i; j < m_nc; ++j)
            for (int k = 0; k < m_ns; ++k)
                A(i,j) += m_B(k,i)*m_B(k,j)*N(k);
    
    // dW/dNbar
    for (int i = 0; i < m_nc; ++i)
        for (int k = 0; k < m_ns; ++k)
            A(i,m_phase(k)+m_nc) += m_B(k,i)*N(k);*/
            
    
    // nc*ns*(2+2*nc), saves nc*nc*ns operations
    double temp;
    for (int i = 0; i < m_nc; ++i) {
        for (int k = 0; k < m_ns; ++k) {
            temp = m_B(k,i)*N(k);
            
            for (int j = i; j < m_nc; ++j)
                A(i,j) += m_B(k,j)*temp;
            
            A(i,m_phase(k)+m_nc) += temp;
        }
    }
    
    // Matrix conditioning
    //for (int i = 0; i < m_np; ++i)
    //    if (Nbar(i) < 1.0e-6) A(m_nc+i,m_nc+i) = 1.0;
}

//==============================================================================

void MultiPhaseEquilSolver::computeSpeciesMoles(
    const RealVector& lambda, const RealVector& Nbar, const RealVector& g,
    RealVector& N) const
{
    N = exp(m_B*lambda - g);
    for (int i = 0; i < m_ns; ++i)
        N(i) *= Nbar(m_phase(i));
}

//==============================================================================

void MultiPhaseEquilSolver::computeResidual(
    const RealVector& c, const RealVector& Nbar, const RealVector& N, 
    RealVector& r) const
{
    r(0,m_nc) = N*m_B - c;
    r(m_nc,m_nc+m_np) = -1.0*Nbar;
    for (int i = 0; i < m_ns; ++i)
        r(m_nc+m_phase(i)) += N(i);
}

//==============================================================================

    } // namespace Thermodynamics
} // namespace Mutation
