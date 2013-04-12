
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
    : m_thermo(thermo), m_phase(m_thermo.nSpecies()), m_c(m_thermo.nElements())
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

void MultiPhaseEquilSolver::addConstraint(
    const double *const p_A, const double c)
{
    // Save constraints
    m_constraints.push_back(asVector(p_A, m_ns));
    
    // Update the B matrix
    m_B = zeros<double>(m_ns, ++m_nc);
    m_B(0,m_ns,0,m_ne) = m_thermo.elementMatrix();
    
    for (int i = 0; i < m_nc-m_ne; ++i)
        m_B.col(m_ne+i) = m_constraints[i];
    
    // Update the constraint vector
    RealVector temp(m_nc);
    temp(m_ne,m_nc-1) = m_c(m_ne,m_nc-1);
    temp(m_nc-1) = c;
    m_c = temp;
}

//==============================================================================
        
void MultiPhaseEquilSolver::clearConstraints()
{
    m_constraints.clear();
    m_B = m_thermo.elementMatrix();
    RealVector temp(m_ne);
    m_c = temp;
    m_nc = m_ne;
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
    double T, double P, const double *const p_cv, double *const p_sv)
{
    static RealVector dg;
    static RealVector g;
    static RealVector N(m_ns);
    static RealVector Nbar(m_np, 0.0);
    RealVector lambda(m_nc);
    RealSymMat A(m_nc+m_np);
    RealVector r(m_nc+m_np);
    
    // Compute the unitless Gibbs function for each species
    // (temporarily store in dg)
    m_thermo.speciesGOverRT(T, P, p_sv);
    dg = asVector(p_sv, m_ns);
    
    // Compute the maxmin composition
    m_c(0,m_ne) = asVector(p_cv, m_ne);

    // Compute the initial conditions for the integration
    initialConditions(dg, lambda, Nbar, g);
    dg -= g;
    
    double s = 0.0;
    double ds = 0.01;
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
        ds = std::min(ds*1.1, 1.0-s);
        
        computeSpeciesMoles(lambda, Nbar, g, N);
        computeResidual(Nbar, N, r);
        
        newt += newton(lambda, Nbar, g, N, r, A, 1.0e-3);
    }
    
    // Use Newton iterations to improve residual for final answer
    newt += newton(lambda, Nbar, g, N, r, A, 1.0e-12);
    
    // Compute the species mole fractions
    double moles = N.sum();
    for (int i = 0; i < m_ns; ++i)
        p_sv[i] = N(i) / moles;
    
    return std::make_pair(iter, newt);
}

//==============================================================================

int MultiPhaseEquilSolver::newton(
    RealVector& lambda, RealVector& Nbar, const RealVector& g, RealVector& N, 
    RealVector& r, RealSymMat& A, const double tol)
{
    int iterations = 0;
    
    while (r.norm2() > tol) {
        iterations++;
        formSystemMatrix(N, A);
        minres(A, r, -1.0*r);
        
        lambda += r(0,m_nc);
        Nbar *= exp(r(m_nc,m_nc+m_np));
        
        computeSpeciesMoles(lambda, Nbar, g, N);
        computeResidual(Nbar, N, r);
    }
    
    return iterations;
}

//==============================================================================

void MultiPhaseEquilSolver::initialConditions(
    const RealVector& gtp, RealVector& lambda, RealVector& Nbar, RealVector& g) 
    const
{
    static RealVector Nmm(m_ns);
    static RealVector Nmg(m_ns);
    static RealVector N(m_ns);
    
    // Initial conditions on N (N = a*Nmg + (1-a)*Nmm)
    maxmin(m_c, Nmm);
    ming(m_c, gtp, Nmg);
    
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
    QRP<double> qrB(m_B);
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
    A = 0.0;
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
    const RealVector& Nbar, const RealVector& N, RealVector& r) const
{
    r(0,m_nc) = N*m_B - m_c;
    r(m_nc,m_nc+m_np) = -1.0*Nbar;
    for (int i = 0; i < m_ns; ++i)
        r(m_nc+m_phase(i)) += N(i);
}

//==============================================================================

    } // namespace Thermodynamics
} // namespace Mutation
