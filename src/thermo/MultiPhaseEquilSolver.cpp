
#include "MultiPhaseEquilSolver.h"
#include "minres.h"
#include "gmres.h"

#include <set>
#include <utility>

//#define VERBOSE

using std::cout;
using std::endl;
using std::setw;

using namespace Mutation::Numerics;

namespace Mutation {
    namespace Thermodynamics {

//==============================================================================

MultiPhaseEquilSolver::MultiPhaseEquilSolver(const Thermodynamics& thermo)
    : m_thermo(thermo), m_phase(m_thermo.nSpecies()), m_c(m_thermo.nElements())
{
    // Sizing information
    m_ns  = m_thermo.nSpecies();
    m_nsr = m_ns;
    m_ne  = m_thermo.nElements();
    m_nc  = m_ne;
    m_ncr = m_nc;
    
    // System element matrix
    m_B  = m_thermo.elementMatrix();
    m_Br = m_B;
    
    // Compute phase information
    initPhases();
}

//==============================================================================

void MultiPhaseEquilSolver::addConstraint(const double *const p_A)
{
    // Save constraints
    m_constraints.push_back(asVector(p_A, m_ns));
    
    // Update the B matrix
    m_B = zeros<double>(m_ns, ++m_nc);
    m_B(0,m_ns,0,m_ne) = m_thermo.elementMatrix();
    
    for (int i = 0; i < m_nc-m_ne; ++i)
        m_B.col(m_ne+i) = m_constraints[i];
}

//==============================================================================
        
void MultiPhaseEquilSolver::clearConstraints()
{
    m_constraints.clear();
    m_B = m_thermo.elementMatrix();
    m_nc = m_ne;
}

//==============================================================================
        
void MultiPhaseEquilSolver::dXdT(double *const p_dxdt) const
{
    // Compute the dg/dT = [-H/RT]/T values
    RealVector dgdt(m_nsr);
    m_thermo.speciesHOverRT(m_T, p_dxdt);
    for (int j = 0; j < m_nsr; ++j)
        dgdt(j) = -p_dxdt[m_sjr(j)] / m_T;
    
    // RHS
    double temp;
    RealVector rhs(m_ncr+m_np, 0.0);
    RealVector sol = rhs;
    
    for (int k = 0; k < m_nsr; k++) {
        temp = m_N(k)*dgdt(k);
        for (int i = 0; i < m_ncr; ++i)
            rhs(i) += temp * m_Br(k,i);
        rhs(m_ncr+m_phase(m_sjr(k))) += temp;
    }
    
    // Now compute the system solution
    RealSymMat A(m_ncr+m_np);
    formSystemMatrix(A);
    
    // Solve the system for dlambda/dT
    LeastSquares<double> ls(A);
    ls.solve(sol, rhs);
    
    // Compute dXj/dT = Xj*[-dgj/dT + sum_i Bji dlambdai/dT]
    RealVector Nbar(m_np, 0.0);
    for (int j = 0; j < m_nsr; ++j)
        Nbar(m_phase(m_sjr(j))) += m_N(j);
    
    std::fill(p_dxdt, p_dxdt+m_ns, 0.0);
    
    for (int j = 0; j < m_nsr; ++j) {
        p_dxdt[m_sjr(j)] = -dgdt(j);
        for (int i = 0; i < m_ncr; ++i)
            p_dxdt[m_sjr(j)] += m_Br(j,i) * sol(i);
        p_dxdt[m_sjr(j)] *= m_N(j) / Nbar(m_phase(m_sjr(j)));
    }
}

void MultiPhaseEquilSolver::dXdP(double *const p_dxdp) const
{
    // Compute the dg/dP = 1/P for gas species
    //RealVector dgdt(m_nsr, 1.0 / m_P);
    const double pinv = 1.0 / m_P;
    //for (int j = 0; j < m_nsr; ++j)
    //    dgdt = pinv;
    
    // RHS
    double temp;
    RealVector rhs(m_ncr+m_np, 0.0);
    RealVector sol = rhs;
    
    for (int k = 0; k < m_nsr; k++) {
        temp = m_N(k)*pinv;//dgdt(k);
        for (int i = 0; i < m_ncr; ++i)
            rhs(i) += temp * m_Br(k,i);
        rhs(m_ncr+m_phase(m_sjr(k))) += temp;
    }
    
    // Now compute the system solution
    RealSymMat A(m_ncr+m_np);
    formSystemMatrix(A);
    
    // Solve the system for dlambda/dT
    LeastSquares<double> ls(A);
    ls.solve(sol, rhs);
    
    // Compute dXj/dT = Xj*[-dgj/dT + sum_i Bji dlambdai/dT]
    RealVector Nbar(m_np, 0.0);
    for (int j = 0; j < m_nsr; ++j)
        Nbar(m_phase(m_sjr(j))) += m_N(j);
    
    std::fill(p_dxdp, p_dxdp+m_ns, 0.0);
    
    for (int j = 0; j < m_nsr; ++j) {
        p_dxdp[m_sjr(j)] = -pinv;//dgdt(j);
        for (int i = 0; i < m_ncr; ++i)
            p_dxdp[m_sjr(j)] += m_Br(j,i) * sol(i);
        p_dxdp[m_sjr(j)] *= m_N(j) / Nbar(m_phase(m_sjr(j)));
    }
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
    m_T = T;
    m_P = P;

#ifdef VERBOSE
    std::cout << "equilibrate(" << T << " K, " << P << " Pa,";
    for (int i = 0; i < m_ne; ++i)
        std::cout << m_thermo.elementName(i) << " " << p_cv[i] << (i == m_ne-1 ? ")" : ",");
    std::cout << std::endl;
#endif

    static RealVector dg;
    static RealVector g;
    static RealVector g0;
    
    RealVector Nbar(m_np);
    RealVector Nbar_trial;
    
    //RealVector lambda(m_nc);
    //RealVector lambda_trial(m_nc);
    
    //RealSymMat A(m_nc+m_np);
    //RealVector r(m_nc+m_np);
    //RealVector dx(m_nc+m_np);
    
    RealVector lambda(m_ncr);
    RealVector lambda_trial;
    
    // Begin by reducing zero species
    m_c = asVector(p_cv, m_nc);
    reduceZeroSpecies();
    
    // Resize vectors if need be
    dg = ones<double>(m_nsr);
    
    // Compute the unitless Gibbs function for each nonzero species
    // (temporarily store in dg)
    m_thermo.speciesGOverRT(T, P, p_sv);
    for (int i = 0; i < m_nsr; ++i)
        dg(i) = p_sv[m_sjr(i)];

    // Compute the initial conditions for the integration
    initialConditions(dg, lambda, Nbar, g0);
    dg -= g0;
    //updatePreConditioner();
    
    RealSymMat A(m_ncr+m_np);
    RealVector r(m_ncr+m_np);
    RealVector dx(m_ncr+m_np);
    
    double s  = 0.0;
    double ds = 1.0;
    
    int iter = 0;
    int newt = 0;
    
    bool update_dx = true;
    
    updateSpeciesMoles(lambda, Nbar, g0);
    
    // Integrate quantities from s = 0 to s = 1
    while (s < 1.0) {
        iter++;
        
#ifdef VERBOSE
        std::cout << "Step: iter = " << iter << ", s = " << s << ", ds = " << ds << ", newt = " << newt << endl;
#endif
        if (update_dx) {
            // Compute the system matrix
            formSystemMatrix(A);
            
            // Solve for dlambda/ds and dNbar/ds with given dg/ds
            double temp;
            r = 0.0;
            for (int k = 0; k < m_nsr; k++) {
                temp = m_N(k)*dg(k);
                for (int i = 0; i < m_ncr; ++i)
                    r(i) += temp * m_Br(k,i);
                r(m_ncr+m_phase(m_sjr(k))) += temp;
            }
  
            LeastSquares<double> ls(A);
            ls.solve(dx, r);
        }
        
        // Update trial values
        lambda_trial = lambda + dx(0,m_ncr)*ds;
        Nbar_trial   = Nbar * exp(dx(m_ncr,m_ncr+m_np)*ds);
        g            = g0 + dg*(s+ds);

        // Compute the 
        updateSpeciesMoles(lambda_trial, Nbar_trial, g);
        computeResidual(Nbar_trial, r);
        
        std::pair<int,double> res = 
            newton(lambda_trial, Nbar_trial, g, r, A, 1.0e-3, 10);
        newt += res.first;
        
        // If newton did not converge, reduce the step size and start over
        if (res.second > 1.0e-3) {
#ifdef VERBOSE
            std::cout << "Newton did not converge, reducing step size" << std::endl;
#endif
            ds *= 0.1;
            update_dx = false;
            continue;
        }
        
        // Update the system variables
        lambda = lambda_trial;
        Nbar   = Nbar_trial;
        s += ds;
        
        ds = std::min(ds*2.0, 1.0-s);
        if (ds == 0.0 && s < 1.0) {
            std::cerr << "Continuation step size dropped to zero in equil solver!" << endl;
            std::exit(1);
        }
        
        update_dx = true;
    }
    
    // Use Newton iterations to improve residual for final answer
    std::pair<int,double> res = newton(lambda, Nbar, g, r, A, 1.0e-16, 100);
    newt += res.first;
    
    // Compute the species mole fractions
    double moles = m_N.sum();
    std::fill(p_sv, p_sv+m_ns, 0.0);
    for (int i = 0; i < m_nsr; ++i)
        p_sv[m_sjr(i)] = m_N(i) / moles;
    
    /*N = max(N, 1.0e-99);
    Nbar = 0.0;
    for (int i = 0; i < m_ns; ++i)
        Nbar(m_phase(i)) += N(i);
    for (int i = 0; i < m_ns; ++i)
        p_sv[i] = N(i) / Nbar(m_phase(i));*/
    
    return std::make_pair(iter, newt);
}

//==============================================================================
std::pair<int,double> MultiPhaseEquilSolver::newton(
    RealVector& lambda, RealVector& Nbar, const RealVector& g, 
    RealVector& r, RealSymMat& A, const double tol, const int max_iters)
{
    
    int iterations = 0;
    double rnorm;
    RealVector dx(m_ncr+m_np);
    
    while ((rnorm = r.norm2()) > tol && iterations < max_iters) {
        iterations++;
#ifdef VERBOSE
        std::cout << "Newton: iter = " << iterations << ", r = " << rnorm << endl;
#endif
        formSystemMatrix(A);

        LeastSquares<double> ls(A);
        r = -r;
        ls.solve(dx, r);

        // Return if the MINRES failed to converge
        //if (hist.second > 1.0e-10)
        //    return std::make_pair(iterations,rnorm);
        
        lambda += dx(0,m_ncr);
        Nbar *= exp(dx(m_ncr,m_ncr+m_np));
        
#ifdef VERBOSE
        std::cout << "lambda:" << std::endl;
        std::cout << lambda << std::endl;
        std::cout << "Nbar:" << std::endl;
        std::cout << Nbar << std::endl;
#endif 
        
        //if (!
        updateSpeciesMoles(lambda, Nbar, g);//)
            //return std::make_pair(iterations, rnorm);
            
        computeResidual(Nbar, r);
    }
    
    return std::make_pair(iterations, rnorm);
}

//==============================================================================

void MultiPhaseEquilSolver::reduceZeroSpecies()
{
    bool zero_species    [m_ns];
    bool zero_constraint [m_nc];
    std::fill(zero_species, zero_species+m_ns, false);
    
    // Check for constraints which force species to be zero
    m_ncr = m_nc;
    for (int i = 0; i < m_nc; ++i) {
        zero_constraint[i] = false;
        if (m_c(i) == double(0)) {
            // Determine first non-zero sign
            bool pos;
            int j;
            for (j = 0; j < m_ns; ++j) {
                if (m_B(j,i) != double(0)) {
                    pos = m_B(j,i) > double(0);
                    break;
                }
            }
            
            // Compare all the other signs to the first non-zero sign
            for ( ; j < m_ns; ++j)
                if ((m_B(j,i) != double(0)) && (pos != (m_B(j,i) > double(0))))
                    break;
            
            if (j == m_ns) {
                // Keep track of the species which will be zeroed
                for (j = 0; j < m_ns; ++j)
                    zero_species[j] |= (m_B(j,i) != double(0));
                
                // Keep track of the constraint that is zero
                zero_constraint[i] = true;
                m_ncr--;
            }
        }
    }
    
    // Count the number of nonzero species
    m_nsr = m_ns;
    for (int i = 0; i < m_ns; ++i)
        m_nsr -= (zero_species[i] ? 1 : 0);
    
    // Update list of reduced species indices
    m_sjr = ones<int>(m_nsr);
    int index = 0;
    for (int i = 0; i < m_ns; ++i)
        if (!zero_species[i]) m_sjr(index++) = i;
    
    // Update list of reduced constraint indices
    m_cir = ones<int>(m_ncr);
    index = 0;
    for (int i = 0; i < m_nc; ++i)
        if (!zero_constraint[i]) m_cir(index++) = i;
    
    // Compute reduced B matrix and constraint vector
    if (m_ncr != m_nc) {
        m_Br = zeros<double>(m_nsr, m_ncr);
        for (int j = 0; j < m_nsr; ++j)
            for (int i = 0; i < m_ncr; ++i)
                m_Br(j,i) = m_B(m_sjr(j), m_cir(i));
        
        m_cr = ones<double>(m_ncr);
        for (int i = 0; i < m_ncr; ++i)
            m_cr(i) = m_c(m_cir(i));
    } else {
        m_Br = m_B;
        m_cr = m_c;
    }
}

//==============================================================================

void MultiPhaseEquilSolver::initialConditions(
    const RealVector& gtp, RealVector& lambda, RealVector& Nbar, RealVector& g)
{    
    RealVector Nmm(m_nsr);
    RealVector Nmg(m_nsr);
    m_N = ones<double>(m_nsr);
    
    // Initial conditions on N (N = a*Nmg + (1-a)*Nmm)
    ming(m_cr, gtp, Nmg);
    perturb(Nmm);
    
    double alpha = 0.999;
    m_N = max(alpha*Nmg + (1.0-alpha)*Nmm, 1.0e-200);

    
    // Initial phase moles
    Nbar = 0.0;
    for (int i = 0; i < m_nsr; ++i)
        Nbar(m_phase(m_sjr(i))) += m_N(i);
    
    // Mole fractions
    for (int i = 0; i < m_nsr; ++i)
        m_N(i) /= Nbar(m_phase(m_sjr(i)));
    
    // Initial lambda
    SVD<double> svd(m_Br);
    g = log(m_N) + gtp;
    lambda = ones<double>(m_ncr);
    svd.solve(lambda, g);
    
    // Initial g
    g = m_Br*lambda - log(m_N);
}

void MultiPhaseEquilSolver::perturb(RealVector& nmm)
{
    // Possibly perturb the constraints
    double ne_max = m_cr.maxValue();
    double ne_low = 1.0e-200 * ne_max;    
    m_cr = Numerics::max(m_cr, ne_low);
    
    // Determine upper bound on moles of undetermined species (ie a species
    // contains all the moles of a given element and thus cannot be any bigger)
    RealVector nk_max(m_nsr, 0.0);
    for (int i = 0; i < m_nsr; ++i) {
        for (int j = 0; j < m_ncr; ++j)
            nk_max(i) = std::max(nk_max(i), (m_cr(j) == 0.0 ? 0.0 : (m_Br(i,j) / m_cr(j))));
        nk_max(i) = 1.0 / nk_max(i);
    }

    // Compute the max-min composition of the undetermined species using the
    // perturbed constraints
    maxmin(m_cr, nmm);
    
    // Compute a (possibly) perturbed vector of undetermined species moles
    // (just storing in nk_max)
    nk_max = max(nmm, nk_max * 1.0e-200);
    
    // Now back out the (possibly) perturbed constraint vector
    m_cr = nk_max * m_Br;
}

//==============================================================================

void MultiPhaseEquilSolver::maxmin(const RealVector& c, RealVector& Nmm)
{
    // Setup the tableau for input to lp()
    RealMatrix tableau(m_ncr, m_nsr+1);
    
    // B'
    tableau(0,m_ncr,0,m_nsr) = m_Br.transpose();    
    
    // sum(B)'
    for (size_t i = 0; i < m_ncr; ++i)
        for (size_t j = 0; j < m_nsr; ++j)
            tableau(i,m_nsr) += tableau(i,j);
    
    // Objective function
    RealVector f(m_nsr+1, 0.0);
    f(m_nsr) = 1.0;
    
    // Temporary solution vector (to make room for zmin)
    RealVector x(m_nsr+1);
    
    // Now solve the LP problem
    double zmin;
    int* iposv;
    LpResult lp_result = lp(f, MAXIMIZE, tableau, c, 0, 0, x, zmin, &iposv);
    
    if (lp_result != SOLUTION_FOUND && lp_result != LINEARLY_DEPENDENT) {
        /*if (lp_result == LINEARLY_DEPENDENT) {
            cout << "The following constraints are linearly dependent on the others:" << endl;
            for (int i = 0; i < m_nc; ++i) {
                if (iposv[i] > m_ns) {
                    int index = iposv[i]-m_ns-1;
                    if (index < m_ne)
                        cout << m_thermo.elementName(index) << endl;
                    else
                        cout << "Constraint #" << index-m_ne+1 << endl;
                }
            }
        } else {*/
            cout << "Error finding max-min solution!  Either the solution is "
                 << "not bounded or there is no feasible solution with the "
                 << "given constraints. Exiting..." << endl;
            exit(1);
        //}
    }
    
    delete [] iposv;
    
    // Transfer solution
    Nmm = x(0,m_nsr) + zmin;
}

//==============================================================================

void MultiPhaseEquilSolver::ming(
    const RealVector& c, const RealVector& g, RealVector& Nmg)
{
    // Solve the linear programming problem
    double temp;
    int* iposv;
    LpResult lp_result =
        lp(g, MINIMIZE, m_Br.transpose(), c, 0, 0, Nmg, temp, &iposv);
    
    int nld = 0;
    
    if (lp_result != SOLUTION_FOUND) {
        if (lp_result == LINEARLY_DEPENDENT) {
            
            bool lindep [m_ncr];
            std::fill(lindep, lindep+m_ncr, false);
            
            for (int i = 0; i < m_ncr; ++i) {
                if (iposv[i] >= m_nsr) {
                    lindep[iposv[i]-m_nsr] = true;
                    nld++;
                }
            }
            
            RealVector li(m_ncr-nld);
            for (int i = 0, index = 0; i < m_ncr; ++i)
                if (!lindep[i]) li(index++) = m_cir(i);
            
            m_ncr -= nld;
            m_cir = li;
            
            m_Br = zeros<double>(m_nsr, m_ncr);
            for (int j = 0; j < m_nsr; ++j)
                for (int i = 0; i < m_ncr; ++i)
                    m_Br(j,i) = m_B(m_sjr(j), m_cir(i));
            
            m_cr = ones<double>(m_ncr);
            for (int i = 0; i < m_ncr; ++i)
                m_cr(i) = m_c(m_cir(i));
            
        } else {
            cout << "Error finding min-g solution!  Either the solution is "
                 << "not bounded or there is no feasible solution with the "
                 << "given constraints. Exiting..." << endl;
            exit(1);
        }
    }
    
    // Compute the base species
    m_base = ones<int>(m_ncr);
    for (int i = 0, index = 0; i < m_ncr+nld; ++i)
        if (iposv[i] < m_nsr) m_base(index++) = iposv[i];
    
    delete [] iposv;
}

//==============================================================================

void MultiPhaseEquilSolver::formSystemMatrix(RealSymMat& A) const
{
    A = 0.0;
    double temp;
    
    for (int i = 0; i < m_ncr; ++i) {
        for (int k = 0; k < m_nsr; ++k) {
            temp = m_Br(k,i)*m_N(k);
            
            for (int j = i; j < m_ncr; ++j)
                A(i,j) += m_Br(k,j)*temp;
            
            A(i,m_phase(m_sjr(k))+m_ncr) += temp;
        }
    }
}

//==============================================================================

bool MultiPhaseEquilSolver::updateSpeciesMoles(
    const RealVector& lambda, const RealVector& Nbar, const RealVector& g)
{
    m_N = m_Br*lambda - g;
    m_N.exp();
    
    for (int i = 0; i < m_nsr; ++i)
        m_N(i) *= Nbar(m_phase(m_sjr(i)));
    
    return true;
}

//==============================================================================

void MultiPhaseEquilSolver::computeResidual(
    const RealVector& Nbar, RealVector& r) const
{
    r(0,m_ncr) = -m_cr;
    for (int i = 0; i < m_nsr; ++i) {
        if (m_N(i) > 0.0) {
            for (int j = 0; j < m_ncr; ++j)
                r(j) += m_Br(i,j) * m_N(i);
        }
    }
    
    r(m_ncr,m_ncr+m_np) = -Nbar;
    for (int i = 0; i < m_nsr; ++i)
        r(m_ncr+m_phase(m_sjr(i))) += m_N(i);
}

//==============================================================================

/*void MultiPhaseEquilSolver::updatePreConditioner()
{    
    const int nb = m_base.size();
    
    m_P = Identity<double>(nb+m_np);
    
    if (nb == 1) return;
    
    RealMatrix A(nb, nb);
    RealVector b(nb);
    RealVector x(nb);

    // Loop over each row of the pre-conditioner and solve linear system
    for (int n = 0; n < nb; ++n) {
        for (int k = 0; k < nb; ++k) {
            A(n,k) = m_Br(m_base(n),k);
        }
    }
    
    QRP<double> qrp(A);
    
    for (int n = 0; n < nb; ++n) {
        b = 0.0; b(n) = 1.0;
        qrp.solve(x, b);
        
        for (int k = 0; k < nb; ++k)
            m_P(n,k) = 
                (std::abs(x(k)) > NumConst<double>::eps*10.0 ? x(k) : 0.0);
    }
    
    //cout << "Base species" << endl;
    //for (int n = 0; n < nb; ++n)
    //    cout << m_thermo.speciesName(m_sjr(m_base(n))) << endl;
    //
    //cout << "Preconditioner" << endl;
    //cout << m_P;
    
}*/

//==============================================================================

    } // namespace Thermodynamics
} // namespace Mutation
