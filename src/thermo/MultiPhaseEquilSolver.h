/**
 * @file MultiPhaseEquilSolver.h
 *
 * Defines the MultiPhaseEquilSolver class which implements a Continuation-
 * Krylov based multiphase equilibrium solver.
 *
 * @author J.B. Scoggins (scoggins@vki.ac.be)
 */

#ifndef THERMO_MULTI_PHASE_EQUIL_SOLVER
#define THERMO_MULTI_PHASE_EQUIL_SOLVER

#include "Thermodynamics.h"
#include "Numerics.h"

namespace Mutation {
    namespace Thermodynamics {

class MultiPhaseEquilSolver
{
public:
    /**
     * Constructs the equilibrium solver.
     */
    MultiPhaseEquilSolver(const Thermodynamics& thermo);
    
    /**
     * Destructor.
     */
    ~MultiPhaseEquilSolver() {};
    
    /**
     * Returns number of phases considered by the equilibrium solver.
     */
    int nPhases() const {
        return m_np;
    }
    
    /**
     * Computes the equilibrium composition of the mixture.
     */
    std::pair<int,int> equilibrate(
        double T, double P, const double *const p_ev, double *const p_sv);
    
    /**
     * Adds an additional linear constraint to the equilibrium solver.
     */
    void addConstraint(const double *const p_A);
    
    /**
     * Removes all linear constraints from the equilibrium solver.
     */
    void clearConstraints();
    
    /**
     * Returns the partial derivative of the equilibrium mole fractions with
     * respect to temperature.
     */
    void dXdT(double* const p_dxdt) const;
    
    /**
     * Returns the partial derivative of the equilibrium mole fractions with
     * respect to pressure.
     */
    void dXdP(double* const p_dxdp) const;
    
    /**
     * Returns the current element potentials as computed by the equilibrate
     * function.
     */
    void elementPotentials(double* const p_lambda) const {
        std::fill(p_lambda, p_lambda+m_nc, 0.0);
        for (int i = 0; i < m_ncr; ++i) {
            p_lambda[m_cir(i)] = m_lambda(i);
        }
    }
    
    /**
     * Returns the current phase moles as computed by the equilibrate function.
     */
    void phaseMoles(double* const p_moles) const {
        for (int i = 0; i < m_np; ++i)
            p_moles[i] = m_Nbar(i);
    }


private:

    /**
     * Initializes the phase information by computing the total number of phases
     * present in the mixture (np) and assigning each species an integer from
     * 0 to np-1 corresponding to the phase to which the species belongs.
     */
    void initPhases();
    
    /**
     * Computes the partial derivatives dX/dalpha given dg/dalpha.  This method
     * is the common code used in dXdT() and dXdP().
     * @see dXdP()
     * @see dXdT()
     */
    void partialOfX(const Numerics::RealVector& dg, double* const p_dx) const;
    
    void initialConditions(
        const Numerics::RealVector& gtp, Numerics::RealVector& lambda, 
        Numerics::RealVector& Nbar, Numerics::RealVector& g);
    
    void perturb(Numerics::RealVector& nmm);
        
    std::pair<int, double> newton(
        Numerics::RealVector& lambda, Numerics::RealVector& Nbar, 
        const Numerics::RealVector& g, Numerics::RealVector& r,
        Numerics::RealSymMat& A, const double tol,
        const int max_iters = 10);
        
    /**
     * Computes the max-min composition \f$N^{mm}\f$, for the undetermined 
     * species which is defined as the composition that maximizes the minimum
     * species moles while satisfying the reduced constraints.
     *
     * @param c    reduced constraint vector
     * @param nmm  on return, the max-min composition for the constraints
     */
    void maxmin(const Numerics::RealVector& c, Numerics::RealVector& Nmm);
    
    /**
     * Computes the min-g composition, \f$N^g\f$, which is defined as the 
     * solution to
     *
     *     \f[ \min \sum_{k=1}^{n_{su}} N_k^u \tilde{g}_k^u \f]
     *
     * subject to the constraints \f$ \hat{B}^T N^u = \hat{c} \f$ and 
     * \f$ N_k^u \ge 0 \f$.  \f$ \tilde{g}_k^u \f$ is the normalized molar
     * specific Gibbs function 
     *
     *     \f[ \tilde{g}_k^u = 
     *            \frac{h_k}{R_u T} - \frac{S^\circ_k}{R_u} + 
     *            \ln \frac{p}{p_{atm}} \f]
     *
     * and \f$ \hat{B} \f$ and \f$ \hat{c} \f$ are the reduced constraint matrix
     * and reduced constraint vector respectively.
     *
     * @param T   temperature in K
     * @param P   pressure in Pa
     * @param cr  reduced constraint vector
     * @param gu  on return, vector of normalized molar specific Gibbs function 
     *            values for undetermined species
     * @param ng  on return, equals the min-g composition of undetermined
     *            species
     */
    void ming(
        const Numerics::RealVector& c, const Numerics::RealVector& gu, 
        Numerics::RealVector& Nmg);   
    
    void formSystemMatrix(Numerics::RealSymMat& A) const;
    
    bool updateSpeciesMoles(
        const Numerics::RealVector& lambda, const Numerics::RealVector& Nbar,
        const Numerics::RealVector& g);
        
    void computeResidual(
        const Numerics::RealVector& Nbar, Numerics::RealVector& r) const;
    
    //void updatePreConditioner();
    
    void reduceZeroSpecies();
     
private:

    const Thermodynamics& m_thermo;
    
    int m_ns;
    int m_ne;
    int m_nc;
    int m_np;
    
    int m_nsr;
    int m_ncr;
    
    double m_T;
    double m_P;

    Numerics::RealMatrix  m_B;
    Numerics::RealMatrix  m_Br;
    //Numerics::RealMatrix  m_P;
    
    Numerics::Vector<int> m_phase;
    Numerics::Vector<int> m_base;
    Numerics::Vector<int> m_sjr;
    Numerics::Vector<int> m_cir;
    
    Numerics::Vector<double> m_c;
    Numerics::Vector<double> m_cr;
    Numerics::Vector<double> m_N;
    Numerics::Vector<double> m_lambda;
    Numerics::Vector<double> m_Nbar;
    
    std::vector<Numerics::RealVector> m_constraints;

};

    } // namespace Thermodynamics
} // namespace Mutation

# endif // THERMO_MULTI_PHASE_EQUIL_SOLVER
