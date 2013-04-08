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
     * Computes the equilibrium composition of the mixture.
     */
    std::pair<int,int> equilibrate(
        double T, double P, const double *const p_ev, double *const p_sv);

private:

    /**
     * Initializes the phase information by computing the total number of phases
     * present in the mixture (np) and assigning each species an integer from
     * 0 to np-1 corresponding to the phase to which the species belongs.
     */
    void initPhases();
    
    void initialConditions(
        const Numerics::RealVector& c, const Numerics::RealVector& gtp, 
        Numerics::RealVector& lambda, Numerics::RealVector& Nbar, 
        Numerics::RealVector& g) const;
        
    /**
     * Computes the max-min composition \f$N^{mm}\f$, for the undetermined 
     * species which is defined as the composition that maximizes the minimum
     * species moles while satisfying the reduced constraints.
     *
     * @param c    reduced constraint vector
     * @param nmm  on return, the max-min composition for the constraints
     */
    void maxmin(const Numerics::RealVector& c, Numerics::RealVector& Nmm) const;
    
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
        Numerics::RealVector& Nmg) const;   
    
    void formSystemMatrix(
        const Numerics::RealVector& N, Numerics::RealSymMat& A) const;
    
    void computeSpeciesMoles(
        const Numerics::RealVector& lambda, const Numerics::RealVector& Nbar,
        const Numerics::RealVector& g, Numerics::RealVector& N) const;
        
    void computeResidual(
        const Numerics::RealVector& c, const Numerics::RealVector& Nbar,
        const Numerics::RealVector& N, Numerics::RealVector& r) const;
     
private:

    const Thermodynamics& m_thermo;
    
    int m_ns;
    int m_ne;
    int m_nc;
    int m_np;

    Mutation::Numerics::RealMatrix  m_B;
    Mutation::Numerics::Vector<int> m_phase;

};

    } // namespace Thermodynamics
} // namespace Mutation

# endif // THERMO_MULTI_PHASE_EQUIL_SOLVER
