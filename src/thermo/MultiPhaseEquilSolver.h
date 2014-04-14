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
    

//==============================================================================

/**
 * Simple class used to represent a reduced constraint matrix by pointing to the
 * full one in order to avoid having to allocate and copy the reduced matrix
 * when taking the SVD.
 */
class ReducedMatrix : public Numerics::MatExpr<double, ReducedMatrix>
{
public:
    ReducedMatrix(
        const Numerics::RealMatrix& A,
        const int* const ir, const int* const jr, int nr, int nc)
        : m_A(A), mp_ir(ir), mp_jr(jr), m_nr(nr), m_nc(nc)
    { }
    
    size_t size() const { return m_nr*m_nc; }
    size_t rows() const { return m_nr; }
    size_t cols() const { return m_nc; }
    
    double operator()(size_t i, size_t j) const {
        return m_A(mp_ir[i], mp_jr[j]);
    }

private:
    const Numerics::RealMatrix& m_A;
    const int* const mp_ir;
    const int* const mp_jr;
    size_t m_nr;
    size_t m_nc;
};

//==============================================================================

/**
 * Represents the H matrix which is the reduced form of diag(y)*B.
 */
class HMatrix : public Numerics::MatExpr<double, HMatrix>
{
public:
    HMatrix(
        const Numerics::RealMatrix& A, const double* const p_y,
        const int* const ir, const int* const jr, int nr, int nc)
        : m_Br(A, ir, jr, nr, nc), mp_y(p_y)
    { }
    
    size_t size() const { return m_Br.size(); }
    size_t rows() const { return m_Br.rows(); }
    size_t cols() const { return m_Br.cols(); }
    
    double operator()(size_t i, size_t j) const {
        return m_Br(i,j)*mp_y[i];
    }

private:
    const ReducedMatrix m_Br;
    const double* const mp_y;
};



//==============================================================================


/**
 * @class MultiPhaseEquilSolver
 * @brief Implements the Mutli-phase Gibbs function continuation method.
 *
 * @internal
 * Remaining issues that need to be implemented or fixed:
 * \li Implement change of basis set which will improve accuracy and prevent
 *     singular Jacobian matrices in some cases.  Could also improve Newton
 *     convergence.
 * \li Figure out how to correctly handle multiple phases in such a way that 
 *     infinite loops do not occur due to adding/removing the same phase over
 *     over.
 */
class MultiPhaseEquilSolver
{
public:
    /**
     * Constructs the equilibrium solver.
     */
    MultiPhaseEquilSolver(
        const Thermodynamics& thermo, const bool pure_condensed = true);
    
    /**
     * Destructor.
     */
    ~MultiPhaseEquilSolver();
    
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
    void dNdT(double* const p_dxdt) const;
    
    /**
     * Returns the partial derivative of the equilibrium mole fractions with
     * respect to pressure.
     */
    void dXdP(double* const p_dxdp) const;
    
    /**
     * Computes the partial derivatives of the equilibrium species moles with
     * respect to some change in the species Gibbs energies.  For instance, to
     * compute dN/dT, supply the dg/dT vector.
     */
    void dNdg(const double* const p_dg, double* const p_dNdg);

    void dXdg(const double* const p_dg, double* const p_dXdg);

    void dSoldg(const double* const p_dg, Numerics::RealVector& dx);

    /**
     * Returns the current element potentials as computed by the equilibrate
     * function.
     */
    void elementPotentials(double* const p_lambda) const
    {
        for (int i = 0; i < m_solution.ncr(); ++i)
            p_lambda[m_solution.cir()[i]] = m_solution.lambda()[i];
        for (int i = m_solution.ncr(); i < m_nc; ++i)
            p_lambda[m_solution.cir()[i]] = 0.0;
    }
    
    /**
     * Returns the current phase moles as computed by the equilibrate function.
     */
    void phaseMoles(double* const p_moles) const {
//        for (int m = 0; m < m_npr; ++m)
//            p_moles[m] = std::exp(m_solution.lnNbar(m));
    }
    
    /*
     * Returns the current species moles vector as computed the by the last
     * call to equilibrate().
     */
    void speciesMoles(double* const p_moles) const
    {
        double temp;
        for (int j = 0; j < m_solution.nsr(); ++j) {
            temp = m_solution.y()[j];
            p_moles[m_solution.sjr()[j]] = temp*temp;
        }

        for (int j = m_solution.nsr(); j < m_ns; ++j)
            p_moles[m_solution.sjr()[j]] = 0.0;
    }

    /**
     * Returns the number of continuation steps on the last call to equilibrate.
     */
    int nSteps() const {
        return m_niters;
    }
    
    /**
     * Returns the total number of newton iterations on the last call to
     * equilibrate.
     */
    int nNewtons() const {
        return m_nnewts;
    }

    /**
     * Computes the partial derivatives dN/dalpha given dg/dalpha.  This method
     * is the common code used in dNdT() and dNdP().  Note that it is safe to
     * use the same vector to pass in dg/dalpha and receive dN/dalpha.
     * @see dXdP()
     * @see dXdT()
     */
    //void partialOfN(const double* const p_dg, double* const p_dN) const;


private:

    class Solution
    {
    public:
    
        /**
         * Empty constructor.  All sizes are zero and no allocations are made.
         */
        Solution(const Thermodynamics& thermo)
            : m_np(0), m_nc(0), m_ns(0), m_npr(0), m_ncr(0), m_nsr(0),
              m_dsize(0), m_isize(0), mp_ddata(NULL), mp_idata(NULL),
              m_thermo(thermo)
        { }
        
        /**
         * Destructor.
         */
        ~Solution()
        {
            if (mp_ddata != NULL) delete [] mp_ddata;
            if (mp_idata != NULL) delete [] mp_idata;
            
            mp_g      = NULL;
            mp_y      = NULL;
            mp_lnNbar = NULL;
            mp_lambda = NULL;
            
            mp_sizes = NULL;
            mp_sjr   = NULL;
            mp_cir   = NULL;
        }
        
        int npr() const { return m_npr; }
        int ncr() const { return m_ncr; }
        int nsr() const { return m_nsr; }
        
        /**
         * Returns the reduced species list as a set of indices into the global
         * species array.
         */
        const int* const sjr() const { return mp_sjr; }
        
        /**
         * Returns the reduced constraints list as a set of indices into the 
         * global constraints array.
         */
        const int* const cir() const { return mp_cir; }
        
        /**
         * Returns the sizing vector.
         */
        const int* const sizes() const { return mp_sizes; }
        
        /**
         * Returns the current y = sqrt(N) vector.
         */
        const double* const y() const { return mp_y; }
        
        /**
         * Returns the current log(Nbar) vector.
         */
        const double* const lnNbar() const { return mp_lnNbar; }
        
        /**
         * Returns the current g vector associated with this solution.
         */
        const double* const g() const { return mp_g; }
        
        /**
         * Returns the current lambda vector associated with this solution.
         */
        const double* const lambda() const { return mp_lambda; }
        
        /**
         * Returns a reduced constraint matrix based on the solution ordering.
         */
        const ReducedMatrix reducedMatrix(const Numerics::RealMatrix& B) const {
            return ReducedMatrix(B, mp_sjr, mp_cir, m_nsr, m_ncr);
        }
        
        /**
         * Returns a reduced H matrix which is the reduced form of diag(y)*B
         * the solution's ordering.
         */
        const HMatrix hMatrix(const Numerics::RealMatrix& B) const {
            return HMatrix(B, mp_y, mp_sjr, mp_cir, m_nsr, m_ncr);
        }
        
        /**
         * Initializes data storage for the solution.  Solution and ordering are
         * just filled with zeros.
         */
        void initialize(int np, int nc, int ns)
        {
            // Init sizes
            m_np = np;
            m_nc = nc;
            m_ns = ns;
            m_npr = np;
            m_ncr = nc;
            m_nsr = ns;
            
            int dsize = 2*ns+np+nc;
            int isize = ns+nc+np+2;
            
            // Allocate data storage if necessary
            if (dsize > m_dsize) {
                if (mp_ddata != NULL)
                    delete [] mp_ddata;
                mp_ddata = new double [dsize];
            }
            m_dsize = dsize;
            
            if (isize > m_isize) {
                if (mp_idata != NULL)
                    delete [] mp_idata;
                mp_idata = new int [isize];
            }
            m_isize = isize;
            
            // Setup the pointer locations
            mp_g      = mp_ddata;
            mp_y      = mp_g+ns;
            mp_lnNbar = mp_y+ns;
            mp_lambda = mp_lnNbar+np;
            
            mp_sizes = mp_idata;
            mp_sjr   = mp_sizes+np+2;
            mp_cir   = mp_sjr+ns;
            
            // Just fill all data with 0
            std::fill(mp_ddata, mp_ddata+m_dsize, 0.0);
            std::fill(mp_idata, mp_idata+m_isize, 0);
        }
        
        /**
         * Initializes the ordering based on defined species "groups" and zero
         * constraints.
         */
        void setupOrdering(int* species_group, bool* zero_constraint, bool pure)
        {
            // Count the number of species in each group and order the species such that
            // the groups are contiguous
            for (int i = 0; i < m_np+2; ++i)
                mp_sizes[i] = 0;

            for (int i = 0; i < m_ns; ++i)
                mp_sizes[species_group[i]+1]++;
            
            for (int i = 1; i < m_np+1; ++i)
                mp_sizes[i] += mp_sizes[i-1];
            
            for (int i = 0; i < m_ns; ++i)
                mp_sjr[mp_sizes[species_group[i]]++] = i;
            
            for (int i = m_np+1; i > 0; --i)
                mp_sizes[i] = mp_sizes[i-1];
            mp_sizes[0] = 0;
            
            m_nsr = mp_sizes[m_np];
            
            // Update list of reduced constraint indices
            int index = 0;
            for (int i = 0; i < m_nc; ++i)
                if (!zero_constraint[i]) mp_cir[index++] = i;
            
            m_ncr = index;
            m_npr = m_np;
            int i = m_np-1;
            while (i > 0) {
                if (mp_sizes[i+1] - mp_sizes[i] == 0)
                    removePhase(i);
                i--;
            }
        }
        
        /**
         * Updates the solution by adding dx to lambda and dlnNbar and
         * recomputes the y vector.  Note that g remains fixed.
         */
        template <typename E>
        void update(
            const Numerics::VecExpr<double, E>& dx,
            const Numerics::RealMatrix& B)
        {
            for (int i = 0; i < m_ncr; ++i)
                mp_lambda[i] += dx(i);
            for (int m = 0; m < m_npr; ++m)
                mp_lnNbar[m] = std::min(300.0, mp_lnNbar[m]+dx(m+m_ncr));
                
            updateY(B);
        }
        
        /**
         * Sets the g vector associated with the solution.  The solution is not
         * updated however.
         */
        void setG(const double* const p_g0, const double* const p_g, double s)
        {
            const double s1 = 1.0 - s;
            for (int j = 0; j < m_ns; ++j)
                mp_g[j] = s1*p_g0[j] + s*p_g[j];
        }
        
        /**
         * Sets the solution vector.
         */
        void setSolution(
            const double* const p_lambda, const double* const p_Nbar,
            const Numerics::RealMatrix& B)
        {
            for (int i = 0; i < m_ncr; ++i)
                mp_lambda[i] = p_lambda[i];
            for (int m = 0; m < m_npr; ++m)
                mp_lnNbar[m] = std::log(p_Nbar[m]);
            
            updateY(B);
        }
        
        /**
         * Updates the y vector based on the current solution variables and g
         * vector.
         */
        void updateY(const Numerics::RealMatrix& B)
        {
            int jk;
            for (int m = 0; m < m_npr; ++m) {
                for (int j = mp_sizes[m]; j < mp_sizes[m+1]; ++j) {
                    jk = mp_sjr[j];
                    mp_y[j] = mp_lnNbar[m] - mp_g[jk];
                    for (int i = 0; i < m_ncr; ++i)
                        mp_y[j] += B(jk, mp_cir[i])*mp_lambda[i];
                    mp_y[j] = std::exp(0.5*std::min(mp_y[j], 300.0));
                }
            }
        }
        
        /**
         * Equality operator.  Reallocates only if the data size required is 
         * larger than the current capacity.
         */
        Solution& operator=(const Solution& state)
        {
            m_np  = state.m_np;
            m_ns  = state.m_ns;
            m_npr = state.m_npr;
            m_ncr = state.m_ncr;
            m_nsr = state.m_nsr;
            
            if (state.m_dsize > m_dsize) {
                if (mp_ddata != NULL) delete mp_ddata;
                mp_ddata  = new double [state.m_dsize];
                mp_g      = mp_ddata;
                mp_y      = mp_g+m_ns;
                mp_lnNbar = mp_y+m_ns;
                mp_lambda = mp_lnNbar+m_np;
            }
            m_dsize = state.m_dsize;
            
            if (state.m_isize > m_isize) {
                if (mp_idata != NULL) delete mp_idata;
                mp_idata = new int [state.m_isize];
                mp_sizes = mp_idata;
                mp_sjr   = mp_sizes+m_np+2;
                mp_cir   = mp_sjr+m_ns;
            }
            m_isize = state.m_isize;
            
            std::copy(state.mp_ddata, state.mp_ddata+m_dsize, mp_ddata);
            std::copy(state.mp_idata, state.mp_idata+m_isize, mp_idata);
            
            return *this;
        }
        
        /**
         * Removes the phase from the solution and returns the new phase number
         * associated with this phase.  The species ordering, and solution 
         * variables are all updated to correspond to the new solution.
         */
        int removePhase(int phase)
        {
//            cout << "removing phase:";
//            for (int j = mp_sizes[phase]; j < mp_sizes[phase+1]; ++j)
//                cout << " " << m_thermo.speciesName(mp_sjr[j]);
//            cout << endl;
            
            // Check that the phase number is feasible
            assert(phase >= 0);
            assert(phase < m_npr);
            assert(m_npr > 1);
            
            // If the phase is not the last non-empty phase, we need to shift it
            // to the right to make the non-empty species list contiguous
            const int size = mp_sizes[phase+1] - mp_sizes[phase];
            if (phase != m_npr-1) {
                int temp [size];
                
                // First copy the phase to be removed into temporary array
                int* ip = temp;
                for (int j = mp_sizes[phase]; j < mp_sizes[phase+1]; ++j)
                    *ip++ = mp_sjr[j];
                
                // Shift all remaining phases to fill the space
                ip = mp_sjr + mp_sizes[phase];
                for (int j = mp_sizes[phase+1]; j < mp_sizes[m_np]; ++j)
                    *ip++ = mp_sjr[j];
                
                // Place the removed phase at the end of the list (before
                // determined species)
                for (int m = phase+1; m < m_np; ++m)
                    mp_sizes[m] = mp_sizes[m+1] - size;
                ip = temp;
                for (int j = mp_sizes[m_np-1]; j < mp_sizes[m_np]; ++j)
                    mp_sjr[j] = *ip++;
                
                // Shift the solution to match the new phase ordering
                for (int j = mp_sizes[phase]; j < m_nsr-size; ++j)
                    mp_y[j] = mp_y[j+size];
                for (int m = phase; m < m_npr-1; ++m)
                    mp_lnNbar[m] = mp_lnNbar[m+1];
                
                // Set the phase number to the new number to return
                phase = m_np-1;
            }
            
            // Update reduced species and phase sizes
            m_npr--;
            m_nsr -= size;
            
            return phase;
        }
        
        /**
         * Moves an inactive phase to the active list and returns the new phase
         * number for that phase.  Note that the solution is not updated.
         */
        int addPhase(int phase)
        {
//            cout << "adding phase:";
//            for (int j = mp_sizes[phase]; j < mp_sizes[phase+1]; ++j)
//                cout << " " << m_thermo.speciesName(mp_sjr[j]);
//            cout << endl;
        
            assert(phase >= m_npr);
            assert(phase < m_np);
            
            // If the phase is not the first empty phase, then we need to shift
            // it to the front
            int size = mp_sizes[phase+1] - mp_sizes[phase];
            if (phase > m_npr) {
                int temp [size];
                
                // First copy the phase to be added into temporary array
                int* p = temp;
                for (int i = mp_sizes[phase]; i < mp_sizes[phase+1]; ++i)
                    *p++ = mp_sjr[i];
                
                // Shift phases that come before to the back
                p = mp_sjr + mp_sizes[phase+1]-1;
                for (int i = mp_sizes[phase]-1; i >= m_nsr; --i)
                    *p-- = mp_sjr[i];
                
                // Place the added phase at the end of the included list
                for (int i = phase+1; i > m_npr; --i)
                    mp_sizes[i] = mp_sizes[i-1] + size;
                p = temp;
                for (int i = mp_sizes[m_npr]; i < mp_sizes[m_npr+1]; ++i)
                    mp_sjr[i] = *p++;
            }
            
            // Update reduced species and phase sizes
            m_npr++;
            m_nsr += size;
            
            return (m_npr-1);
        }
        
        /**
         * Returns the index of the phase which will reduce the Gibbs energy the
         * most with the current solution.  If no phase will reduce further than
         * -1 is returned.
         */
        int checkCondensedPhase(const Numerics::RealMatrix& B)
        {
            if (m_np <= m_ncr)
                return -1;

            int    min_m   = -1;
            double min_sum = 0.0;

            for (int m = m_npr; m < m_np; ++m) {
                for (int j = mp_sizes[m]; j < mp_sizes[m+1]; ++j) {
                    double sum = mp_g[mp_sjr[j]];
                    for (int i = 0; i < m_ncr; ++i)
                        sum -= mp_lambda[i]*B(mp_sjr[j], mp_cir[i]);

                    if (sum < min_sum) {
                        min_m   = m;
                        min_sum = sum;
                    }
                }
            }

            return min_m;
        }

        void printOrder()
        {
            cout << "Species order:" << endl;
            cout << "  Active Phases:" << endl;
            for (int m = 0; m < m_npr; ++m) {
                cout << "    " << m << ":";
                for (int j = mp_sizes[m]; j < mp_sizes[m+1]; ++j)
                    cout << " " << m_thermo.speciesName(mp_sjr[j]);
                cout << endl;
            }
            cout << "  Inactive Phases:" << endl;
            for (int m = m_npr; m < m_np; ++m) {
                cout << "    " << m << ":";
                for (int j = mp_sizes[m]; j < mp_sizes[m+1]; ++j)
                    cout << " " << m_thermo.speciesName(mp_sjr[j]);
                cout << endl;
            }
            cout << "  Determined Species:" << endl;
            cout << "   ";
            for (int j = mp_sizes[m_np]; j < m_ns; ++j)
                cout << " " << m_thermo.speciesName(mp_sjr[j]);
            cout << endl;
        }
        
        void printSolution()
        {
            cout << "Solution:" << endl;
            cout << "  lambda = " << endl;
            for (int i = 0; i < m_ncr; ++i)
                cout << "    " << mp_lambda[i] << endl;
            cout << "  Nbar = " << endl;
            for (int m = 0; m < m_npr; ++m)
                cout << "    " << std::exp(mp_lnNbar[m]) << endl;
            cout << "  N = " << endl;
            for (int j = 0; j < m_nsr; ++j)
                cout << "   " << setw(12) << m_thermo.speciesName(mp_sjr[j])
                     << ": " << mp_y[j]*mp_y[j] << endl;
        }
        
        void printG()
        {
            cout << "Current G vector" << endl;
            for (int j = 0; j < m_ns; ++j)
                cout << setw(12) << m_thermo.speciesName(j) << ": "
                     << mp_g[j] << endl;
        }
        
    private:
    
        int m_np;
        int m_nc;
        int m_ns;
        int m_npr;
        int m_ncr;
        int m_nsr;
        int m_dsize;
        int m_isize;
        
        double* mp_ddata;
        int*    mp_idata;
        
        double* mp_g;
        double* mp_y;
        double* mp_lnNbar;
        double* mp_lambda;
        
        int* mp_sizes;
        int* mp_sjr;
        int* mp_cir;
        
        const Thermodynamics& m_thermo;
    };

    /**
     * Initializes the phase information by computing the total number of phases
     * present in the mixture (np) and assigning each species an integer from
     * 0 to np-1 corresponding to the phase to which the species belongs.
     */
    void initPhases();
    
    /**
     * Figure out "determined species" and setup the species ordering data
     * structure.
     */
    void checkForDeterminedSpecies();
    
    void initialConditions();
    void rates(Numerics::RealVector& dx);
    double newton();
    void updateMinGSolution(const double* const p_g);
    void updateMaxMinSolution();
    void formSystemMatrix(Numerics::RealSymMat& A) const;
    void computeResidual(Numerics::RealVector& r) const;

    /**
     * Adjusts the phase ordering according to the vapor pressure test and phase
     * rule based on the current (converged) solution parameters.
     *
     * @return True if the phase ordering changed, false otherwise
     */
    bool phaseRedistribution();
     
private:

    static const double ms_eps_rel;
    static const double ms_eps_abs;
    static const double ms_ds_inc;
    static const double ms_ds_dec;
    static const double ms_max_ds;

    const Thermodynamics& m_thermo;
    
    int m_ns;
    int m_ne;
    int m_nc;
    int m_np;
    
    double m_T;
    double m_P;

    Numerics::RealMatrix  m_B;
    
//    std::vector<Numerics::RealVector> m_constraints;
    
    // ** new data **
    Solution m_solution;
    
    size_t  m_tableau_capacity;
    double* mp_tableau;
    double* mp_ming;
    double* mp_maxmin;
    
    int* mp_phase;
    double* mp_g;
    double* mp_g0;
    double* mp_c;
    
    bool m_pure_condensed;
    
    int m_niters;
    int m_nnewts;

};

    } // namespace Thermodynamics
} // namespace Mutation

# endif // THERMO_MULTI_PHASE_EQUIL_SOLVER
