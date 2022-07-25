/**
 * @file MultiPhaseEquilSolver.h
 *
 * @brief Defines the MultiPhaseEquilSolver class which implements the
 * multiphase Gibbs function continuation (MPGFC) method.
 */

/*
 * Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Mutation++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Mutation++.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef THERMO_MULTI_PHASE_EQUIL_SOLVER
#define THERMO_MULTI_PHASE_EQUIL_SOLVER

#include <vector>
#include <Eigen/Dense>



namespace Mutation {
    namespace Thermodynamics {

class Thermodynamics;

    
/**
 * Enumerates the two possibilities for defining mole fractions.
 */
enum MoleFracDef {
    IN_PHASE,
    GLOBAL
};


//==============================================================================

/**
 * Simple class used to represent a reduced constraint matrix by pointing to the
 * full one in order to avoid having to allocate and copy the reduced matrix
 * when taking the SVD.
 */
//class ReducedMatrix : public Eigen::MatrixBase<ReducedMatrix>
//{
//public:
//    ReducedMatrix(
//        const Numerics::RealMatrix& A,
//        const int* const ir, const int* const jr, int nr, int nc)
//        : m_A(A), mp_ir(ir), mp_jr(jr), m_nr(nr), m_nc(nc)
//    { }
//
//    size_t size() const { return m_nr*m_nc; }
//    size_t rows() const { return m_nr; }
//    size_t cols() const { return m_nc; }
//
//    double operator()(size_t i, size_t j) const {
//        return m_A(mp_ir[i], mp_jr[j]);
//    }
//
//private:
//    const Numerics::RealMatrix& m_A;
//    const int* const mp_ir;
//    const int* const mp_jr;
//    size_t m_nr;
//    size_t m_nc;
//};

//==============================================================================

/**
 * Represents the H matrix which is the reduced form of diag(y)*B.
 */
//class HMatrix : public Numerics::MatExpr<double, HMatrix>
//{
//public:
//    HMatrix(
//        const Numerics::RealMatrix& A, const double* const p_y,
//        const int* const ir, const int* const jr, int nr, int nc)
//        : m_Br(A, ir, jr, nr, nc), mp_y(p_y)
//    { }
//
//    size_t size() const { return m_Br.size(); }
//    size_t rows() const { return m_Br.rows(); }
//    size_t cols() const { return m_Br.cols(); }
//
//    double operator()(size_t i, size_t j) const {
//        return m_Br(i,j)*mp_y[i];
//    }
//
//private:
//    const ReducedMatrix m_Br;
//    const double* const mp_y;
//};



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
     * Computes the equilibrium composition of the mixture.  The returned
     * species mole-fraction array is computed based on the value of the
     * MoleFracDef parameter.  The default method defines the mole fractions
     * globally, though they can also be defined in-phase.
     *
     * @return A pair representing the number of iterations and Newton
     * iterations used to compute the solution.  If the solution procedure fails
     * because the feasible region is empty or unbounded, then the iterations
     * are negative.
     */
    std::pair<int,int> equilibrate(
        double T, double P, const double *const p_ev, double *const p_sv,
        MoleFracDef mfd = GLOBAL);
    
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
    void dNdg(const double* const p_dg, double* const p_dNdg) const;

    void dXdg(const double* const p_dg, double* const p_dXdg) const;

    void dXdc(int i, double* const p_dxdc);

    void dSoldg(const double* const p_dg, Eigen::VectorXd& dx) const;

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
        for (int m = 0; m < m_np; ++m)
            p_moles[m] = 0.0;
        for (int mr = 0; mr < m_solution.npr(); ++mr) {
            int m = mp_phase[m_solution.sjr()[m_solution.sizes()[mr]]];
            p_moles[m] = std::exp(m_solution.lnNbar()[mr]);
        }
    }
    
    /*
     * Returns the current species moles vector as computed by the last
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
        //const ReducedMatrix reducedMatrix(const Numerics::RealMatrix& B) const {
        //    return ReducedMatrix(B, mp_sjr, mp_cir, m_nsr, m_ncr);
        //}
        Eigen::MatrixXd& reducedMatrix(
            const Eigen::MatrixXd& B, Eigen::MatrixXd& Br) const
        {
            Br.resize(m_nsr, m_ncr);
            for (int j = 0; j < m_nsr; ++j)
                for (int i = 0; i < m_ncr; ++i)
                    Br(j,i) = B(mp_sjr[j], mp_cir[i]);
            return Br;
        }
        
        /**
         * Returns a reduced H matrix which is the reduced form of diag(y)*B
         * the solution's ordering.
         */
        //const HMatrix hMatrix(const Numerics::RealMatrix& B) const {
        //    return HMatrix(B, mp_y, mp_sjr, mp_cir, m_nsr, m_ncr);
        //}
        
        /**
         * Initializes data storage for the solution.  Solution and ordering are
         * just filled with zeros.
         */
        void initialize(int np, int nc, int ns);
        
        /**
         * Initializes the ordering based on defined species "groups" and zero
         * constraints.
         * @retval true if internal ordering information changed after this call
         */
        bool setupOrdering(int* species_group, bool* zero_constraint);
        
        /**
         * Updates the solution by adding dx to lambda and dlnNbar and
         * recomputes the y vector.  Note that g remains fixed.
         */
        template <typename E>
        void update(
            const Eigen::MatrixBase<E>& dx, const Eigen::MatrixXd& B)
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
        void setG(const double* const p_g0, const double* const p_g, double s);
        
        /**
         * Sets the solution vector.
         */
        void setSolution(
            const double* const p_lambda, const double* const p_Nbar,
            const Eigen::MatrixXd& B);
        
        /**
         * Updates the y vector based on the current solution variables and g
         * vector.
         */
        void updateY(const Eigen::MatrixXd& B);
        
        /**
         * Equality operator.  Reallocates only if the data size required is 
         * larger than the current capacity.
         */
        Solution& operator=(const Solution& state);
        
        /**
         * Removes the phase from the solution and returns the new phase number
         * associated with this phase.  The species ordering, and solution 
         * variables are all updated to correspond to the new solution.
         */
        int removePhase(int phase);
        
        /**
         * Moves an inactive phase to the active list and returns the new phase
         * number for that phase.  Note that the solution is not updated.
         */
        int addPhase(int phase);
        
        /**
         * Returns the index of the phase which will reduce the Gibbs energy the
         * most with the current solution.  If no phase will reduce further than
         * -1 is returned.
         */
        int checkCondensedPhase(const Eigen::MatrixXd& B);

        void printOrder();
        
        void printSolution();
        
        void printG();
        
        /**
         * Unpacks the current solution into an output array of mole fractions.
         */
        void unpackMoleFractions(double* const p_x, const MoleFracDef mfd);

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
    bool checkForDeterminedSpecies();
    
    /**
     * Sets up the initial conditions of the equilibrium problem.
     * @return true if a set of initial conditions could be computed, false
     * otherwise
     */
    bool initialConditions(
        const double T, const double P, const double* const p_c);

    /**
     * Computes \f$\bar{N}(0)\f$, \f$\lambda(0)\f$, and \f$\tilde{g}(0)\f$ given
     * an initial set of species moles.
     */
    void initZeroResidualSolution(
    	double* const p_N, double* const p_Nbar, double* const p_lambda);

    /**
     * Computes the solution derivative with respect to \f$s\f$ at constant
     * residual.
     */
    void rates(Eigen::VectorXd& dx, bool save = false);

    void rates_ci(int i, Eigen::VectorXd& dx);

    /**
     * Tries to reduce the solution residual (at a fixed value of \f$s\f$) below
     * a given tolerance using Newton's method.
     */
    double newton();

    /**
     * Updates the Min-G solution based on the current temperature and
     * constraints.
     * @return true if Min-G solution exists, false otherwise
     */
    bool updateMinGSolution(const double* const p_g);

    /**
     * Updates the Max-Min solution based on the current constraints.
     * @return true if Max-Min solution exists, false otherwise
     */
    bool updateMaxMinSolution();

    /**
     * Computes the Jacobian of hte system for the Newton's method.
     */
    void formSystemMatrix(Eigen::MatrixXd& A) const;

    /**
     * Computes the current value of the residual vector.
     */
    void computeResidual(Eigen::VectorXd& r) const;

    /**
     * Adjusts the phase ordering according to the vapor pressure test and phase
     * rule based on the current (converged) solution parameters.  If a phase
     * is added, then the solution is reinitialized based after adding some
     * small amount of moles of the new phase.
     *
     * @return true if a phase was added, false otherwise
     */
    bool phaseRedistribution();
     
private:

    static const double ms_eps_rel;
    static const double ms_eps_abs;
    static const double ms_ds_inc;
    static const double ms_ds_dec;
    static const double ms_max_ds;
    static const double ms_temp_change_tol;
    static const double ms_pres_change_tol;

    const Thermodynamics& m_thermo;
    
    int m_ns;
    int m_ne;
    int m_nc;
    int m_np;
    
    double m_T;
    double m_P;

    //Numerics::RealMatrix m_B;
    Eigen::MatrixXd m_B;
    Eigen::MatrixXd m_Br;
    
    std::vector<Eigen::VectorXd> m_constraints;
    
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
