#ifndef THERMO_GFC_EQUIL_SOLVER
#define THERMO_GFC_EQUIL_SOLVER

#include "Thermodynamics.h"
#include "Numerics.h"

#include <utility>

/**
 * Equilibrium solver which uses the Gibbs function continuation method by
 * S. B. Pope.  This solver is guaranteed to work for all well posed problems
 * unlike the basic equilibrium solver.
 */
class GfcEquilSolver
{
public:

    /**
     * Constructs the equilibrium solver.
     */
    GfcEquilSolver(const Thermodynamics& thermo);
    
    /**
     * Destructor.
     */
    ~GfcEquilSolver() {
        if (mp_ls_Br != NULL) {
            delete mp_ls_Br;
            mp_ls_Br = NULL;
        }
    }
    
    /**
     * Computes the equilibrium composition of the mixture.
     */
    std::pair<int,int> equilibrate(
        const double &T, const double &P, const double *const p_ev, 
        double *const p_sv);

private:

    bool newton(
        const Numerics::RealVector& cr, const Numerics::RealVector& gu, 
        const Numerics::RealVector& lam0, const Numerics::RealVector& Q, 
        Numerics::RealVector& lam, Numerics::RealVector& y, double& error,
        std::pair<int,int>& stats) const;

    void perturb(Numerics::RealVector& cr, Numerics::RealVector& nmm) const;
    
    void lamgInit(
        const Numerics::RealVector& zd, const Numerics::RealVector& zu0, 
        const Numerics::RealVector& gu, Numerics::RealVector& lam0, 
        Numerics::RealVector& gu0) const;
    
    /**
     * Computes the max-min composition \f$N^{mm}\f$, for the undetermined 
     * species which is defined as the composition that maximizes the minimum
     * species moles while satisfying the reduced constraints.
     *
     * @param c    reduced constraint vector
     * @param nmm  on return, the max-min composition for the constraints
     */
    void maxmin(const Numerics::RealVector& c, Numerics::RealVector& nmm) const;
    
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
        const double& T, const double& P, const Numerics::RealVector& cr, 
        Numerics::RealVector& gu, Numerics::RealVector& ng) const;
        
    void solveFixedT(
        const Numerics::RealVector& gu, const Numerics::RealVector& lam0, 
        const Numerics::RealVector& gu0, const double& res0, 
        const Numerics::RealVector& Q, const Numerics::RealVector& cr,
        Numerics::RealVector& zu, std::pair<int,int>& stats) const;
    
    bool computeRate(
        const Numerics::RealVector& lam, const Numerics::RealVector& gu, 
        const Numerics::RealVector& dgudt, const Numerics::RealVector& Q, 
        Numerics::RealVector& dlamdt) const;
    
    /**
     * Computes the vector \f$ y=\sqrt{\exp(-\bar{g}+\hat{B}\bar{\lambda})} \f$
     * by avoiding the square-root and also checks for possible overflow.
     * Returns false if overflow will occur, true on success.
     */
    bool formy(
        const Numerics::RealVector& gu, const Numerics::RealVector& lam,
        Numerics::RealVector& y) const;
private:
    
    const Thermodynamics& m_thermo;
    
    int m_ns;  ///< number of species
    int m_nsd; ///< number of determined species
    int m_nsu; ///< number of undetermined species
    int m_ne;  ///< number of elements
    int m_ned; ///< number of determined elements
    int m_neu; ///< number of undetermined elements
    int m_nb;  ///< number of linearly independent constraints
    int m_nrc; ///< number of reduced constraints
    
    std::vector<int> m_species_order;
    std::vector<int> m_element_order;
    
    Numerics::RealMatrix m_B;    ///< Basic constraint matrix
    Numerics::RealMatrix m_B_r;  ///< Reduced constraint matrix (RCM)
    Numerics::RealMatrix m_B_rt; ///< Transpose of reduced constraint matrix
    Numerics::RealMatrix m_A;    ///< Modified constraint xform matrix
    
    Numerics::RealVector m_au;   ///< Number of elements in molecule of und. sp.
    
    Numerics::LeastSquares<double>* mp_ls_Br; ///< Least-squares of RCM

private:

    static const double sm_max_min_frac;
    static const double sm_logy_lim;
    static const double sm_res_tol;
    static const double sm_ds_inc;
    static const double sm_ds_dec;
    static const double sm_dec_min;
    static const double sm_err_huge;
    static const double sm_eps_e;
    static const double sm_eps_s;
    
};

#endif // THERMO_GFC_EQUIL_SOLVER
