#include "Numerics.h"
#include "Thermodynamics.h"

namespace Mutation {
    namespace Thermodynamics {

class BPrimeProvider
{
public:
    BPrimeProvider(const Thermodynamics& thermo)
        : m_thermo(thermo), m_g(thermo.nSpecies()), 
          m_lambda(thermo.nElements()), m_delta(thermo.nElements()+2),
          m_res(thermo.nElements()+2), m_Ykc(thermo.nElements()), 
          m_Et(m_thermo.elementMatrix().transpose()),
          m_J(thermo.nElements()+2,thermo.nElements()+2),
          m_mwk(thermo.nElements())
    {
        m_ce_index = m_thermo.elementIndex("C");
        m_cs_index = m_thermo.speciesIndex("C");
        
        if (m_ce_index < 0 || m_cs_index < 0) {
            cout << "Error: B' tables can only be made for mixtures "
                 << "containing Carbon!";
            std::exit(1);
        }
        
        m_Ykc = 0.0;
        m_Ykc(m_ce_index) = 1.0;
        
        for (int k = 0; k < m_thermo.nElements(); k++)
            m_mwk(k) = m_thermo.atomicMass(k);
    }
    
    /**
     * Solves the coupled equilibrium mass balance equations for the 
     * non-dimensional ablation mass flux at the surface of an ablator.
     * 
     * @param Yke - elemental mass fractions at boundary layer edge
     * @param Ykg - elemental mass fractions of pyrolysis gas at the surface
     * @param Bg  - non-dimensional pyrolysis gas mass flux
     * @param T   - temperature at the surface
     * @param P   - pressure at the surface
     * @param Bc  - on return, non-dimensional ablation mass flux
     * @param X   - on return, the species mole fractions at the surface
     */
    void massBalance(
        const double* const Yke, const double* const Ykg, double Bg, 
        const double T, const double P, double& Bc, double* const p_X)
    {
        const int ne = m_thermo.nElements();
        const int ns = m_thermo.nSpecies();
        const double tol = 1.0E-8;
        
        // Element matrix E_ij = number of atoms of element j in species i
        const RealMatrix& E = m_thermo.elementMatrix();
        
        // Write out the problem
        cout << "Mass Balance:" << endl;
        cout << "Yke:" << endl;
        for (int i = 0; i < ne; i++)
            cout << Yke[i] << endl;
        cout << "Ykg:" << endl;
        for (int i = 0; i < ne; i++)
            cout << Ykg[i] << endl;
        cout << "Ykc:" << endl;
        cout << m_Ykc << endl;
        cout << "T:  " << T << endl;
        cout << "P:  " << P << endl;
        cout << "Bg: " << Bg << endl;
        cout << endl;
        
        // Compute species Gibbs free energies
        m_thermo.speciesGOverRT(T, m_g);
        m_lambda = 1.0;
        m_Bc     = 0.0;
        m_M      = 1.0;
        m_Kcgr   = 1.0 / 0.0001; // for now
        computeX(ns, p_X);
        
        // Compute the residual
        double resini = residual(Yke, Ykg, Bg, p_X);
        double res = resini;
        double alpha = 0.5;
        
        cout << "Initial residual: " << resini << endl;
        
        // Newton method
        int iter = 0;
        while (res > tol && iter < 10000) {
            cout << "iter: " << iter << endl;
            cout << "res:  " << res << endl;
            cout << "lambda: " << endl;
            cout << m_lambda << endl;
            cout << "M: " << m_M << endl;
            cout << "Bc: " << m_Bc << endl;
            // Compute the Jacobian
            jacobian(Yke, Ykg, Bg, p_X);
            cout << "J:" << endl;
            cout << m_J << endl;
            cout << "R: " << endl;
            cout << m_res << endl;
            cout << "X" << endl;
            for (int i = 0; i < ns; i++)
                cout << m_thermo.speciesName(i) << " " << p_X[i] << endl;
            
            // Solve linear system for solution update
            Numerics::SVD<double> svd(m_J);
            svd.solve(m_delta, m_res);
            
            // Update the solution
            m_lambda -= alpha * m_delta(0,ne);
            m_M      -= alpha * m_delta(ne);
            m_Bc     -= alpha * m_delta(ne+1);
            
            // Now recompute X and the residual
            computeX(ns, p_X);
            res = residual(Yke, Ykg, Bg, p_X);
            iter++;
        }
        
        // Hopefully we are converged at this point
        Bc = m_Bc;        
    }
    
private:

    /**
     * Computes the species mole fractions from the current values of g and
     * lambda.
     */
    void computeX(const int ns, double* const p_X) const
    {
        Numerics::VectorWrapper<double> X(p_X, ns);
        X = Numerics::exp(m_thermo.elementMatrix() * m_lambda - m_g);        
    }
    
    /**
     * Updates the residual.
     */
    double residual(
        const double* const Yke, const double* const Ykg, double Bg, 
        double* const p_X)
    {
        const int ne = m_thermo.nElements();
        const int ns = m_thermo.nSpecies();
        
        Numerics::VectorWrapper<double> X(p_X, ns);
        
        m_res(0,ne) = m_Et * X;
        double factor = m_M / (1.0 + Bg + m_Bc);
        for (int i = 0; i < ne; i++)
            m_res(i) -= 
                factor / m_mwk(i) * (Yke[i] + Bg*Ykg[i] + m_Bc*m_Ykc(i));
        
        m_res(ne) = X.sum() - 1.0;
        m_res(ne+1) = X(m_cs_index) - 1.0 / m_Kcgr;
        
        return m_res.norm2();
    }
    
    /**
     * Updates the system Jacobian.
     */
    void jacobian(
        const double* const Yke, const double* const Ykg, double Bg, 
        double* const p_X)
    {
        const int ne = m_thermo.nElements();
        const int ns = m_thermo.nSpecies();
        
        Numerics::VectorWrapper<double> X(p_X, ns);
        
        // Compute d R1_k/d lambda_j
        for (int k = 0; k < ne; k++) {
            for (int j = 0; j < ne; j++) {
                m_J(k,j) = 0.0;
                for (int i = 0; i < ns; i++)
                    m_J(k,j) += m_Et(k,i) * X(i) * m_Et(j,i);
            }
        }
        
        // d R1_k /d M and d R1_k /d Bc
        double fac = 1.0 / (1.0 + Bg + m_Bc);
        for (int k = 0; k < ne; k++) {
            m_J(k,ne) = 
                -fac / m_mwk(k) * (Yke[k] + Bg*Ykg[k] + m_Bc*m_Ykc(k));
            m_J(k,ne+1) = 
                -m_M * fac * fac / m_mwk(k) * 
                ((m_Ykc(k) - Yke[k]) + Bg*(m_Ykc(k) - Ykg[k]));
        }
        
        // R2 derivatives
        m_J(ne,ne+1,0,ne) = m_Et * X;
        m_J(ne,ne+1,ne,ne+2) = 0.0;
        
        // R3 derivatives
        m_J(ne+1,ne+2,0,ne+2) = 0.0;
        m_J(ne+1,m_ce_index) = X(m_cs_index);
        
        cout << "J" << endl;
        cout << m_J << endl;
    }

private:

    const Thermodynamics& m_thermo;
    Numerics::RealVector m_g;
    Numerics::RealVector m_lambda;
    Numerics::RealVector m_delta;
    Numerics::RealVector m_res;
    Numerics::RealVector m_Ykc;
    Numerics::RealMatrix m_Et;
    Numerics::RealMatrix m_J;
    Numerics::RealVector m_mwk;
    
    double m_M;
    double m_Bc;
    double m_Kcgr;
    int    m_ce_index;
    int    m_cs_index;
};

} // namespace Thermodynamics
} // namespace Mutation

