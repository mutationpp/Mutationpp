
#include "Transport.h"
#include "Constants.h"
#include "Numerics.h"

#include <iostream>

using namespace Mutation::Numerics;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities;

namespace Mutation {
    namespace Transport {

using Mutation::Thermodynamics::Thermodynamics;

//==============================================================================

Transport::Transport(
    Thermodynamics& thermo, const std::string& viscosity, const std::string& lambda, const bool load_data)
    : m_thermo(thermo),
      mp_collisions(NULL),
      mp_viscosity(NULL),
      mp_thermal_conductivity(NULL),
      mp_diffusion_matrix(NULL),
      mp_wrk1(NULL),
      mp_tag(NULL)
{ 
    if (!load_data)
        return;

    // Load the collision integral data
    mp_collisions = new CollisionDB(thermo);

    // Load the viscosity calculator
    mp_viscosity = 
        Config::Factory<ViscosityAlgorithm>::create(viscosity, *mp_collisions);
    
    // Load the thermal conductivity calculator
    mp_thermal_conductivity = 
        Config::Factory<ThermalConductivityAlgorithm>::create(
            lambda, *mp_collisions);
    
    // Load the diffusion matrix calculator
    mp_diffusion_matrix =
        new Ramshaw(thermo, *mp_collisions);
    
    // Allocate work array storage
    mp_wrk1 = new double [m_thermo.nSpecies()*3];
    mp_wrk2 = mp_wrk1 + m_thermo.nSpecies();
    mp_wrk3 = mp_wrk2 + m_thermo.nSpecies();
    if(m_thermo.nEnergyEqns() > 1) {
        mp_tag  = new int [m_thermo.nEnergyEqns()*5];
        m_thermo.getTagModes(mp_tag);
    }
    
    //thermo.stateModel()->notifyOnUpdate(this);
}
    
//==============================================================================
    
Transport::~Transport()
{
    delete mp_collisions;
    delete mp_viscosity;
    delete mp_thermal_conductivity;
    delete mp_diffusion_matrix;
    
    delete [] mp_wrk1;
    delete [] mp_tag;
}

//==============================================================================

void Transport::omega11ii(double* const p_omega)
{
    if (mp_collisions == NULL) {
        cout << "Error! Trying to use transport without loading collision integrals!!" << endl;
        return;
    }

    const double ns  = m_thermo.nSpecies();
    const double Th  = m_thermo.T();
    const double Te  = m_thermo.Te();
    const double nd  = m_thermo.numberDensity();
    const double *const X = m_thermo.X();
    const RealSymMat& Q11 = mp_collisions->Q11(Th, Te, nd, X);
    
    for (int i = 0; i < ns; ++i)
        p_omega[i] = Q11(i,i);
}

//==============================================================================

void Transport::omega22ii(double* const p_omega)
{
    if (mp_collisions == NULL) {
        cout << "Error! Trying to use transport without loading collision integrals!!" << endl;
        return;
    }

    const double ns  = m_thermo.nSpecies();
    const double Th  = m_thermo.T();
    const double Te  = m_thermo.Te();
    const double nd  = m_thermo.numberDensity();
    const double *const X = m_thermo.X();
    const RealSymMat& Q22 = mp_collisions->Q22(Th, Te, nd, X);
    
    for (int i = 0; i < ns; ++i)
        p_omega[i] = Q22(i,i);
}

//==============================================================================

void Transport::frozenThermalConductivityVector(double* const p_lambda)
{

    const int neq = m_thermo.nEnergyEqns();
    double lambda_th, lambda_te, lambda_rot, lambda_vib, lambda_elec;

    if(neq > 1) {
        lambda_th   = heavyThermalConductivity(); 
        lambda_te   = electronThermalConductivity(); 
        lambda_rot  = rotationalThermalConductivity();
        lambda_vib  = vibrationalThermalConductivity();
        lambda_elec = electronicThermalConductivity();

        for (int i_eq = 0; i_eq < neq; ++i_eq) {
           p_lambda[i_eq] = mp_tag[i_eq*5+0] * lambda_th 
                          + mp_tag[i_eq*5+1] * lambda_te 
                          + mp_tag[i_eq*5+2] * lambda_rot
                          + mp_tag[i_eq*5+3] * lambda_vib
                          + mp_tag[i_eq*5+4] * lambda_elec;
        }
    } else {
        p_lambda[0] = frozenThermalConductivity();
    }

}

//==============================================================================

double Transport::electronThermalConductivity()
{
    if (mp_collisions == NULL) {
        cout << "Error! Trying to use transport without loading collision integrals!!" << endl;
        return 0.0;
    }

    if (!m_thermo.hasElectrons() || m_thermo.X()[0] == 0.0)
        return 0.0;

    // Get thermodynamic properties
    const int ns     = m_thermo.nSpecies();
    const double Th  = m_thermo.T();
    const double Te  = m_thermo.Te();
    const double nd  = m_thermo.numberDensity();
    const double me  = m_thermo.speciesMw(0)/NA;
    const double *const X = m_thermo.X();
    
    // Get collision integral information
    const RealSymMat& Q11   = mp_collisions->Q11(Th, Te, nd, X);
    const RealSymMat& Q22   = mp_collisions->Q22(Th, Te, nd, X);
    const RealSymMat& B     = mp_collisions->Bstar(Th, Te, nd, X);
    const RealVector& Q12ei = mp_collisions->Q12ei(Th, Te, nd, X);
    const RealVector& Q13ei = mp_collisions->Q13ei(Th, Te, nd, X);
    const RealVector& Q14ei = mp_collisions->Q14ei(Th, Te, nd, X);
    const RealVector& Q15ei = mp_collisions->Q15ei(Th, Te, nd, X);
    const double      Q23ee = mp_collisions->Q23ee(Th, Te, nd, X);
    const double      Q24ee = mp_collisions->Q24ee(Th, Te, nd, X);
    
    // Compute the lambdas
    double fac;
    double lam11 = 0.0;
    double lam12 = 0.0;
    double lam22 = 0.0;
    
    for (int i = 1; i < ns; ++i) {
        lam11 += X[i]*Q11(i)*(25.0/4.0-3.0*B(i));
        lam12 += X[i]*(175.0/16.0*Q11(i)-315.0/8.0*Q12ei(i)+57.0*Q13ei(i)-30.0*
            Q14ei(i));
        lam22 += X[i]*(1225.0/64.0*Q11(i)-735.0/8.0*Q12ei(i)+399.0/2.0*Q13ei(i)-
            210.0*Q14ei(i)+90.0*Q15ei(i));
    }
    
    fac = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0];
    lam11 = fac*(lam11 + SQRT2*X[0]*Q22(0));
    lam12 = fac*(lam12 + SQRT2*X[0]*(7.0/4.0*Q22(0)-2.0*Q23ee));
    lam22 = fac*(lam22 + SQRT2*X[0]*(77.0/16.0*Q22(0)-7.0*Q23ee+5.0*Q24ee));
    
    if (lam11 <= 0.0) {
        for (int i = 1; i < ns; ++i)
            cout << m_thermo.speciesName(i) << " " << 6.25 - 3.0*B(i) << " " << Q11(i)*X[i] << endl;
    }
    assert(lam11 > 0.0);
    //assert(lam22 > 0.0);
    //assert(lam11*lam22 > lam12*lam12);

    // 2nd order solution
    return (X[0]*X[0]/lam11);
    
    // 3rd order solution
    //return (X[0]*X[0]*lam22/(lam11*lam22-lam12*lam12));
}

//==============================================================================

double Transport::internalThermalConductivity()
{
    if (mp_collisions == NULL) {
        cout << "Error! Trying to use transport without loading collision integrals!!" << endl;
        return 0.0;
    }

    const int ns     = m_thermo.nSpecies();
    const int nh     = m_thermo.nHeavy();
    const double Th  = m_thermo.T();
    const double Te  = m_thermo.Te();
    const double Tr  = m_thermo.Tr();
    const double Tv  = m_thermo.Tv();
    const double Tel = m_thermo.Tel();
    const double nd  = m_thermo.numberDensity(); 
    const double *const X = m_thermo.X();
    
    const RealSymMat& nDij = mp_collisions->nDij(Th, Te, nd, X);
    
    m_thermo.speciesCpOverR(
        Th, Te, Tr, Tv, Tel, NULL, NULL, mp_wrk1, mp_wrk2, mp_wrk3);
    
    double lambda = 0.0;
    for (int i = ns-nh; i < ns; ++i) {
        double sum = 0.0;
        for (int j = 0; j < ns; ++j)
            sum += X[j] / nDij(i,j);
        
        lambda += X[i] * (mp_wrk1[i] + mp_wrk2[i] + mp_wrk3[i]) / sum;
    }
    
    return (lambda * KB);
}

//==============================================================================

double Transport::rotationalThermalConductivity()
{
    if (mp_collisions == NULL) {
        cout << "Error! Trying to use transport without loading collision integrals!!" << endl;
        return 0.0;
    }

    const int ns     = m_thermo.nSpecies();
    const int nh     = m_thermo.nHeavy();
    const double Th  = m_thermo.T();
    const double Te  = m_thermo.Te();
    const double Tr  = m_thermo.Tr();
    const double Tv  = m_thermo.Tv();
    const double Tel = m_thermo.Tel();
    const double nd  = m_thermo.numberDensity(); 
    const double *const X = m_thermo.X();
    
    const RealSymMat& nDij = mp_collisions->nDij(Th, Te, nd, X);
    
    m_thermo.speciesCpOverR(
        Th, Te, Tr, Tv, Tel, NULL, NULL, mp_wrk1, NULL, NULL);
    
    double lambda = 0.0;
    for (int i = ns-nh; i < ns; ++i) {
        double sum = 0.0;
        for (int j = 0; j < ns; ++j)
            sum += X[j] / nDij(i,j);
        
        lambda += X[i] * mp_wrk1[i] / sum;
    }
    
    return (lambda * KB);
}

//==============================================================================

double Transport::vibrationalThermalConductivity()
{
    if (mp_collisions == NULL) {
        cout << "Error! Trying to use transport without loading collision integrals!!" << endl;
        return 0.0;
    }

    const int ns     = m_thermo.nSpecies();
    const int nh     = m_thermo.nHeavy();
    const double Th  = m_thermo.T();
    const double Te  = m_thermo.Te();
    const double Tr  = m_thermo.Tr();
    const double Tv  = m_thermo.Tv();
    const double Tel = m_thermo.Tel();
    const double nd  = m_thermo.numberDensity(); 
    const double *const X = m_thermo.X();
    
    const RealSymMat& nDij = mp_collisions->nDij(Th, Te, nd, X);
    
    m_thermo.speciesCpOverR(
        Th, Te, Tr, Tv, Tel, NULL, NULL, NULL, mp_wrk1, NULL);
    
    double lambda = 0.0;
    for (int i = ns-nh; i < ns; ++i) {
        double sum = 0.0;
        for (int j = 0; j < ns; ++j)
            sum += X[j] / nDij(i,j);
        
        lambda += X[i] * mp_wrk1[i] / sum;
    }
    
    return (lambda * KB);
}

//==============================================================================

double Transport::electronicThermalConductivity()
{
    if (mp_collisions == NULL) {
        cout << "Error! Trying to use transport without loading collision integrals!!" << endl;
        return 0.0;
    }

    const int ns     = m_thermo.nSpecies();
    const int nh     = m_thermo.nHeavy();
    const double Th  = m_thermo.T();
    const double Te  = m_thermo.Te();
    const double Tr  = m_thermo.Tr();
    const double Tv  = m_thermo.Tv();
    const double Tel = m_thermo.Tel();
    const double nd  = m_thermo.numberDensity(); 
    const double *const X = m_thermo.X();
    
    const RealSymMat& nDij = mp_collisions->nDij(Th, Te, nd, X);
    
    m_thermo.speciesCpOverR(
        Th, Te, Tr, Tv, Tel, NULL, NULL, NULL, NULL, mp_wrk1);
    
    double lambda = 0.0;
    for (int i = ns-nh; i < ns; ++i) {
        double sum = 0.0;
        for (int j = 0; j < ns; ++j)
            sum += X[j] / nDij(i,j);
        
        lambda += X[i] * mp_wrk1[i] / sum;
    }
    
    return (lambda * KB);
}

//==============================================================================
/*double Transport::reactiveThermalConductivity()
{
    ///
    /// @todo This function fails if there is no electric field.  Does not take
    /// into account what happens when there are no electrons.
    ///

    // Compute dX_i/dT
    m_thermo.dXidT(m_thermo.T(), m_thermo.P(), m_thermo.X(), mp_wrk1);
    
    // Compute the multicomponent diffusion coefficient matrix
    const int ns = m_thermo.nSpecies();
    const RealMatrix& Dij = diffusionMatrix();
    const double* const X = m_thermo.X();
    const double* const Y = m_thermo.Y();
    
    // Store the species charges and determine the mixture charge
    double q = 0.0;
    for (int i = 0; i < ns; ++i) {
        mp_wrk3[i] = m_thermo.speciesCharge(i);
        q += mp_wrk3[i] * X[i];
    }
    
    for (int i = 0; i < ns; ++i) {
        mp_wrk2[i] = 0.0;
        for (int j = 0; j < ns; ++j) {
            mp_wrk2[i] += Dij(i,j)*(mp_wrk1[j]-(X[j]*mp_wrk3[j]-Y[j]*q)/
                (X[0]*mp_wrk3[0]-Y[0]*q)*mp_wrk1[0]);
        }
    }
    
    // Compute the species enthalpies per unit mass
    m_thermo.speciesHOverRT(mp_wrk1);
    
    double lambda = 0.0;
    for (int i = 0; i < ns; ++i)
        lambda += mp_wrk1[i] * mp_wrk2[i] * X[i];
    
    return (m_thermo.P()*lambda);
}*/

double Transport::reactiveThermalConductivity()
{
    if (mp_collisions == NULL) {
        cout << "Error! Trying to use transport without loading collision integrals!!" << endl;
        return 0.0;
    }

    // Compute dX_i/dT
    m_thermo.dXidT(mp_wrk1);
    
    // Compute the thermal diffusion ratios
    thermalDiffusionRatios(mp_wrk2);
    
    // Combine to get the driving forces
    for (int i = 0; i < m_thermo.nSpecies(); i++)
        mp_wrk1[i] += mp_wrk2[i] / m_thermo.T();
    
    // Compute the diffusion velocities
    double E;
    stefanMaxwell(mp_wrk1, mp_wrk2, E);
    
    // Unitless enthalpies
    m_thermo.speciesHOverRT(mp_wrk1);
    
    double lambda = 0.0;
    //const double* const X = m_thermo.X();
    const double rho = m_thermo.density();
    const double* const Y = m_thermo.Y();
    for (int i = 0; i < m_thermo.nSpecies(); ++i)
        lambda -= mp_wrk1[i] / m_thermo.speciesMw(i) * mp_wrk2[i] * Y[i] * rho;
    
    return (RU * m_thermo.T() * lambda);
}

double Transport::soretThermalConductivity()
{
    if (mp_collisions == NULL) {
        cout << "Error! Trying to use transport without loading collision integrals!!" << endl;
        return 0.0;
    }

    // Compute dX_i/dT
    m_thermo.dXidT(mp_wrk1);
    
    // Compute the thermal diffusion ratios
    thermalDiffusionRatios(mp_wrk2);
    
    // Combine to get the driving forces
    for (int i = 0; i < m_thermo.nSpecies(); i++)
        mp_wrk1[i] += mp_wrk2[i] / m_thermo.T();
    
    // Compute the diffusion velocities
    double E;
    stefanMaxwell(mp_wrk1, mp_wrk1, E);
    
    double lambda = 0.0;
    for (int i = 0; i < m_thermo.nSpecies(); i++)
        lambda -= mp_wrk2[i]*mp_wrk1[i];
    
    return (m_thermo.P()*lambda);
}

//==============================================================================

void Transport::averageDiffusionCoeffs(double *const p_Di)
{
    if (mp_collisions == NULL) {
        cout << "Error! Trying to use transport without loading collision integrals!!" << endl;
        return;
    }

    const int ns = m_thermo.nSpecies();
    const double Th = m_thermo.T();
    const double Te = m_thermo.Te();
    const double nd = m_thermo.numberDensity();
    const double* const p_X = m_thermo.X();
    const RealSymMat& nDij = mp_collisions->nDij(Th, Te, nd, p_X);
    
    for (int i = 0; i < ns; ++i)
        p_Di[i] = 0.0;
    
    for (int i = 0, index = 1; i < ns; ++i, ++index) {
        for (int j = i + 1; j < ns; ++j, ++index) {
            p_Di[i] += p_X[j] / nDij(index);
            p_Di[j] += p_X[i] / nDij(index);
        }
    }
    
    for (int i = 0; i < ns; ++i)
        p_Di[i] = (1.0 - p_X[i]) / (p_Di[i] * nd);
}

void Transport::equilibriumFickP(double* const p_F)
{
    if (mp_collisions == NULL) {
        cout << "Error! Trying to use transport without loading collision integrals!!" << endl;
        return;
    }

    // Get some state data
    const int ns = m_thermo.nSpecies();
    const int ne = m_thermo.nElements();
    const double* const p_Y = m_thermo.Y();
    const double* const p_X = m_thermo.X();
    const double rho = m_thermo.density();
    const double p   = m_thermo.P();

    const RealMatrix& Dij = diffusionMatrix();
    const RealMatrix& nu  = m_thermo.elementMatrix();

    // Compute the dXj/dP term
    m_thermo.dXidP(mp_wrk1);

    for (int i = 0; i < ns; ++i) {
        mp_wrk2[i] = 0.0;
        for (int j = 0; j < ns; ++j)
            mp_wrk2[i] += Dij(i,j)*(p_X[j]/p + mp_wrk1[j]);
        mp_wrk2[i] *= -rho*p_Y[i];
    }

    for (int k = 0; k < ne; ++k) {
        p_F[k] = 0.0;
        for (int i = 0; i < ns; ++i)
            p_F[k] += nu(i,k)*mp_wrk2[i];
    }
}

void Transport::equilibriumFickT(double* const p_F)
{
    if (mp_collisions == NULL) {
        cout << "Error! Trying to use transport without loading collision integrals!!" << endl;
        return;
    }

    // Get some state data
    const int ns = m_thermo.nSpecies();
    const int ne = m_thermo.nElements();
    const double* const p_Y = m_thermo.Y();
    const double* const p_X = m_thermo.X();
    const double rho = m_thermo.density();

    const RealMatrix& Dij = diffusionMatrix();
    const RealMatrix& nu  = m_thermo.elementMatrix();

    // Compute the dXj/dP term
    m_thermo.dXidT(mp_wrk1);

    for (int i = 0; i < ns; ++i) {
        mp_wrk2[i] = 0.0;
        for (int j = 0; j < ns; ++j)
            mp_wrk2[i] += Dij(i,j)*mp_wrk1[j];
        mp_wrk2[i] *= -rho*p_Y[i];
    }

    for (int k = 0; k < ne; ++k) {
        p_F[k] = 0.0;
        for (int i = 0; i < ns; ++i)
            p_F[k] += nu(i,k)*mp_wrk2[i];
    }
}

void Transport::equilibriumFickXe(double* const p_F)
{
    if (mp_collisions == NULL) {
        cout << "Error! Trying to use transport without loading collision integrals!!" << endl;
        return;
    }

    // Get some state data
    const int ns = m_thermo.nSpecies();
    const int ne = m_thermo.nElements();
    const double* const p_Y = m_thermo.Y();
    const double* const p_X = m_thermo.X();
    const double rho = m_thermo.density();
    const double T   = m_thermo.T();
    const double p   = m_thermo.P();

    const RealMatrix& Dij = diffusionMatrix();
    const RealMatrix& nu  = m_thermo.elementMatrix();

    for (int l = 0; l < ne; ++l) {
       // Compute the dXj/dZl term using a finite difference
       m_thermo.elementFractions(p_X, mp_wrk1);
       double h = std::max(mp_wrk1[l]*1.0e-6, 1.0e-6);
       mp_wrk1[l] += h;
       m_thermo.equilibriumComposition(T, p, mp_wrk1, mp_wrk2);

       for (int i = 0; i < ns; ++i)
           mp_wrk1[i] = (mp_wrk2[i]-p_X[i])/h;

       for (int i = 0; i < ns; ++i) {
           mp_wrk2[i] = 0.0;
           for (int j = 0; j < ns; ++j)
               mp_wrk2[i] += Dij(i,j)*mp_wrk1[j];
           mp_wrk2[i] *= -rho*p_Y[i];
       }

       for (int k = 0; k < ne; ++k) {
           double& Fkl = p_F[l*ne+k];
           Fkl = 0.0;
           for (int i = 0; i < ns; ++i)
               Fkl += nu(i,k)*mp_wrk2[i];
       }
    }

   // Be sure to set the state back in the equilibrium solver in case other
   // calculations rely on the correct element potential values
   m_thermo.elementFractions(p_X, mp_wrk1);
   m_thermo.equilibriumComposition(T, p, mp_wrk1, mp_wrk2);
}

//==============================================================================

void Transport::stefanMaxwell(
    const double* const p_dp, double* const p_V, double& E)
{
    if (mp_collisions == NULL) {
        cout << "Error! Trying to use transport without loading collision integrals!!" << endl;
        return;
    }

    const int ns = m_thermo.nSpecies();
    const double Th = m_thermo.T();
    const double Te = m_thermo.Te();
    const double nd = m_thermo.numberDensity();
    
    // Need to place a tolerance on X and Y
    const double tol = 1.0e-16;
    static std::vector<double> X(ns);
    static std::vector<double> Y(ns);
    for (int i = 0; i < ns; ++i)
        X[i] = std::max(tol, m_thermo.X()[i]);
    m_thermo.convert<X_TO_Y>(&X[0], &Y[0]);

    // Get reference to binary diffusion coefficients
    const RealSymMat& nDij = mp_collisions->nDij(Th, Te, nd, &X[0]);

    // Electron offset
    const int k = m_thermo.hasElectrons() ? 1 : 0;

    // Compute mixture charge
    double q = 0.0;
    for (int i = 0; i < ns; ++i) {
        mp_wrk3[i] = m_thermo.speciesCharge(i);
        q += mp_wrk3[i] * X[i];
    }

    // Compute kappas
    const double one_over_kbt = 1.0/(KB*Th);
    for (int i = 0; i < ns; ++i)
        p_V[i] = (X[i]*mp_wrk3[i] - Y[i]*q)*one_over_kbt;

    // Compute ambipolar electric field
    E = (k > 0 ? p_dp[0]/p_V[0] : 0.0);

    // Compute right hand side
    static RealVector b(ns-k);
    for (int i = k; i < ns; ++i)
        b(i-k) = -p_dp[i] + p_V[i]*E;

    // Compute singular system matrix
    double fac;
    static RealSymMat G(ns-k);
    for (int i = k; i < ns; ++i)
        G(i-k,i-k) = 0.0;
    for (int i = k; i < ns; ++i) {
        for (int j = i+1; j < ns; ++j) {
            fac = X[i]*X[j]/nDij(i,j)*nd;
            G(i-k,i-k) += fac;
            G(j-k,j-k) += fac;
            G(i-k,j-k) = -fac;
        }
    }

    // Compute average binary diffusion coefficient
    double a = 0.0;
    for (int i = 0; i < nDij.size(); ++i)
        a += nDij(i);
    a /= ((double)nDij.size()*nd);
    a = 1.0/a;

    // Compute the constraints that are applied to the matrix
    fac = (k > 0 ? Y[0]/(X[0]*mp_wrk3[0]) : 0.0);
    for (int i = k; i < ns; ++i)
        p_V[i] = Y[i] - X[i]*mp_wrk3[i]*fac;

    // Add mass balance relation to make make matrix nonsigular
    for (int i = k; i < ns; ++i)
        for (int j = i; j < ns; ++j)
            G(i-k,j-k) += a*p_V[i]*p_V[j];

    // Solve for the velocities
    static RealVector V(ns-k);
    static LDLT<double> ldlt;
    ldlt.setMatrix(G);
    ldlt.solve(V, b);

    // Copy to output vector
    for (int i = k; i < ns; ++i)
        p_V[i] = V(i-k);

    // Compute electron diffusion velocity if need be
    if (k == 0)
        return;

    p_V[0] = 0.0;
    for (int i = k; i < ns; ++i)
        p_V[0] += X[i]*mp_wrk3[i]*p_V[i];
    p_V[0] /= -mp_wrk3[0]*X[0];
}

// For now assume that we have ions and solve the full system in thermal
// nonequilibrium with ambipolar assumption
/*void Transport::stefanMaxwell(
    const double* const p_dp, double* const p_V, double& E)
{
    // Determine some constants
    const int ns = m_thermo.nSpecies();
    const double Th = m_thermo.T();
    const double Te = m_thermo.Te();
    const double nd = m_thermo.numberDensity();

    // Place a tolerance on X and Y
    const double tol = 1.0e-16;
    static std::vector<double> xy(2*ns);
    double *X = &xy[0], *Y = &xy[ns];
    double sum = 0.0;
    for (int i = 0; i < ns; ++i) {
        X[i] = tol + m_thermo.X()[i];
        sum += X[i];
    }
    for (int i = 0; i < ns; ++i)
        X[i] /= sum;
    m_thermo.convert<X_TO_Y>(X, Y);

    // Get reference to binary diffusion coefficients
    const RealSymMat& nDij = mp_collisions->nDij(Th, Te, nd, &X[0]);

    // Compute mixture charge (store species' charges in mp_wrk3 work array)
    double q = 0.0;
    for (int i = 0; i < ns; ++i) {
        mp_wrk3[i] = m_thermo.speciesCharge(i);
        q += mp_wrk3[i] * X[i];
    }

    // Compute kappas (stored in p_V to avoid creating another array)
    // also compute the 2-norm of the kappa vector
    const double one_over_kbt = 1.0/(KB*Th);
    double s = 0.0;
    for (int i = 0; i < ns; ++i) {
        p_V[i] = (X[i]*mp_wrk3[i] - Y[i]*q)*one_over_kbt;
        s += p_V[i]*p_V[i];
    }
    s = std::sqrt(s);

    // Compute the RHS vector
    static RealVector b(ns+1);
    for (int i = 0; i < ns; ++i)
        b(i) = -p_dp[i];
    b(0) *= Th/Te;
    b(ns) = 0.0;
    //cout << b << endl;
    // Compute system matrix
    static RealMatrix G(ns+1,ns+1);

    // - electron contribution
    double fac1 = Te/Th;
    double fac2 = fac1*X[0]*nd;
    G(0,0) = 0.0;
    for (int j = 1; j < ns; ++j) {
        G(0,j) =  -fac2*X[j]/nDij(0,j);
        G(0,0) += G(0,j);
        G(j,0) =  G(0,j);
        G(j,j) = -fac1*G(j,0);
    }
    G(0,0) /= -fac1;
    G(0,ns) = -p_V[0]/(fac1*s);

    // - heavy contribution
    for (int i = 1; i < ns; ++i) {
        for (int j = i+1; j < ns; ++j) {
            G(i,j) =  -X[i]*X[j]/nDij(i,j)*nd;
            G(i,i) -= G(i,j);
            G(j,i) =  G(i,j);
            G(j,j) -= G(j,i);
        }
        G(i,ns) = -p_V[i]/s;
    }

    // - ambipolar constraint
    for (int j = 0; j < ns; ++j)
        G(ns,j) = -p_V[j]/s;
    G(ns,ns) = 0.0;

    //cout << G << endl;

    // Finally solve the system for Vi and E
    static RealVector x(ns+1);
    std::pair<int, double> ret = Numerics::gmres(
        G, x, b, Numerics::DiagonalPreconditioner<double>(G));

    // Compute mass constraint projector
    static RealVector R(ns,1.0);
    R(0) /= fac1;
    double r = 0.0;
    for (int i = 0; i < ns; ++i)
        r += R(i)*Y[i];
    double p = 0.0;
    for (int i = 0; i < ns; ++i)
        p += x(i)*Y[i];
    p /= r;

    // Retrieve the solution
    for (int i = 0; i < ns; ++i)
        p_V[i] = x(i) - p*R(i);
    E = b(ns)/s;
}*/

//==============================================================================

double Transport::sigma() 
{
    if (mp_collisions == NULL) {
        cout << "Error! Trying to use transport without loading collision integrals!!" << endl;
        return 0.0;
    }

    if (!m_thermo.hasElectrons() || m_thermo.X()[0] < 1.0e-30) 
        return 0.0;

    const int ns = m_thermo.nSpecies();
    const double Th = m_thermo.T();
    const double Te = m_thermo.Te();
    const double nd = m_thermo.numberDensity();
    const double* const X = m_thermo.X();
    const double me = m_thermo.speciesMw(0)/NA;

    const RealSymMat& Q11 = mp_collisions->Q11(Th, Te, nd, X);
    const RealSymMat& Q22 = mp_collisions->Q22(Th, Te, nd, X);
    const RealSymMat& B   = mp_collisions->Bstar(Th, Te, nd, X);
    const RealVector& Cei = mp_collisions->Cstei(Th, Te, nd, X);
    
    // Compute lambdas
    double lam00 = 0.0;
    double lam01 = 0.0;
    double lam11 = 0.0;
    double fac   = 0.0;
    
    for (int i = 1; i < ns; ++i) {
        fac = X[i]*Q11(i);
        lam00 += fac;
        lam01 += fac*(2.5-3.0*Cei(i));
        lam11 += fac*(25.0/4.0-3.0*B(i));
    }
    
    fac = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0];
    lam00 = fac*lam00;
    lam01 = fac*lam01;
    lam11 = fac*(lam11 + SQRT2*X[0]*Q22(0));
    
    // First order
    //return (4.0/25.0*(X[0]*X[0]*QE*QE)/(KB*KB*Te*lam00));
    
    // Second order
    return (4.0/25.0*(X[0]*X[0]*QE*QE)/(KB*KB*Te*(lam00-lam01*lam01/lam11)));
}

//==============================================================================

    } // namespace Transport
} // namespace Mutation

