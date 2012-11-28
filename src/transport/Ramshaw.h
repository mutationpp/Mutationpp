#ifndef TRANSPORT_RANSHAW_H
#define TRANSPORT_RANSHAW_H

#include "CollisionDB.h"
#include "Thermodynamics.h"
#include "Numerics.h"

class Ramshaw
{
public:
    
    Ramshaw(const Thermodynamics& thermo, CollisionDB& collisions)
        : m_thermo(thermo), m_collisions(collisions), 
          m_Di(thermo.nSpecies()), m_D(thermo.nSpecies(), thermo.nSpecies())
    { }
    
    /**
     * Computes the binary diffusion coefficient matrix.
     *
     * @todo This method currently breaks the symmetry property of the diffusion
     * matrix.  It should be reformulated to return a symmetric matrix which
     * would be faster to compute and use.
     */
    const RealMatrix& diffusionMatrix(
        const double T, const double nd, const double *const p_x)
    {
        const int ns = m_thermo.nSpecies();
    
        // First step is to compute nDij and Y
        const Numerics::RealSymMat nDij = m_collisions.nDij(T, T, nd, p_x);
        double *Y = new double [ns];
        m_thermo.convert<X_TO_Y>(p_x, Y);
        
        // Now we can compute the mixture averaged diffusion coefficient for the
        // ith species
        m_Di = 0.0;
        for (int i = 0, index = 1; i < ns; ++i, ++index) {
            for (int j = i + 1; j < ns; ++j, ++index) {
                m_Di(i) += p_x[j] / nDij(index);
                m_Di(j) += p_x[i] / nDij(index);
            }
        }
        
        for (int i = 0; i < ns; ++i)
            m_Di(i) = (1.0 - Y[i]) / (m_Di(i) * p_x[i] * nd);
        
        // Form the diffusion matrix
        for (int i = 0; i < ns; ++i) {
            for (int j = 0; j < ns; ++j) {
                m_D(i,j) = -Y[j] * m_Di(j);
            }
            m_D(i,i) += m_Di(i);
        }
        
        delete [] Y;
        return m_D;
    }
    
private:

    const Thermodynamics& m_thermo;
    CollisionDB& m_collisions;
    Numerics::RealVector m_Di;
    Numerics::RealMatrix m_D;
        
}; // class Ramshaw

#endif // TRANSPORT_RANSHAW_H

