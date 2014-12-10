/**
 * @file ViscosityGuptaYos.cpp
 *
 * @brief Provides ViscosityGuptaYos class.
 */

/*
 * Copyright 2014 von Karman Institute for Fluid Dynamics (VKI)
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

#include "ViscosityAlgorithm.h"
#include "GuptaYos.h"
#include "Numerics.h"

using namespace Mutation::Numerics;
using namespace Mutation::Utilities;

namespace Mutation {
    namespace Transport {

/**
 * Uses the Gupta-Yos formula to compute viscosity.
 */
class ViscosityGuptaYos : public ViscosityAlgorithm, public GuptaYos
{
public:

    ViscosityGuptaYos(CollisionDB& collisions) 
        : ViscosityAlgorithm(collisions),
          m_shift(collisions.nSpecies()-collisions.nHeavy()),
          a(collisions.nSpecies() - m_shift),
          B(collisions.nSpecies() - m_shift),
          A(collisions.nSpecies() - m_shift)
    { }

    /**
     * Returns the viscosity of the mixture in Pa-s.
     */
    double viscosity(const double T, const double nd, const double *const p_x)
    {
        const int ns = m_collisions.nSpecies();
        
        const RealSymMat& nDij  = m_collisions.nDij(T, T, nd, p_x);
        const RealSymMat& Astar = m_collisions.Astar(T, T, nd, p_x);
        const RealVector& mass  = m_collisions.mass();
        
        // Compute the symmetric matrix a and vector A
        for (int i = m_shift; i < ns; ++i) {
            for (int j = i; j < ns; ++j) {
                a(i-m_shift,j-m_shift) = (2.0 - 1.2 * Astar(i,j)) /
                    ((mass(i) + mass(j)) * nDij(i,j));
                B(i-m_shift,j-m_shift) = Astar(i,j) / nDij(i,j);
            }
        }
                
        // Leave B as symmetric and postpone dividing by m(i) until after B*x
        // which uses fewer divides and allows for sym-matrix-vector product
        A = B * asVector(p_x, ns)(m_shift,ns);
        A = 1.2 * A / mass(m_shift,ns);
            
        // Now compute the viscosity using Gupta-Yos
        return guptaYos(a, A, asVector(p_x, ns)(m_shift,ns));
    }
    
private:
    
    int m_shift;

    RealSymMat a;
    RealSymMat B;
    RealVector A;
};

Config::ObjectProvider<ViscosityGuptaYos, ViscosityAlgorithm> 
    viscosityGuptaYos("Gupta-Yos");
    

    } // namespace Transport
} // namespace Mutation


