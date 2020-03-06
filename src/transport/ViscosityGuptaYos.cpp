/**
 * @file ViscosityGuptaYos.cpp
 *
 * @brief Provides ViscosityGuptaYos class.
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

#include "AutoRegistration.h"
#include "CollisionDB.h"
#include "Thermodynamics.h"
#include "ViscosityAlgorithm.h"
#include "GuptaYos.h"

using namespace Mutation::Utilities;
using namespace Eigen;

namespace Mutation {
    namespace Transport {

/**
 * Uses the Gupta-Yos formula to compute viscosity.
 */
class ViscosityGuptaYos : public ViscosityAlgorithm, public GuptaYos
{
public:

    ViscosityGuptaYos(ViscosityAlgorithm::ARGS collisions)
        : ViscosityAlgorithm(collisions),
          A(collisions.nHeavy(), collisions.nHeavy()),
          B(collisions.nHeavy(), collisions.nHeavy()),
          a(collisions.nHeavy())
    { }

    /**
     * Returns the viscosity of the mixture in Pa-s.
     */
   double viscosity()
    {
        const int ns = m_collisions.nSpecies();
        const int nh = m_collisions.nHeavy();
        const int k  = ns - nh;
        
        const ArrayXd& nDij = m_collisions.nDij();
        const ArrayXd& Ast  = m_collisions.Astij();
        const ArrayXd& mi   = m_collisions.mass();
        const Map<const ArrayXd> x(m_collisions.thermo().X()+k, nh);
        
        // Compute the symmetric matrix A and vector a
        for (int j = 0, index = 0; j < nh; ++j) {
            for (int i = j; i < nh; ++i, ++index) {
                A(i,j) = (2.0-1.2*Ast(index))/((mi(i+k)+mi(j+k))*nDij(index));
                B(i,j) = Ast(index) / nDij(index);
            }
        }
                
        // Leave B as symmetric and postpone dividing by m(i) until after B*x
        // which uses fewer divides and allows for sym-matrix-vector product
        a.matrix() = B.selfadjointView<Lower>() * x.matrix();
        a *= 1.2 / mi.tail(nh);
            
        // Now compute the viscosity using Gupta-Yos
        return guptaYos(A, a, x);
    }
    
private:
    
    ArrayXXd A;
    MatrixXd B;
    ArrayXd  a;
};

Config::ObjectProvider<ViscosityGuptaYos, ViscosityAlgorithm> 
    viscosityGuptaYos("Gupta-Yos");
    

    } // namespace Transport
} // namespace Mutation


