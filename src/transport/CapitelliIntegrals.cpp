/**
 * @file CapitelliIntegrals.cpp
 *
 * Provides various curve-fits from the group of Capitelli.
 *
 * @author J.B. Scoggins
 */

/*
 * Copyright 2016 von Karman Institute for Fluid Dynamics (VKI)
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
#include "Constants.h"
#include "CollisionIntegral.h"
#include "XMLite.h"

#include <Eigen/Dense>
#include <string>

namespace Mutation {
    namespace Transport {

//==============================================================================

class CapitelliIntegral : public CollisionIntegral
{
public:

    /**
     * Self registering constructor.
     */
    CapitelliIntegral(CollisionIntegral::ARGS args) :
        CollisionIntegral(args)
    {
        setFactor(PI);
        setUnits("A-A");
    }

    bool canTabulate() const { return true; }
};

//==============================================================================

/**
 * Represents collision integrals computed using the curve-fit found in Eq. (11)
 * of \cite{Bruno2010}.
 */
class BrunoEq11ColInt : public CapitelliIntegral
{
public:

    /**
     * Self registering constructor.
     */
    BrunoEq11ColInt(CollisionIntegral::ARGS args) :
        CapitelliIntegral(args)
    {
        // Load the coefficients
        std::istringstream ss(args.xml.text());
        for (int i = 0; i < 7; ++i)
            if(!(ss >> m_a[i]))
                args.xml.parseError("Must provide 7 coefficients.");
    }

private:

    double compute_(double T)
    {
        double x  = std::log(T);
        double e1 = std::exp((x - m_a[2])/m_a[3]);
        double e2 = std::exp((x - m_a[5])/m_a[6]);

        return (m_a[0]+x*m_a[1])*e1/(e1+1.0/e1) + m_a[4]*e2/(e2+1.0/e2);
    }

    /**
     * Returns true if the coefficients are the same.
     */
    bool isEqual(const CollisionIntegral& ci) const {
        const BrunoEq11ColInt& compare =
            dynamic_cast<const BrunoEq11ColInt&>(ci);
        return (m_a == compare.m_a);
    }

private:

    Eigen::Matrix<double,7,1> m_a;

}; // class BrunoEq11ColInt

// Register the "Bruno-Eq(11)" CollisionIntegral
Mutation::Utilities::Config::ObjectProvider<
    BrunoEq11ColInt, CollisionIntegral> bruno11_ci("Bruno-Eq(11)");

//==============================================================================

/**
 * Represents collision integrals computed using the curve-fit found in Eq. (17)
 * of \cite{Bruno2010}.
 */
class BrunoEq17ColInt : public CapitelliIntegral
{
public:

    /**
     * Self registering constructor.
     */
    BrunoEq17ColInt(CollisionIntegral::ARGS args) :
        CapitelliIntegral(args)
    {
        // Load the coefficients
        std::istringstream ss(args.xml.text());
        for (int i = 0; i < 3; ++i)
            if(!(ss >> m_d[i]))
                args.xml.parseError("Must provide 3 coefficients.");
    }

    virtual bool canTabulate() const { return true; }

private:

    double compute_(double T)
    {
        double x  = std::log(T);
        return m_d[0] + x*(m_d[1] + x*m_d[2]);
    }

    /**
     * Returns true if the coefficients are the same.
     */
    bool isEqual(const CollisionIntegral& ci) const {
        const BrunoEq17ColInt& compare =
            dynamic_cast<const BrunoEq17ColInt&>(ci);
        return (m_d == compare.m_d);
    }

private:

    Eigen::Vector3d m_d;

}; // class BrunoEq17ColInt

// Register the "Bruno-Eq(17)" CollisionIntegral
Mutation::Utilities::Config::ObjectProvider<
    BrunoEq17ColInt, CollisionIntegral> bruno17_ci("Bruno-Eq(17)");

//==============================================================================

    } // namespace Transport
} // namespace Mutation
