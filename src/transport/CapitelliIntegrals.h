/**
 * @file CapitelliIntegrals.h
 *
 * Declaration of collision integral types coming from the Capitelli group.
 *
 * @author J.B. Scoggins
 */

/*
 * Copyright 2017-2020 von Karman Institute for Fluid Dynamics (VKI)
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

#ifndef TRANSPORT_CAPITELLI_INTEGRALS_H
#define TRANSPORT_CAPITELLI_INTEGRALS_H

#include "CollisionIntegral.h"
#include "Constants.h"

namespace Mutation {
namespace Transport {

/**
 * Base class for all collision integrals from the Capitelli group.
 */
class CapitelliIntegral : public CollisionIntegral
{
public:
    CapitelliIntegral(CollisionIntegral::ARGS args) :
        CollisionIntegral(args)
    {
        setFactor(PI);
        setUnits("Å-Å");
    }
    bool canTabulate() const { return true; }
};


/**
 * Represents collision integrals computed using the curve-fit found in Eq. (11)
 * of \cite{Bruno2010}.
 */
class BrunoEq11ColInt : public CapitelliIntegral
{
public:
    BrunoEq11ColInt(CollisionIntegral::ARGS args);
private:
    double compute_(double T);

    /**
     * Returns true if the coefficients are the same.
     */
    bool isEqual(const CollisionIntegral& ci) const;

private:
    Eigen::Matrix<double,7,1> m_a;
}; // class BrunoEq11ColInt


/**
 * Represents collision integrals computed using the curve-fit found in Eq. (17)
 * of \cite{Bruno2010}.
 */
class BrunoEq17ColInt : public CapitelliIntegral
{
public:
    BrunoEq17ColInt(CollisionIntegral::ARGS args);
    virtual bool canTabulate() const { return true; }
private:
    double compute_(double T);

    /**
     * Returns true if the coefficients are the same.
     */
    bool isEqual(const CollisionIntegral& ci) const;

private:
    Eigen::Vector3d m_d;
}; // class BrunoEq17ColInt


/**
 * Represents collision integrals computed using the curve-fit found in Eq. (19)
 * of \cite{Bruno2010}.
 */
class BrunoEq19ColInt : public CapitelliIntegral
{
public:
    BrunoEq19ColInt(CollisionIntegral::ARGS args);
    virtual bool canTabulate() const { return true; }
private:
    double compute_(double T);

    /**
     * Returns true if the coefficients are the same.
     */
    bool isEqual(const CollisionIntegral& ci) const;

private:
    Eigen::Matrix<double,8,1> m_g;
}; // class BrunoEq19ColInt


/**
 * Computes the collision integral associated with the Pirani potential, using
 * the approximate curve-fits of Laricchiuta et al. @cite Laricchiuta2007.
 */
class PiraniColInt : public CapitelliIntegral
{
public:
    PiraniColInt(CollisionIntegral::ARGS args);
    virtual bool loaded() const { return m_loaded; }

private:
    double compute_(double T);

    /**
     * Returns true if the coefficients are the same.
     */
    bool isEqual(const CollisionIntegral& ci) const;

    /**
     * Searches the database and puts together the potential parameters, if
     * possible.
     */
    void loadPotentialParameters(
        CollisionIntegral::ARGS args, double& beta, double& re, double& phi0);
private:
    bool m_loaded;
    Eigen::Matrix<double,7,1> m_a;
    double m_phi0;
    double m_sig2;
    static std::map<std::string, Eigen::Matrix<double,7,3> > sm_c4;
    static std::map<std::string, Eigen::Matrix<double,7,3> > sm_c6;
}; // class PiraniColInt

} // namespace Transport
} // namespace Mutation

#endif // TRANSPORT_CAPITELLI_INTEGRALS_H
