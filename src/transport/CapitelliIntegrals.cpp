/**
 * @file CapitelliIntegrals.cpp
 *
 * Provides various curve-fits from the group of Capitelli.
 *
 * @author J.B. Scoggins
 */

/*
 * Copyright 2016-2017 von Karman Institute for Fluid Dynamics (VKI)
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
#include "CollisionPair.h"
#include "XMLite.h"

#include <eigen3/Eigen/Dense>
#include <string>

using namespace Mutation::Utilities;

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

/**
 * Represents collision integrals computed using the curve-fit found in Eq. (19)
 * of \cite{Bruno2010}.
 */
class BrunoEq19ColInt : public CapitelliIntegral
{
public:

    /**
     * Self registering constructor.
     */
    BrunoEq19ColInt(CollisionIntegral::ARGS args) :
        CapitelliIntegral(args)
    {
        // Load the coefficients
        std::istringstream ss(args.xml.text());
        for (int i = 0; i < 8; ++i)
            if(!(ss >> m_g[i]))
                args.xml.parseError("Must provide 8 coefficients.");
    }

    virtual bool canTabulate() const { return true; }

private:

    double compute_(double T)
    {
        double x  = std::log(T);
        double e1 = std::exp((x-m_g[0])/m_g[1]);
        double f1 = (x-m_g[6])/m_g[7];
        return m_g[2]*std::pow(x,m_g[4])*e1/(e1+1.0/e1)+m_g[5]*std::exp(-f1*f1)+m_g[3];
    }

    /**
     * Returns true if the coefficients are the same.
     */
    bool isEqual(const CollisionIntegral& ci) const {
        const BrunoEq19ColInt& compare =
            dynamic_cast<const BrunoEq19ColInt&>(ci);
        return (m_g == compare.m_g);
    }

private:

    Eigen::Matrix<double,8,1> m_g;

}; // class BrunoEq19ColInt

// Register the "Bruno-Eq(19)" CollisionIntegral
Mutation::Utilities::Config::ObjectProvider<
    BrunoEq19ColInt, CollisionIntegral> bruno19_ci("Bruno-Eq(19)");

//==============================================================================

class LaricchiutaEq15ColInt : public CapitelliIntegral
{
public:

    /**
     * Self registering constructor.
     */
    LaricchiutaEq15ColInt(CollisionIntegral::ARGS args) :
        CapitelliIntegral(args), m_phi0(0.0), m_sig2(0.0)
    {
        // Find Laricchiuta-Eq(15)-Data element
        IO::XmlElement::const_iterator pair = args.pair.findPair();
        if (pair == args.xml.document()->end())
            args.xml.parseError(
                "Must provide Laricchiuta-Eq(15)-Data to use Laricchiuta-Eq(15) integral.");

        IO::XmlElement::const_iterator data =
            pair->findTag("Laricchiuta-Eq(15)-Data");
        if (data == pair->end())
            args.xml.parseError(
                "Must provide Laricchiuta-Eq(15)-Data to use Laricchiuta-Eq(15) integral.");

        // Load the parameters associated with this pair
        double beta, re;
        data->getAttribute("beta", beta, "must have a 'beta' attribute.");
        data->getAttribute("re", re, "must have an 're' attribute.");
        data->getAttribute("phi0", m_phi0, "must have a 'phi0' attribute.");
        m_phi0 = Units("meV").convertToBase(m_phi0);
        std::string ref; data->getAttribute("ref", ref);
        setReference(ref);

        // Lastly, need to compute the fitting parameters based on this type
        Eigen::Matrix<double, 3,1> b; b << 1.0, beta, beta*beta;
        std::map<std::string, Eigen::Matrix<double,7,3> >& map = sm_c4;

        switch (args.pair.type()) {
        case ION_NEUTRAL:
            map = sm_c4;
            m_sig2 = 0.7564*std::pow(beta, 2.0*0.064605)*re*re;
            break;
        case NEUTRAL_NEUTRAL:
            map = sm_c6;
            m_sig2 = 0.8002*std::pow(beta, 2.0*0.049256)*re*re;
            break;
        default:
            args.xml.parseError(
                "Laricchiuta-Eq(15) only provided for ion-neutral and neutral-neutral collisions.");
        }

        std::map<std::string, Eigen::Matrix<double,7,3> >::const_iterator
            iter = map.find(args.kind);
        if (iter == map.end())
            args.xml.parseError(
                "Laricchiuta-Eq(15) is not supported for " + args.kind + " integral.");

        m_a = iter->second * b;
    }

private:

    double compute_(double T)
    {
        double x  = std::log(KB*T/m_phi0);
        double e1 = std::exp((x - m_a[2])/m_a[3]);
        double e2 = std::exp((x - m_a[5])/m_a[6]);

        return m_sig2*std::exp((m_a[0]+x*m_a[1])*e1/(e1+1.0/e1)+m_a[4]*e2/(e2+1.0/e2));
    }

    /**
     * Returns true if the coefficients are the same.
     */
    bool isEqual(const CollisionIntegral& ci) const {
        const LaricchiutaEq15ColInt& compare =
            dynamic_cast<const LaricchiutaEq15ColInt&>(ci);
        return (m_a == compare.m_a &&
                m_phi0 == compare.m_phi0 &&
                m_sig2 == compare.m_sig2);
    }

private:

    Eigen::Matrix<double,7,1> m_a;

    double m_phi0;
    double m_sig2;

    static std::map<std::string, Eigen::Matrix<double,7,3> > sm_c4;
    static std::map<std::string, Eigen::Matrix<double,7,3> > sm_c6;

}; // class LaricchiutaEq15ColInt

// Initialization of the ci coefficients for m = 4, Table 1 in Laricchiuta2007.
std::map<std::string, Eigen::Matrix<double,7,3> > init_c4()
{
    std::map<std::string, Eigen::Matrix<double,7,3> > c4_map;

    c4_map["Q11"] <<
         9.851755e-1, -2.870704e-2,  0.0,
        -4.737800e-1, -1.370344e-3,  0.0,
         7.080799e-1,  4.575312e-3,  0.0,
        -1.239439e+0, -3.683605e-2,  0.0,
        -4.638467e+0,  4.418904e-1, -1.220292e-2,
         3.841835e+0,  3.277658e-1, -2.660275e-2,
         2.317342e+0,  3.912768e-1, -3.136223e-2;

    c4_map["Q12"] <<
         8.361751e-1, -3.201292e-2,  0.0,
        -4.707355e-1, -1.783284e-3,  0.0,
         1.771157e-1,  1.172773e-2,  0.0,
        -1.094937e+0, -3.115598e-2,  0.0,
        -4.976384e+0,  4.708074e-1, -1.283818e-2,
         3.645873e+0,  3.699452e-1, -2.988684e-2,
         2.428864e+0,  4.267351e-1, -3.278874e-2;

    c4_map["Q13"] <<
         7.440562e-1, -3.453851e-2,  0.0,
        -4.656306e-1, -2.097901e-3,  0.0,
        -1.465752e-1,  1.446209e-2,  0.0,
        -1.080410e+0, -2.712029e-2,  0.0,
        -5.233907e+0,  4.846691e-1, -1.280346e-2,
         3.489814e+0,  4.140270e-1, -3.250138e-2,
         2.529678e+0,  4.515088e-1, -3.339293e-2;

    c4_map["Q14"] <<
         6.684360e-1, -3.515695e-2,  0.0,
        -4.622014e-1, -2.135808e-3,  0.0,
        -3.464990e-1,  1.336362e-2,  0.0,
        -1.054374e+0, -3.149321e-2,  0.0,
         5.465789e+0,  4.888443e-1, -1.228090e-2,
         3.374614e+0,  4.602468e-1, -3.463073e-2,
         2.648622e+0,  4.677409e-1, -3.339297e-2;

    c4_map["Q15"] <<
         6.299083e-1, -3.720000e-2,  0.0,
        -4.560921e-1, -2.395779e-3,  0.0,
        -5.228598e-1,  1.594610e-2,  0.0,
        -1.124725e+0, -2.862354e-2,  0.0,
        -5.687354e+0,  4.714789e-1, -1.056602e-2,
         3.267709e+0,  5.281419e-1, -3.678869e-2,
         2.784725e+0,  4.840700e-1, -3.265127e-2;

    c4_map["Q22"] <<
         9.124518e-1, -2.398461e-2,  0.0,
        -4.697184e-1, -7.809681e-4,  0.0,
         1.031053e+0,  4.069668e-3,  0.0,
        -1.090782e+0, -2.413508e-2,  0.0,
        -4.127243e+0,  4.302667e-1, -1.352874e-2,
         4.059078e+0,  2.597379e-1, -2.169951e-2,
         2.086906e+0,  2.920310e-1, -2.560437e-2;

    c4_map["Q23"] <<
         8.073459e-1, -2.581232e-2,  0.0,
        -4.663682e-1, -1.030271e-3,  0.0,
         6.256342e-1,  4.086881e-3,  0.0,
        -1.063437e+0, -1.235489e-2,  0.0,
        -4.365989e+0,  4.391454e-1, -1.314615e-2,
         3.854346e+0,  3.219224e-1, -2.587493e-2,
         2.146207e+0,  3.325620e-1, -2.686959e-2;

    c4_map["Q24"] <<
         7.324117e-1, -2.727580e-2,  0.0,
        -4.625614e-1, -1.224292e-3,  0.0,
         3.315871e-1,  7.216776e-3,  0.0,
        -1.055706e+0, -8.585500e-3,  0.0,
        -4.571022e+0,  4.373660e-1, -1.221457e-2,
         3.686006e+0,  3.854493e-1, -2.937568e-2,
         2.217893e+0,  3.641196e-1, -2.763824e-2;

    return c4_map;
}
std::map<std::string, Eigen::Matrix<double,7,3> > LaricchiutaEq15ColInt::sm_c4 = init_c4();

// Initialization of the ci coefficients for m = 6, Table 2 in Laricchiuta2007.
std::map<std::string, Eigen::Matrix<double,7,3> > init_c6()
{
    std::map<std::string, Eigen::Matrix<double,7,3> > c6_map;
    Eigen::Matrix<double,7,3> ci;

    c6_map["Q11"] <<
         7.884756e-1, -2.438494e-2,  0.0,
        -2.952759e-1, -1.744149e-3,  0.0,
         5.020892e-1,  4.316985e-2,  0.0,
        -9.042460e-1, -4.017103e-2,  0.0,
        -3.373058e+0,  2.458538e-1, -4.850047e-3,
         4.161981e+0,  2.202737e-1, -1.718010e-2,
         2.462523e+0,  3.231308e-1, -2.281072e-2;

    c6_map["Q22"] <<
         7.898524e-1, -2.114115e-2,  0.0,
        -2.998325e-1, -1.243977e-3,  0.0,
         7.077103e-1,  3.583907e-2,  0.0,
        -8.946857e-1, -2.473947e-2,  0.0,
        -2.958969e+0,  2.303358e-1, -5.226562e-3,
         4.348412e+0,  1.920321e-1, -1.496557e-2,
         2.205440e+0,  2.567027e-1, -1.861359e-2;

    return c6_map;
}
std::map<std::string, Eigen::Matrix<double,7,3> > LaricchiutaEq15ColInt::sm_c6 = init_c6();

// Register the "Laricchiuta-Eq(15)" CollisionIntegral
Mutation::Utilities::Config::ObjectProvider<
    LaricchiutaEq15ColInt, CollisionIntegral> laricchiuta15_ci("Laricchiuta-Eq(15)");

//==============================================================================

    } // namespace Transport
} // namespace Mutation
