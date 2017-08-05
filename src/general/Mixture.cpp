/**
 * @file Mixture.cpp
 *
 * @brief Mixture class implementation. @see Mutation::Mixture
 */

/*
 * Copyright 2014-2017 von Karman Institute for Fluid Dynamics (VKI)
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

#include "Mixture.h"
#include "StateModel.h"

#include <iostream>
using namespace std;
#include <eigen3/Eigen/Dense>
using namespace Eigen;

using namespace Mutation::Thermodynamics;

namespace Mutation {

//==============================================================================

Mixture::Mixture(const MixtureOptions& options)
    : Thermodynamics::Thermodynamics(
        options.getSpeciesDescriptor(),
        options.getThermodynamicDatabase(),
        options.getStateModel()),
      Transport(
        *this,
        options.getViscosityAlgorithm(),
        options.getThermalConductivityAlgorithm()),
      Kinetics(
        static_cast<const Thermodynamics&>(*this),
        options.getMechanism())
{
    // Add all the compositions given in mixture options to the composition list
    for (int i = 0; i < options.compositions().size(); ++i)
        addComposition(options.compositions()[i]);

    // Set default composition if available
    if (options.hasDefaultComposition())
        setDefaultComposition(m_compositions[options.getDefaultComposition()]);
    
    // Instantiate a new energy transfer model
    state()->initializeTransferModel(*this);
    
}

//==============================================================================

void Mixture::addComposition(const Composition& c, bool make_default)
{
    // Check if all names in the composition are elements
    bool elements = true;
    for (int i = 0; i < c.size(); ++i) {
        if (elementIndex(c[i].name) < 0) {
            elements = false;

            if (speciesIndex(c[i].name) < 0) {
                std::cerr << "Error: composition '" << c.name()
                      << "' has component which is not an element or species"
                      << " belonging to the mixture!" << std::endl;
                std::exit(1);
            }
        }
    }

    // If this composition has species names, then treat all as species and
    // convert to elements
    if (!elements) {
        ArrayXd svals(nSpecies());
        ArrayXd evals(nElements());
        c.getComposition(m_species_indices, svals.data());

        if (c.type() == Composition::MASS)
            convert<Y_TO_X>(svals.data(), svals.data());
        convert<X_TO_XE>(svals.data(), evals.data());

        // Get list of element names
        std::vector<std::string> names;
        for (int i = 0; i < nElements(); ++i)
            names.push_back(elementName(i));

        m_compositions.push_back(
            Composition(names, evals.data(), Composition::MOLE));
    } else
        m_compositions.push_back(c);

    if (make_default)
        setDefaultComposition(m_compositions.back());
}

//==============================================================================

bool Mixture::getComposition(
    const std::string& name, double* const p_vec, Composition::Type type) const
{
    int i = 0;
    for ( ; i < m_compositions.size(); ++i)
        if (m_compositions[i].name() == name) break;

    // Check if there is a composition with the given name
    if (i == m_compositions.size())
        return false;

    // Get the composition
    m_compositions[i].getComposition(m_element_indices, p_vec);

    // Convert if necessary
    if (m_compositions[i].type() != type) {
        if (type == Composition::MOLE)
            convert<YE_TO_XE>(p_vec, p_vec);
        else
            convert<XE_TO_YE>(p_vec, p_vec);
    }

    return true;
}

//==============================================================================

double Mixture::smb(
    double Tw, double Pw, double* const p_rhoi, const double* const p_xip, double dx)
{
    const double eps = 1.0e-8;
    const double tol = 1.0e-6;
    const double alpha = 0.05;
    const int maxits = 20;
    const int ns = nSpecies();

    Map<VectorXd> rhoi(p_rhoi, nSpecies());
    VectorXd rhoit(nSpecies());
    Map<const VectorXd> xip(p_xip, nSpecies());
    VectorXd Tvec(nEnergyEqns()); Tvec.fill(Tw);

    // Set the state
    setState(rhoi.data(), Tvec.data(), 1);

    // Compute the driving forces
    VectorXd dp = (xip - Map<const VectorXd>(X(), nSpecies())) / dx;

    // Compute f
    VectorXd f(ns); smbf(Tw, Pw, dp, f);
    VectorXd fp(ns);
    MatrixXd jac(ns, ns);
    VectorXd d(ns);
    cout << "norm(f) = " << f.norm() << endl;
    cout << "f = \n" << f << endl;

    for (int newt = 0; newt < maxits && f.norm() > tol; ++newt) {
        // Compute the Jacobian
        for (int j = 0; j < ns; ++j) {
            // Perturb mass fraction of species j
            double h = std::max(eps*rhoi[j], 1.0e-20);
            double rhoj = rhoi[j];
            rhoi[j] += h;

            // Now set the state
            setState(rhoi.data(), Tvec.data(), 1);

            // Compute f
            dp = (xip - Map<const VectorXd>(X(), ns)) / dx;
            smbf(Tw, Pw, dp, fp);

            // Column j of jacobian
            jac.col(j) = (fp - f) / h;
            jac(ns-1,j) = 1.0/speciesMw(j);

            // Reset yj
            rhoi[j] = rhoj;
        }

        cout << "jac = \n" << jac << endl;
        cout << "rank = " << jac.jacobiSvd().rank() << endl;

        // Compute update direction
        d = jac.colPivHouseholderQr().solve(f);
        rhoit = (rhoi - d).array().max(0.0);

        // Update state and recompute f
        setState(rhoit.data(), Tvec.data(), 1);
        dp = (xip - Map<const VectorXd>(X(), nSpecies())) / dx;
        smbf(Tw, Pw, dp, fp);

        // Armijo line search
        double lambda = 1.0;
        while (fp.norm() > (1.0 - alpha*lambda)*f.norm()) {
            lambda *= 0.5;
            rhoit = (rhoi - lambda*d).array().max(0.0);

            // Update state and recompute f
            setState(rhoit.data(), Tvec.data(), 1);
            dp = (xip - Map<const VectorXd>(X(), nSpecies())) / dx;
            smbf(Tw, Pw, dp, fp);
        }

        f = fp;
        rhoi = rhoit;
        cout << "norm(f) = " << f.norm() << endl;
    }

    return 0.0;
}

void Mixture::smbf(double Tw, double Pw, const VectorXd& dx, VectorXd& f)
{
    const double ratio = 1.0;
    Map<const ArrayXd> yiw(Y(), nSpecies());

    // Compute diffusion fluxes
    double E; stefanMaxwell(dx.data(), f.data(), E);
    f.array() *= yiw * density();

    // Compute the surface source terms
    ArrayXd mdotci(nSpecies()); computeSurfaceSource(mdotci);
    double mdotc = mdotci.sum();
    double mdotg = (ratio-1.0)*mdotc;

    //cout << "mdotci = \n";
    //for (int i = 0; i < nSpecies(); ++i)
    //    cout << setw(20) << speciesName(i) << " " << mdotci[i] << endl;

    // Compute the mass fractions of pyrolysis gas
    ArrayXd xeg(nElements()); getComposition("pyro", xeg.data());
    ArrayXd yig(nSpecies());
    equilibriumComposition(Tw, Pw, xeg.data(), yig.data());
    convert<X_TO_Y>(yig.data(), yig.data());

    // Now put it all together to get the function
    f.array() += (ratio*mdotc)*yiw - mdotci - mdotg*yig;
    f[nSpecies()-1] = density()/mixtureMw() - Pw/(Tw*RU);
}

void Mixture::computeSurfaceSource(ArrayXd& mci)
{
    const int ns = nSpecies();
    const int nh = nHeavy();
    const int ne = nElements();

    mci.fill(0.0);
    double rho = density();

    // Ion recombination
    if (ns > nh) {
        for (int i = 1, j; i < ns; ++i) {
            // If this is an ion, find parent neutral
            if (species()[i].charge() != 0) {
                for (j = 1; j < ns; ++j) {
                    if (i == j) continue;
                    if (elementMatrix().row(j).tail(ne-1) ==
                        elementMatrix().row(i).tail(ne-1)) break;
                }

                if (j == ns)
                    cout << "Warning: ion without parent neutral in smb." << endl;
                else {
                    //cout << "ion = " << speciesName(i) << ", neutral = " << speciesName(j) << endl;
                    // Add mass flux contributions
                    double mwi = speciesMw(i);
                    double flux = std::sqrt(RU*T()/(TWOPI*mwi))*rho*Y()[i];
                    mci[i] -= flux; // ion
                    mci[j] += flux*speciesMw(j)/mwi; // neutral
                }
            }
        }

        // Electrons go away
        mci[0] = -std::sqrt(RU*T()/TWOPI*speciesMw(0))*rho*Y()[0];
    }

    // Surface reactions
    int iO  = speciesIndex("O");
    int iO2 = speciesIndex("O2");
    int iCO = speciesIndex("CO");
    int iN  = speciesIndex("N");
    int iCN = speciesIndex("CN");
    int iC3 = speciesIndex("C3");

    double mdot, fac = std::sqrt(RU*T()/TWOPI)*rho;
    // O + C(s) -> CO
    mdot = fac/std::sqrt(speciesMw(iO))*Y()[iO]*0.63*std::exp(-1160.0/T());
    mci[iO]  += mdot;
    mci[iCO] -= mdot*speciesMw(iCO)/speciesMw(iO);

    // O2 + 2C(s) -> 2CO
    mdot = fac/std::sqrt(speciesMw(iO2))*Y()[iO2];
    mci[iO2] += mdot;
    mci[iCO] -= 2.0*speciesMw(iCO)/speciesMw(iO2);

    // N + C(s) -> CN
    mdot = fac/std::sqrt(speciesMw(iN))*Y()[iN]*3.0e-3;
    mci[iN]  += mdot;
    mci[iCN] -= mdot*speciesMw(iCN)/speciesMw(iN);

    // 3C(s) -> C3
    double yeq = 5.19e14*std::exp(-90845.0/T())*speciesMw(iC3)/(RU*T()*rho);
    mci[iC3] -= fac/std::sqrt(speciesMw(iC3))*(yeq - Y()[iC3]);

    //mci.fill(0.0);
}

//==============================================================================

} // namespace Mutation
