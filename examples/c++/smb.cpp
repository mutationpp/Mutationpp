/**
 * @file smb.cpp
 *
 * @brief Surface Mass Balance example program.
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


// Must include this header file to use the Mutation++ library
#include "mutation++.h"
#include <iostream>
#include <iomanip>
using namespace std;
#include <Eigen/Dense>
using namespace Eigen;

using namespace Mutation;
using namespace Mutation::Thermodynamics;

int main()
{
    // Generate the default options for the air11 mixture
    Mixture mix("Johnston22");
    double Tw = 2500.0;
    double Pw = 10000.0;

    // Compute equilibrium air composition
    mix.equilibrate(Tw, Pw);
    VectorXd rhoi(mix.nSpecies());
    double rho = mix.density();
    for (int i = 0; i < rhoi.size(); ++i)
        rhoi[i] = rho * mix.Y()[i];

    // Physical cell next to wall will be at similar conditions
    mix.equilibrate(Tw + 500.0, Pw);
    VectorXd xip(mix.nSpecies());
    for (int i = 0; i < rhoi.size(); ++i)
        xip[i] = mix.X()[i];

    mix.smb(Tw, Pw, rhoi.data(), xip.data(), 0.0001);
    
    cout << "xip = \n";
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(20) << mix.speciesName(i) << " " << xip[i] << endl;

    cout << "xiw = \n";
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(20) << mix.speciesName(i) << " " << mix.X()[i] << endl;

    cout << "mdotci = \n";
    ArrayXd mdotci(mix.nSpecies());
    mix.computeSurfaceSource(mdotci);
    cout << mdotci << endl;
    cout << mdotci.sum() << endl;
    
    return 0;
}
/// [example_code]


