/**
 * @file equilibrium_air.cpp
 *
 * @brief Equilibrium Air example program.
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

/**
 * @page example_equilibrium_air Equilibrium Air
 *
 * @tableofcontents
 *
 * @section equil_air_example_intro Introduction
 * This example program makes use of the Mutation++ library to compute
 * equilibrium mole fractions, mixture frozen specific heat, enthalpy, and
 * entropy for an air mixture containing 20% O and 80% N at a pressure of 1 atm
 * and a range of temperatures from 300 to 6,000 K using the NASA 9-coefficient
 * polynomial thermodynamic database.  The output is a property table written to
 * standard output as a function of temperature.
 *
 * The example demonstrates the following features:
 * - Using the Mutation::MixtureOptions class to modify the basic options in a
 *   [mixture](@ref mixtures)
 * - Adding a default elemental composition to a Mutation::Mixture object with
 *   Mutation::Mixture::addComposition()
 * - Using the [EquilTP](@ref Mutation::Thermodynamics::EquilTPStateModel)
 *   state model for thermochemical equilibrium mixtures
 * - Setting the state of a mixture
 * - Retrieving various properties of a mixture
 *
 * @section equil_air_example_code Example Code
 * @snippet examples/c++/equilibrium_air.cpp example_code
 *
@section equil_air_example_output Expected Output
@verbatim
  T[K]          X_N          X_O         X_NO         X_N2         X_O2   Cp[J/kg-K]      H[J/kg]    S[J/kg-K]
    300  4.97208e-80  2.21196e-41  2.54688e-16          0.8          0.2       1012.1      1872.32      6888.78
    400  1.87774e-59  1.60532e-30  2.34838e-12          0.8          0.2      1019.76       103405      7180.82
    500  4.32695e-47  5.46135e-24  5.61479e-10          0.8          0.2      1034.71       206072      7409.85
    600  7.73459e-39   1.2622e-19  2.16196e-08          0.8          0.2      1055.22       310534      7600.25
    700  6.16772e-33  1.67571e-16  2.93277e-07          0.8          0.2      1078.31       417202      7764.64
    800  1.66538e-28  3.72075e-14  2.07298e-06     0.799999     0.199999      1101.53       526206      7910.16
    900  4.69776e-25  2.50363e-12  9.48874e-06     0.799995     0.199995      1123.39       637491      8041.21
   1000   2.7227e-22  7.29609e-11  3.20434e-05     0.799984     0.199984      1143.19       750910      8160.68
   1100  4.98744e-20  1.15605e-09  8.67361e-05     0.799957     0.199957      1160.75       866298      8270.65
   1200  3.84849e-18  1.15917e-08  0.000198877     0.799901     0.199901      1176.15       983515      8372.63
   1300   1.5264e-16  8.17101e-08  0.000401321     0.799799     0.199799      1189.62  1.10246e+06      8467.83
   1400  3.58744e-15  4.36527e-07  0.000732423     0.799634     0.199634       1201.4  1.22307e+06       8557.2
   1500   5.5458e-14  1.86768e-06    0.0012333     0.799383     0.199382      1211.74  1.34534e+06      8641.55
   1600  6.09844e-13  6.66995e-06   0.00194494     0.799025     0.199024      1220.86  1.46927e+06      8721.53
   1700  5.06489e-12  2.05218e-05   0.00290559     0.798539     0.198535      1228.94  1.59493e+06       8797.7
   1800  3.32865e-11  5.57536e-05    0.0041486     0.797903     0.197892      1236.16  1.72244e+06      8870.58
   1900    1.796e-10  0.000136373   0.00570074     0.797095     0.197068      1242.63    1.852e+06      8940.63
   2000  8.19371e-10  0.000305015   0.00758099     0.796087     0.196026      1248.47  1.98401e+06      9008.33
   2100   3.2371e-09  0.000631657    0.0097995     0.794848     0.194721      1253.78  2.11904e+06       9074.2
   2200  1.12929e-08   0.00122357    0.0123566     0.793332     0.193088      1258.61  2.25801e+06      9138.85
   2300  3.53538e-08   0.00223545    0.0152414     0.791485     0.191038      1263.05  2.40224e+06      9202.95
   2400  1.00663e-07   0.00387833    0.0184305     0.789233     0.188458      1267.15  2.55353e+06      9267.32
   2500  2.63635e-07   0.00642503    0.0218851     0.786487     0.185202      1270.95  2.71421e+06       9332.9
   2600  6.41152e-07    0.0102097    0.0255491     0.783141       0.1811       1274.5   2.8872e+06      9400.73
   2700  1.45968e-06    0.0156179    0.0293462     0.779078     0.175956      1277.85  3.07589e+06      9471.93
   2800  3.13254e-06    0.0230645    0.0331777     0.774183     0.169572      1281.04  3.28407e+06      9547.61
   2900  6.37482e-06    0.0329551    0.0369218     0.768351     0.161765      1284.11   3.5156e+06      9628.83
   3000  1.23657e-05    0.0456316    0.0404363     0.761518     0.152402      1287.09  3.77399e+06       9716.4
   3100  2.29681e-05    0.0613007    0.0435635     0.753677     0.141435      1290.01  4.06185e+06      9810.76
   3200  4.10141e-05    0.0799578    0.0461425     0.744909      0.12895      1292.89  4.38015e+06      9911.79
   3300  7.06659e-05     0.101316    0.0480244     0.735398     0.115191      1295.76  4.72742e+06      10018.6
   3400  0.000117865     0.124772    0.0490931     0.725439     0.100579      1298.61  5.09914e+06      10129.6
   3500  0.000190878     0.149422    0.0492862     0.715417    0.0856849      1301.43  5.48752e+06      10242.2
   3600  0.000300966     0.174165     0.048611     0.705758    0.0711652      1304.21  5.88206e+06      10353.3
   3700   0.00046317     0.197873    0.0471498     0.696859    0.0576551      1306.91  6.27109e+06      10459.9
   3800  0.000697231     0.219579    0.0450484     0.689017    0.0456587      1309.51  6.64383e+06      10559.3
   3900   0.00102859     0.238639    0.0424907     0.682374    0.0354687      1311.97  6.99248e+06      10649.9
   4000   0.00148945     0.254782    0.0396671     0.676913    0.0271482       1314.3  7.31341e+06      10731.2
   4100   0.00211986      0.26807    0.0367463     0.672491    0.0205729      1316.51  7.60702e+06      10803.7
   4200   0.00296883      0.27878    0.0338587     0.668887    0.0155055       1318.6   7.8768e+06      10868.7
   4300   0.00409528     0.287293    0.0310934      0.66585     0.011668      1320.61  8.12794e+06      10927.8
   4400    0.0055693     0.293998     0.028502     0.663137   0.00879327      1322.58  8.36629e+06      10982.6
   4500   0.00747252     0.299245    0.0261089     0.660522   0.00665143      1324.54  8.59768e+06      11034.6
   4600   0.00989947     0.303322    0.0239186     0.657802   0.00505758      1326.53  8.82762e+06      11085.1
   4700    0.0129577     0.306454    0.0219245     0.654794   0.00386935       1328.6  9.06125e+06      11135.4
   4800    0.0167676     0.308811    0.0201133     0.651328   0.00298003      1330.81  9.30343e+06      11186.4
   4900    0.0214626     0.310514    0.0184687     0.647244   0.00231082       1333.2  9.55886e+06        11239
   5000    0.0271883      0.31165    0.0169736     0.642384   0.00180409      1335.84  9.83219e+06      11294.2
   5100    0.0341008     0.312277    0.0156114     0.636593   0.00141776      1338.81  1.01281e+07      11352.8
   5200    0.0423636     0.312432    0.0143665     0.629717   0.00112117      1342.18  1.04516e+07      11415.6
   5300    0.0521444     0.312135    0.0132251     0.621603   0.00089185      1346.04  1.08075e+07      11483.4
   5400    0.0636139     0.311397    0.0121746     0.612102  0.000713304       1350.5  1.12013e+07        11557
   5500    0.0769405     0.310218     0.011204     0.601064  0.000573338      1355.68  1.16385e+07      11637.2
   5600    0.0922694     0.308597    0.0103042     0.588367  0.000462931      1361.69  1.21246e+07      11724.8
   5700     0.109741     0.306529   0.00946703     0.573889  0.000375297      1368.66  1.26656e+07      11820.5
   5800     0.129457      0.30401    0.0086859     0.557542  0.000305351      1376.74  1.32672e+07      11925.2
   5900     0.151482     0.301041   0.00795539     0.539272  0.000249238      1386.08   1.3935e+07      12039.3
   6000     0.175836     0.297628   0.00727106     0.519061  0.000204012      1396.81  1.46744e+07      12163.5
@endverbatim
*/


/// [example_code]
// Must include this header file to use the Mutation++ library
#include "mutation++.h"
#include <iostream>

using namespace Mutation;
using namespace Mutation::Thermodynamics;

int main()
{
    // Generate the default options for the air11 mixture
    MixtureOptions opts("air_5");
    opts.setStateModel("EquilTP");

    // Change from default thermodynamic database (RRHO) to NASA 9-coefficient
    // polynomial database
    opts.setThermodynamicDatabase("RRHO");

    // Load the mixture with the new options
    Mixture mix(opts);

    // Setup the default composition
    mix.addComposition("N:0.8, O:0.2", true);

    // Write a header line for the table
    std::cout << std::setw(7) << "T[K]";
    for (int i = 0; i < mix.nSpecies(); ++i)
        std::cout << std::setw(13) << "X_" + mix.speciesName(i);
    std::cout << std::setw(13) << "Cp[J/kg-K]";
    std::cout << std::setw(13) << "H[J/kg]";
    std::cout << std::setw(13) << "S[J/kg-K]";
    std::cout << std::endl;

    // Loop over range of temperatures and compute equilibrium values at 1 atm
    double P = ONEATM;
    for (int i = 0; i < 58; ++i) {
        // Compute the temperature
        double T = 300.0 + static_cast<double>(i) * 100.0;

        // Set the mixture state equal to the equilibrium state for the given
        // temperature and pressure
        mix.setState(&T, &P);

        // Temperature
        std::cout << std::setw(7) << mix.T();

        // Species mole fractions
        for (int j = 0; j < mix.nSpecies(); ++j)
            std::cout << std::setw(13) << mix.X()[j];

        // Other properties
        std::cout << std::setw(13) << mix.mixtureFrozenCpMass(); // Cp [J/kg-K]
        std::cout << std::setw(13) << mix.mixtureHMass();        // H  [J/kg]
        std::cout << std::setw(13) << mix.mixtureSMass();        // S  [J/kg-K]

        std::cout << std::endl;
    }

    return 0;
}
/// [example_code]


