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

#ifndef TEST_MACROS_H
#define TEST_MACROS_H

/**
 * Loops over mixtures, loads each mixture and runs __CODE__ in Catch SECTION
 * which is the name of the mixture.
 */
#define MIXTURE_LOOP(__CODE__)\
std::string _names_[10] = {\
    "air5_RRHO_ChemNonEq1T",\
    "air5_RRHO_ChemNonEqTTv",\
    "air5_NASA-7_ChemNonEq1T",\
    "air5_NASA-9_ChemNonEq1T",\
    "air11_RRHO_ChemNonEq1T",\
    "air11_RRHO_ChemNonEqTTv",\
    "air11_NASA-7_ChemNonEq1T",\
    "air11_NASA-9_ChemNonEq1T",\
    "argon_CR_ChemNonEq1T",\
    "argon_CR_ChemNonEqTTv"\
};\
Mutation::GlobalOptions::workingDirectory(TEST_DATA_FOLDER);\
for (int i = 0; i < 10; ++i) {\
    SECTION(_names_[i]) {\
        Mixture mix(_names_[i]);\
        __CODE__ ;\
    }\
}


/**
 * Loops over mixtures for Gas-Surface Interaction Mass Balance,
 * loads each mixture and runs __CODE__ in Catch SECTION
 * which is the name of the mixture.
 */
#define MIXTURE_GSI_MASS_LOOP(__CODE__)\
std::string _names_[5] = {\
    "smb_air5_RRHO_ChemNonEq1T",\
    "smb_oxidation_NASA9_ChemNonEq1T",\
    "smb_sublimation_NASA9_ChemNonEq1T",\
    "smb_aircarbon11_RRHO_ChemNonEq1T",\
    "smb_o2_RRHO_ChemNonEq1T"\
};\
Mutation::GlobalOptions::workingDirectory(TEST_DATA_FOLDER);\
for (int i = 0; i < 5; ++i) {\
    SECTION(_names_[i]) {\
        Mixture mix(_names_[i]);\
        __CODE__ ;\
    }\
}

/**
 * Loops over 10 temperatures and 10 pressures between 1000-10000 K (linear) and
 * 10-100000 Pa (log) and equilibrates the mixture at each (T,P) combination
 * before running the given code snippet.
 */
#define EQUILIBRATE_LOOP(__CODE__)\
for (int ip = 0.0; ip < 10; ++ip) {\
    double P = std::exp(ip/9.0*std::log(100000.0)+std::log(10.0));\
    for (int it = 0.0; it < 10; ++it) {\
        double T = 1000.0*it + 1000.0;\
        mix.equilibrate(T, P);\
        __CODE__ ;\
    }\
}

#endif // TEST_MACROS_H
