/**
 * @file MatrixHelper.h
 *
 * @brief Provides some helper functions for matrices.
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

#ifndef TESTS_MATRIX_HELPER_H
#define TESTS_MATRIX_HELPER_H

#include <fstream>
#include <string>
#include <Eigen/Dense>

/**
 * Reads a matrix in from a file.
 */
Eigen::MatrixXd readMatrix(std::ifstream& infile)
{
    int cols = 0, rows = 0;
    double buff[((int) 1e6)];

    // Read numbers from file into buffer.
    while (!infile.eof()) {
        std::string line;
        getline(infile, line);
        if (line.empty()) break;

        int temp_cols = 0;
        std::stringstream stream(line);
        while(!stream.eof())
            stream >> buff[cols*rows+temp_cols++];

        if (temp_cols == 0)
            continue;

        if (cols == 0)
            cols = temp_cols;

        rows++;
    }

    // Populate matrix with numbers.
    Eigen::MatrixXd result(rows,cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result(i,j) = buff[cols*i+j];

    return result;
};

#endif // TESTS_MATRIX_HELPER_H
