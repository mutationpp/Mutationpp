/**
 * @file MillikanWhite.h
 *
 * @brief Declaration of classes related to Millikan and White model.
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


#ifndef MILLIKAN_WHITE_H
#define MILLIKAN_WHITE_H

#include <Eigen/Dense>

#include <cmath>
#include <memory>
#include <string>
#include <vector>

namespace Mutation { namespace Thermodynamics { class Thermodynamics; } }

namespace Mutation {
    namespace Transfer {

/**
 * Convenience class for building necessary data used in the Millikan-White
 * vibration-translation energy relaxation model.
 */
class MillikanWhiteModelData
{
public:
    /**
     * Constructs default data for Millikan-White model from thermo model,
     * species index, and characteristic vibrational temperature.  Default
     * values for the A and B parameters in the Millikan-White model are taken
     * from the original formulas in the paper:
     * 
     * \f{align*}{
     * &A_{ij} = 1.16e-3 \mu_{ij}^{1/2} \theta_{i}^{4/3} \\
     * &B_{ij} = 0.015 \mu_{ij}^{1/4}
     * \f}
     * 
     * where \f$\mu_{ij}\f$ is the reduced mass of the vibrator and colliding
     * partner in g/mol, and \f$\theta_{i}\f$ is the characteristic vibrational
     * temperature of the vibrator.  The default value for Park's reference 
     * cross section is taken as 10^-16 cm^2, based on recommendation in
     * @cite Gnoffo1989. The high-temperature limiting cross section is based
     * on @cite Park1993.
     */
    MillikanWhiteModelData(
        const class Thermodynamics::Thermodynamics&, size_t, double);

    /// Copy constructor.
    MillikanWhiteModelData(const MillikanWhiteModelData& data);

    /// Copy assignment operator.
    MillikanWhiteModelData& operator=(MillikanWhiteModelData data);

    /// Move constructor.
    MillikanWhiteModelData(MillikanWhiteModelData&& data) noexcept;

    /// Move assignment operator.
    MillikanWhiteModelData& operator=(MillikanWhiteModelData&& data) noexcept;

    /// Destructor.
    ~MillikanWhiteModelData();

    /// Index of species in mixture.
    size_t speciesIndex() const;

    /// Molecular weight of the vibrator in kg/mol
    double molecularWeight() const;

    /// Returns the array of A model parameters.
    Eigen::ArrayXd& a();

    /// Returns the array of A model parameters.
    const Eigen::ArrayXd& a() const;

    /// Returns the array of B model parameters.
    Eigen::ArrayXd& b();

    /// Returns the array of B model parameters.
    const Eigen::ArrayXd& b() const;

    /// Sets the reference cross section in m^2.
    MillikanWhiteModelData& setReferenceCrossSection(double omegav);

    /// Returns the reference cross section for Park's correction in m^2.
    double referenceCrossSection() const;

    /// Returns the limiting cross section for Park's correction in m^2.
    double limitingCrossSection(const double& T) const;

private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};


/**
 * @brief Millikan-White relaxation model with Park correction.
 * 
 * Implements the Millikan and White relaxation time model proposed in 
 * @cite Millikan1963 with the high temperature correction of Park @cite Park1993.
 * For an overview, see Eqs. (55-58) in @cite Gnoffo1989.
 */
class MillikanWhiteModel
{
public:
    /// Constructs a MillikanWhiteModel from data
    MillikanWhiteModel(const MillikanWhiteModelData&);

    /// Index of species in mixture.
    size_t speciesIndex() const { return m_data.speciesIndex(); }

    /// Molecular weight of the vibrator in kg/mol
    double molecularWeight() const { return m_data.molecularWeight(); }

    /// Relaxation time of the vibrator.
    double relaxationTime(const Mutation::Thermodynamics::Thermodynamics&) const;

    const MillikanWhiteModelData& data() const { return m_data; }

private:
    MillikanWhiteModelData m_data;
};


/// Builder for MillikanWhiteModel based on XML database.
class MillikanWhiteModelDB
{
public:
    /// Constructs the builder using default file path for "VT.xml".
    MillikanWhiteModelDB(
        const Mutation::Thermodynamics::Thermodynamics&);

    /// Constructs the builder using given file path for database.
    MillikanWhiteModelDB(
        const Mutation::Thermodynamics::Thermodynamics&, std::string file_path);

    /// Creates a MillikanWhiteModel given the species name and characteristic
    /// vibrational temperature in K.
    MillikanWhiteModel create(std::string species_name, double theta);

private:
    struct Data;
    std::shared_ptr<Data> m_data;
};

    } // namespace Kinetics
} // namespace Mutation

#endif // MILLIKAN_WHITE_H
