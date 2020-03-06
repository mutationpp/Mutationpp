/**
 * @file TransferModel.h
 *
 * @brief Energy Transfer Models
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

#ifndef TRANSFER_TRANSFER_MODEL_H
#define TRANSFER_TRANSFER_MODEL_H

namespace Mutation {

	// Forward declaration of Mixture type
    class Mixture;

    namespace Transfer {

/**
 * @defgroup transfermodels Energy Transfer Models
 * @{
 */

/**
 * Base class for all energy transfer models.  Provides the abstract interface
 * which must be implemented by concrete models.  Also enables self-registration
 * of TransferModel objects.
 */
class TransferModel
{
public:
    // Allows self-registration of TransferModel objects
    typedef Mutation::Mixture& ARGS;

/**
 *@brief Returns name of this type.
 */

    static std::string typeName() { return "TransferModel"; }

    TransferModel(ARGS mix)
        : m_mixture(mix)
    { }

/**
 *@brief Transfer Model destructor
 */
    virtual ~TransferModel() { }

/**
 *@brief Purely virtual function to be called to
 * to solve the source term
 *
 * @return energy source term
 */
    virtual double source() = 0;

protected:
    Mutation::Mixture& m_mixture;
};

/// @}

    } // namespace Transfer
} // namespace Mutation

#endif // TRANSFER_TRANSFER_MODEL_H
