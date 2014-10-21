/**
 * @file Mixture.h
 *
 * @copyright
 * Copyright 2014 James B. Scoggins
 *
 * @license
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

#ifndef MUTATION_MIXTURE_H
#define MUTATION_MIXTURE_H

#include "Thermodynamics.h"
#include "Kinetics.h"
#include "TransferModel.h"
#include "Transport.h"
#include "MixtureOptions.h"
#include "StateModel.h"

namespace Mutation {

/**
 * Packages all of the Mutation++ functionality pertaining to a mixture into a
 * single class.  A Mixture object simply inherits all functionalities from the
 * Thermodynamics, Transport, and Kinetics classes.  A Mixture object can be
 * constructed using a mixture file name or a MixtureOptions object.
 *
 * @see Thermodynamics
 * @see Transport
 * @see Kinetics
 */
class Mixture
    : public Thermodynamics::Thermodynamics, 
      public Transport::Transport, 
      public Kinetics::Kinetics
{
public:

    /**
     * Constructs a Mixture using the options given.  Note that a mixture file
     * name can also be given instead and the options will be loaded 
     * accordingly.
     *
     * @see MixtureOptions
     */
    Mixture(const MixtureOptions& options);
    
    /** 
     * Destructor.
     */
    ~Mixture(){}
//        if (mp_transfer != NULL) delete mp_transfer;
    
    /**
     * Provides energy transfer source terms based on the current state of the
     * mixture.
     */
    void energyTransferSource(double* const p_source) {
         state()->energyTransferSource(p_source);
    }

//private:

//    Transfer::TransferModel* mp_transfer;
    
}; // class Mixture

} // namespace Mutation

#endif // MUTATION_MIXTURE_H

