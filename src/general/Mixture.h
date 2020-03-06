/**
 * @file Mixture.h
 *
 * @brief Provides Mixture class declaration. @see Mutation::Mixture
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

#ifndef MUTATION_MIXTURE_H
#define MUTATION_MIXTURE_H

#include <vector>

#include "Thermodynamics.h"
#include "Kinetics.h"
#include "Transport.h"
#include "MixtureOptions.h"
#include "StateModel.h"
#include "Composition.h"
#include "GasSurfaceInteraction.h"

namespace Mutation {

/**
 * Packages all of the Mutation++ functionality pertaining to a mixture into a
 * single class.  A Mixture object simply inherits all functionalities from the
 * Thermodynamics, Transport, and Kinetics classes.  A Mixture object can be
 * constructed using a mixture file name or a MixtureOptions object.
 *
 * @see Thermodynamics::Thermodynamics
 * @see Transport::Transport
 * @see Kinetics::Kinetics
 */
class Mixture
    : public Thermodynamics::Thermodynamics, 
      public Transport::Transport, 
      public Kinetics::Kinetics,
      public GasSurfaceInteraction::GasSurfaceInteraction
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
    ~Mixture() { }
    
    /**
     * Provides energy transfer source terms based on the current state of the
     * mixture.
     */
    void energyTransferSource(double* const p_source) {
         state()->energyTransferSource(p_source);
    }

    /**
     * Add a named element composition to the mixture which may be retrieved
     * with getComposition().
     */
    void addComposition(
        const Mutation::Thermodynamics::Composition& c,
        bool make_default = false);

    /**
     * Gets the element mole or mass fractions associated with a named
     * composition in the mixture.  Returns false if the composition does not
     * exist in the list, true otherwise.
     */
    bool getComposition(
        const std::string& name, double* const p_vec,
        Mutation::Thermodynamics::Composition::Type type =
            Mutation::Thermodynamics::Composition::MOLE) const;


private:

    std::vector<Mutation::Thermodynamics::Composition> m_compositions;

}; // class Mixture

} // namespace Mutation

#endif // MUTATION_MIXTURE_H

