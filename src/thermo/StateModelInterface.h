/**
 * @file StateModelInterface.h
 *
 * @brief Definition of the StateModelInterface class.
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

#ifndef THERMO_STATE_MODEL_INTERFACE_H
#define THERMO_STATE_MODEL_INTERFACE_H

namespace Mutation {
    namespace Thermodynamics {

class StateModelInterface
{
public:
    StateModelInterface()
        : mp_model(0)

    void setStateModel(StateModel* const p_model) {
        mp_model = p_model;
    }

    StateModel* const stateModel() {
        assert(mp_model);
        return mp_model;
    }

private:
    SharedPointer<StateModel> mp_model;
};

    } // namespace Thermodynamics
} // namespace Mutation

#endif // THERMO_STATE_MODEL_INTERFACE_H
