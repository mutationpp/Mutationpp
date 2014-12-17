/**
 * @file StateModelInterface.cpp
 *
 * @brief Implementation of the StateModelInterface class.
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

#include "StateModelInterface.h"
#include "StateModel.h"

namespace Mutation {
    namespace Thermodynamics {

//==============================================================================
StateModelInterface::~StateModelInterface()
{
	if (mp_model)
		delete mp_model;
}

//==============================================================================
void StateModelInterface::setStateModel(StateModel* const p_model)
{
	assert(p_model);
	if (mp_model != p_model) {
		if (mp_model) delete mp_model;
		mp_model = p_model;
	}
}

//==============================================================================
    } // namespace Thermodynamics
} // namespace Mutation
