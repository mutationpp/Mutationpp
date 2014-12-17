/**
 * @file StateModel.h
 *
 * @brief Implementation of the StateModel class.
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

#include "StateModel.h"
#include "TransferModel.h"
#include "Mixture.h"

using namespace Mutation::Transfer;
using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace Thermodynamics {

//==============================================================================
StateModel::StateModel(ARGS mix, const int nenergy, const int nmass)
	: m_mix(mix), m_ns(mix.nSpecies()), m_nenergy(nenergy), m_nmass(nmass)
{
	m_T = m_Tr = m_Tv = m_Tel = m_Te = 300.0;
	m_P = 0.0;
	mp_X = new double [m_ns];
	for (int i = 0; i < m_ns; ++i)
		mp_X[i] = 0.0;
}

//==============================================================================
StateModel::~StateModel()
{
	// Delete data pointers
	delete [] mp_X;

	// Delete transfer models
	for (int i = 0; i < m_transfer_models.size(); ++i)
		delete m_transfer_models[i].second;
}

//==============================================================================
void StateModel::addTransferTerm(int i, const std::string& model)
{
	assert(i >= 0);  assert(i < m_nenergy-1);
	m_transfer_models.push_back(std::make_pair(i,
		Factory<TransferModel>::create(model, m_mix)));
}

//==============================================================================
void StateModel::energyTransferSource(double* const p_omega)
{
	for (int i = 0; i < m_nenergy-1; ++i)
		p_omega[i] = 0.0;

	for (int i = 0; i < m_transfer_models.size(); ++i)
		p_omega[m_transfer_models[i].first] +=
			m_transfer_models[i].second->source();
}

//==============================================================================
    } // namespace Thermodynamics
} // namespace Mutation
