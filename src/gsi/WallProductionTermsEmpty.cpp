/**
 * @file WallProductionTermsEmpty.cpp
 *
 * @brief Class which should be used when no surface processes occur,
 */

/*
 * Copyright 2014-2018 von Karman Institute for Fluid Dynamics (VKI)
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


#include "AutoRegistration.h"
#include "Transport.h"
#include "Utilities.h"

#include "WallProductionTerms.h"

using namespace Eigen;

using namespace Mutation::Utilities;

namespace Mutation {
    namespace GasSurfaceInteraction {

class WallProductionsTermsEmpty : public WallProductionTerms
{
public:
    WallProductionsTermsEmpty(ARGS args)
                         : WallProductionTerms(args),
                           m_tag("empty") {}

//==============================================================================

    ~WallProductionsTermsEmpty(){}

//==============================================================================

    void productionRate(VectorXd& v_mass_prod_rate){
        v_mass_prod_rate.setZero();
    }

//==============================================================================

    const std::string& getWallProductionTermTag() const { return m_tag; }

private:
    const std::string m_tag;

};

Config::ObjectProvider<
    WallProductionsTermsEmpty, WallProductionTerms>
    wall_productions_terms_empty("empty");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
