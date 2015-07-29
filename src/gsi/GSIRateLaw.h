#ifndef GSIRATELAW_H
#define GSIRATELAW_H

#include "DataGSIRateLaw.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLaw {

public:
    typedef const DataGSIRateLaw& ARGS;

    GSIRateLaw( ARGS l_data_gsi_rate_law )
              : m_thermo( l_data_gsi_rate_law.s_thermo )
    { }
    virtual ~GSIRateLaw(){ }

//=============================================================================================================

    virtual double forwardReactionRate( const Mutation::Numerics::RealVector& v_rhoi, const Mutation::Numerics::RealVector& v_Twall ) const = 0;

//=============================================================================================================

protected:
    const Mutation::Thermodynamics::Thermodynamics& m_thermo;

    inline double computeAverageThermalSpeedforSpeciesI( const int& l_index_species, const Mutation::Numerics::RealVector& v_Twall ) const {

        return sqrt( ( 8.E0 * Mutation::RU * v_Twall(0) ) / ( Mutation::PI * m_thermo.speciesMw( l_index_species ) ) ) ; /** @todo 0 -> position_translational_temperature */

    }

//=============================================================================================================

};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // GSIRATELAW_H
