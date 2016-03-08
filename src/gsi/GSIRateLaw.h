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

    virtual double forwardReactionRate( const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Twall ) const = 0;

//=============================================================================================================

protected:
    const Mutation::Thermodynamics::Thermodynamics& m_thermo;

    inline double computeAverageThermalSpeedforSpeciesI( const int& l_index_species, const Eigen::VectorXd& v_Twall ) const {

    	int pos_T_trans = 0;
    	double l_T_trans = v_Twall(pos_T_trans);
        return sqrt( ( 8.E0 * Mutation::RU * l_T_trans ) / ( Mutation::PI * m_thermo.speciesMw( l_index_species ) ) ) ;

    }

//=============================================================================================================

};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // GSIRATELAW_H
