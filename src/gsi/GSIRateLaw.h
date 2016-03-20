#ifndef GSIRATELAW_H
#define GSIRATELAW_H

// JB: why not simply include all the DataX classes directly in the header
// of the class that uses them?  This would remove 7 extra header files and
// improve the readability of X.h headers because you wouldn't have to go look
// for what is included in DataX.h.  Further more, DataX is only used by the X
// class and any classes that use X, so you don't save anything by placing
// DataX in another header.
#include "DataGSIRateLaw.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

// JB: Please comment all classes with Doxygen style comment
class GSIRateLaw {

public:
    typedef const DataGSIRateLaw& ARGS;

    // JB: avoid using unnecessarily long argument names, (ARGS args) works just fine
    GSIRateLaw( ARGS l_data_gsi_rate_law )
              : m_thermo( l_data_gsi_rate_law.s_thermo )
    { }
    virtual ~GSIRateLaw(){ }

//=============================================================================================================

    // JB: please split long function declarations on multiple lines like this
    virtual double forwardReactionRate(
        const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Twall) const = 0;

//=============================================================================================================

protected:

    const Mutation::Thermodynamics::Thermodynamics& m_thermo;

    // JB: this function doesn't really belong in the GSIRateLaw class, I think it could be moved to
    // Transport (if it isn't already there...)
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
