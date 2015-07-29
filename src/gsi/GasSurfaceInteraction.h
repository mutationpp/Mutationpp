#ifndef GASSURFACEINTERACTION_H
#define GASSURFACEINTERACTION_H

#include <string>

#include "Thermodynamics.h"
#include "Transport.h"

#include "SurfaceDescription.h" 
#include "SurfaceBalanceSolver.h" 

namespace Mutation {
    namespace GasSurfaceInteraction {

class GasSurfaceInteraction {

public:
    GasSurfaceInteraction( Mutation::Thermodynamics::Thermodynamics& l_thermo, Mutation::Transport::Transport& l_transport, std::string l_gsi_mechanism_file );
    ~GasSurfaceInteraction();

    void setWallState( const double* const l_mass, const double* const l_energy, const int state_variable );

    void surfaceProductionRates( double * const lp_mass_prod_rate );
    
    void setDiffusionModel( const double* const lp_mole_frac_edge, const double& l_dx );
    void solveSurfaceBalance();

private:
    Mutation::Thermodynamics::Thermodynamics& m_thermo;
    Mutation::Transport::Transport& m_transport;

    SurfaceDescription* mp_surf_descr;
    SurfaceBalanceSolver* mp_surf_solver;

    std::string m_gsi_mechanism;

    Mutation::Numerics::RealVector v_mass_prod_rate;

    inline void locateGSIInputFile( std::string& l_gsi_mechanism_file );
    inline void errorWrongTypeofGSIFile( const std::string& l_gsi_root_tag );
    inline void getXmlPositionPointerSurfPropsDiffModelProdTerm( Mutation::Utilities::IO::XmlElement::const_iterator& l_index_surf_descr, Mutation::Utilities::IO::XmlElement::const_iterator& l_index_diff_model, Mutation::Utilities::IO::XmlElement::const_iterator& l_index_prod_terms, const Mutation::Utilities::IO::XmlElement& root );
    inline void errorInvalidGSIFileProperties( const std::string& l_gsi_option );

    // TEMPORARY
    Mutation::Numerics::RealVector v_mole_frac_edge;

};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // GASSURFACEINTERACTION_H
