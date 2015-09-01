#ifndef DIFFUSIONDRIVINGFORCES_H
#define DIFFUSIONDRIVINGFORCES_H

namespace Mutation {
    namespace GasSurfaceInteraction {

class DiffusionDrivingForces {
public:
    DiffusionDrivingForces( const Mutation::Thermodynamics::Thermodynamics& l_thermo )
                          : v_mole_frac_edge( l_thermo.nSpecies() ),
                            m_diff_model_set( 0 ) { }
    ~DiffusionDrivingForces(){ }

    void computeDrivingForces( const Mutation::Numerics::RealVector& lv_mole_frac, Mutation::Numerics::RealVector& lv_driving_force ){ 


        if ( m_diff_model_set == 0 ) { std::cerr << "Error diffusion model!" << std::endl; exit(1); }

        for ( int i_ns = 0 ; i_ns < lv_mole_frac.size() ; ++i_ns ){
            lv_driving_force( i_ns ) = ( lv_mole_frac( i_ns ) - v_mole_frac_edge( i_ns ) ) / m_dx;
        }

    }

    void setDiffusionCalculator( const Mutation::Numerics::RealVector& lv_mole_frac_edge, const double& l_dx ){

            v_mole_frac_edge = lv_mole_frac_edge; 
            m_dx = l_dx;

            m_diff_model_set = 1;

    }

private:
    Mutation::Numerics::RealVector v_mole_frac_edge;
    double m_dx;

    bool m_diff_model_set;

};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // DIFFUSIONDRIVINGFORCES_H
