#include "AutoRegistration.h"
#include "Utilities.h"

#include "GSIRateLaw.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLawGammaConst : public GSIRateLaw {

public:
    GSIRateLawGammaConst( ARGS l_data_gsi_rate_law )
                        : GSIRateLaw ( l_data_gsi_rate_law ),
                          v_reactants( l_data_gsi_rate_law.s_reactants ) {

        assert( l_data_gsi_rate_law.s_node_rate_law.tag() == "gamma_const" );

        if ( v_reactants.size() > 2 ) { /** errorMoreThanTwoReactants @todo */; }

        int l_diff_reac = 1;
        for ( int i_reac = 0; i_reac < v_reactants.size() - 1; ++i_reac ){
            if (v_reactants[i_reac] != v_reactants[i_reac + 1]) l_diff_reac++;
        }
    
        std::vector<std::string> tokens;
        Mutation::Utilities::String::tokenize( l_data_gsi_rate_law.s_node_rate_law.text(), tokens, ":, " );
        int index_ref = -2;
        for (int i = 0; i < tokens.size(); i+=2) {
            int index = m_thermo.speciesIndex(tokens[i]);
            if(index_ref > index){
                v_gamma.insert(v_gamma.begin(), atof(tokens[i+1].c_str()));
            } else {
                v_gamma.insert(v_gamma.end(), atof(tokens[i+1].c_str()));
            }
            index_ref = index;
        } /** @todo: test **/
    
        if (l_diff_reac != v_gamma.size()){
            std::cerr << "ERROR @ GammaModelConst @ CatalysisRateLaws.cpp. #gammas should be = #different reactants." << std::endl; /** @todo cout better error */
            exit(1);
        }

        v_output_impinging_flux.resize(v_gamma.size());
        v_impinging_flux_per_stoichiometric_coefficient.resize(v_gamma.size());

    }

//=============================================================================================================

    ~GSIRateLawGammaConst( ){ }

//=============================================================================================================

    double forwardReactionRate( const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Twall ) const {

        l_index_v_reactants = 0;
    
        for( int i_gammas = 0 ; i_gammas < v_gamma.size() ; i_gammas++ ){
    
            getSpeciesIndexandStoichiometricCoefficient( l_index_v_reactants, l_species_index, l_stoichiometric_coefficient );
            v_output_impinging_flux[ i_gammas ] = computeWallImpingingMassFlux( l_species_index, v_rhoi, v_Twall );
    
            v_impinging_flux_per_stoichiometric_coefficient[ i_gammas ] = v_output_impinging_flux[ i_gammas ] / l_stoichiometric_coefficient;
            v_output_impinging_flux[ i_gammas ] = v_impinging_flux_per_stoichiometric_coefficient[ i_gammas ] * v_gamma[ i_gammas ];
    
            l_index_v_reactants += l_stoichiometric_coefficient;
    
        }
    
        return getLimitingImpingingMassFlux();

    }

//=============================================================================================================

private:
    mutable int l_index_v_reactants; 
    mutable int l_species_index;
    mutable int l_stoichiometric_coefficient;

    std::vector<double> v_gamma;

    mutable std::vector<double> v_output_impinging_flux;
    mutable std::vector<double> v_impinging_flux_per_stoichiometric_coefficient;

    const std::vector<int>& v_reactants;

//=============================================================================================================

    inline void getSpeciesIndexandStoichiometricCoefficient( int l_index_v_reactants, int& l_species_index, int& l_stoichiometric_coefficient ) const { 

        l_species_index = v_reactants[ l_index_v_reactants ];
        l_stoichiometric_coefficient = 1;
        l_index_v_reactants++;
        
        while( l_index_v_reactants < v_reactants.size() ){
            if ( l_species_index != v_reactants[ l_index_v_reactants ] ) {
                break;
            } 
            l_stoichiometric_coefficient++;
            l_index_v_reactants++;
        }

    }

//=============================================================================================================

    inline double computeWallImpingingMassFlux( const int& l_index_species, const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Twall ) const {
        
        return computeAverageThermalSpeedforSpeciesI( l_index_species, v_Twall ) / ( 4.E0 )
        		* v_rhoi( l_index_species ) / m_thermo.speciesMw( l_index_species );

    }

//=============================================================================================================

    inline double getLimitingImpingingMassFlux() const {

        return v_output_impinging_flux[ min_element( v_impinging_flux_per_stoichiometric_coefficient.begin(),
        		v_impinging_flux_per_stoichiometric_coefficient.end()) - v_impinging_flux_per_stoichiometric_coefficient.begin() ];

    }

//=============================================================================================================

};

Mutation::Utilities::Config::ObjectProvider<GSIRateLawGammaConst, GSIRateLaw> gsi_rate_law_gamma_const("gamma_const");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
