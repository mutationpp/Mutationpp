#ifndef CATRATELAW_H
#define CATRATELAW_H

#include "Utilities.h"
#include "Thermodynamics.h"

namespace Mutation{
    namespace gsi{
  
  
/**
 * Abstract base class for all rate laws which allows owners such as class
 * Reaction to store any rate law polymorphically.
 */
class CatalysisRateLaw{
public:

    CatalysisRateLaw(const Mutation::Thermodynamics::Thermodynamics& thermo):m_thermo(thermo) { }

    virtual ~CatalysisRateLaw() { }
    virtual CatalysisRateLaw* clone() const = 0;
    
    virtual double forwardReactionRate(const double* const lp_rhoi, const double* const lp_Twall) const = 0; 

    /**
     * This function computes the average thermal speed in the species at a given temperature . It has units $[m/s]$.
     */
    inline double computeAverageThermalSpeed( const int& l_sp, const double* const lp_Twall ) const;
    inline double computeAverage2DThermalSpeed( const int& l_sp, const double* const lp_Twall ) const;
    
protected:
    const Mutation::Thermodynamics::Thermodynamics& m_thermo;

};

/**
 * Gamma Model with recombination probability, gamma ($\gamma$), constant.
 * The rate for species $i$ is given by $\f \gamma M_i^\downarrow $\f,
 * with $M_i^\downarrow$ being the impinging flux at the wall.
 */

class GammaModelConst: public CatalysisRateLaw{
  
public:
    GammaModelConst( const Mutation::Utilities::IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics& thermo, std::vector<int> reactants );
    ~GammaModelConst() {}
  
    GammaModelConst* clone() const {
        return new GammaModelConst( *this );
    }
    
    double forwardReactionRate( const double* const lp_rhoi, const double* const lp_Twall ) const;
    
private:

    mutable int l_index_v_reactants; 
    mutable int l_species_index;
    mutable int l_stoichiometric_coefficient;

    std::vector<double> v_gamma;

    mutable std::vector<double> v_output_impinging_flux;
    mutable std::vector<double> v_impinging_flux_per_stoichiometric_coefficient;

    std::vector<int> v_reactants;

    inline void getSpeciesIndexandStoichiometricCoefficient( int lp_reactants_index, int& l_species_index, int& l_stoichiometric_coefficient ) const;
    inline double computeWallImpingingMassFlux( const int& l_sp, const double* const lp_rhoi, const double* const lp_Twall ) const;
    inline double getLimitingImpingingMassFlux() const;

    
};

/**
 * @brief The Physisorption class
 */
class Physisorption: public CatalysisRateLaw{
  
public:
    Physisorption( const Mutation::Utilities::IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics& thermo );
    ~Physisorption(){ }
  
    Physisorption* clone() const {
        return new Physisorption(*this);
    }
    
    double forwardReactionRate(const double* const p_rhoi, const double* const p_rhoie) const;
    
/**
 * Guerta
 * rate_coefficient = steric_factor_for_phys_sites*(1 - fraction_of_surface_covered_with_chemisorption_sites)
 * thermal_velocity / (4 * total_number_of_physisorption_sites) * exp(-activation_energy_for_phys/(RU*Tw)) // (Schwartzentruber)* S_o (Sticking Coefficient)
 */

private:
    double m_S_coef;
    double m_beta;
    double m_E_act;
  
};

/**
 * @brief The ThermalDesorption class
 */
class ThermalDesorption: public CatalysisRateLaw{
  
public:
    ThermalDesorption(const Mutation::Utilities::IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics& thermo);
    ~ThermalDesorption(){ }
    
    ThermalDesorption* clone() const {
        return new ThermalDesorption(*this);
    }
    
    double forwardReactionRate(const double* const p_rhoi, const double* const p_rhoie) const;
/**
 * Guerta 
 * rate_coefficient = frequency_factor * exp(- activation_energy_for_th_desorp/ (RU * Tw))
 */

private:
    double m_steric_factor;
    double m_vib_perp;
    double m_beta;
    double m_E_des;

};

/**
 * @brief The Chemisorption class
 */
class Chemisorption: public CatalysisRateLaw{

public:
    Chemisorption(const Mutation::Utilities::IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics& thermo);
    ~Chemisorption(){ }
    
    Chemisorption* clone() const {
       return new Chemisorption(*this); 
    }
    
    double forwardReactionRate(const double* const p_rhoi, const double* const p_rhoie) const; /** @todo Why should this be a const? **/
/**
 * Guerta
 * rate_coefficient = steric_factor_for_chem_sites * fraction_of_surface_covered_with_chemisorption_sites *
 * thermal_velocity / (4 * total_number_of_chemisorption_sites) * exp(- activation_energy_for_chem/ (RU * Tw))
 */

private: 
    double m_S_coef;
    double m_beta;
    double m_E_act;
  
};

/**
 * @brief The ERRecombination class
 */
class ERRecombination: public CatalysisRateLaw{
public:
    ERRecombination(const Mutation::Utilities::IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics& thermo);
    ~ERRecombination(){ }
    
    ERRecombination* clone() const {
        return new ERRecombination(*this);
    }
  
    double forwardReactionRate(const double* const p_rhoi, const double* const p_rhoie) const;
    
private:
    double m_steric_factor;
    double m_beta;
    double m_E_act_recomb;
  
};

/**
 * @brief The PhysisorptiontoChemisorption class
 */
class PhysisorptiontoChemisorption: public CatalysisRateLaw{
public:
    PhysisorptiontoChemisorption(const Mutation::Utilities::IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics& thermo);
    ~PhysisorptiontoChemisorption(){ }
    
    PhysisorptiontoChemisorption* clone() const {
        return new PhysisorptiontoChemisorption(*this);
    }

    double forwardReactionRate(const double* const p_rhoi, const double* const p_rhoie) const {
     std::cout << "Debug Point: forwardReactionRate: PhysisorptiontoChemisorption" << std::endl;
     return 1.0;
    };
    
private:
    double m_steric_factor;
    double m_vib_par;
    double m_E_diff;
};

/**
 * @brief The LHRecombination class
 */
class LHRecombination: public CatalysisRateLaw{
public:
    LHRecombination(const Mutation::Utilities::IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics& thermo);
    ~LHRecombination(){ }
    
    LHRecombination* clone() const {
        return new LHRecombination(*this);
    }
  
    double forwardReactionRate(const double* const p_rhoi, const double* const p_rhoie) const {
     std::cout << "Debug Point: forwardReactionRate: L-H Recombination" << std::endl;
     return 1.0; 
    };
    
private:
  
};


    } // namespace gsi
} // namespace Mutation

#endif // CATRATELAW_H
