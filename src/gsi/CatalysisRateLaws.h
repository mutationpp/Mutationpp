#ifndef CATRATELAW_H
#define CATRATELAW_H

#include "Utilities.h"

namespace Mutation{
    namespace gsi{
  
  
// Abstract class for the rate laws for a catalytic reaction
class CatalysisRateLaw{
public:

    virtual ~CatalysisRateLaw() { };
    virtual CatalysisRateLaw* clone() const = 0;
    /**
     * Implementing general functions which are necessary for all the reaction laws
     */
    // Thermal Speed
    // Diffusion Speed
};

/**
 * Gamma Model with constant recombination probability. 
 * @todo Add more details for that kind of models.
 */

class GammaModelConst: public CatalysisRateLaw{
  
public:
    GammaModelConst(const Mutation::Utilities::IO::XmlElement& node);
    ~GammaModelConst() { };
  
    GammaModelConst* clone() const {
        return new GammaModelConst(*this);
    }
    
    double getgammaCoefficient();
    
/**
 * This function computes the reaction rate of the reaction. It has units [@todo] and
 * it is given by the formula for no slip and Maxwellian distribution at the wall.
 */

    double* reactionRate();
    
private:
    double m_gamma;
    //matrices nu, mu
    
};

/**
 * @todo Add This kind of models in the future, where gamma is also a function of temperature
 * 
 *class GammaModelT: public CatalysisRateLaw{
 *  
 *public:
 *    GammaModelT(const Mutation::Utilities::IO::XmlElement& node);
 *    ~GammaModelT();
 *  
 *};
 *
 *class GammaModelTP: public CatalysisRateLaw{
 *  
 *public:
 *    GammaModelTP(const Mutation::Utilities::IO::XmlElement& node);
 *    ~GammaModelTP();
 *  
 *};
*/


//class Physisorption: public CatalysisRateLaw{
//  
///**
// * Guerta
// * rate_coefficient = steric_factor_for_phys_sites*(1 - fraction_of_surface_covered_with_chemisorption_sites)
// * thermal_velocity / (4 * total_number_of_physisorption_sites) * exp(-activation_energy_for_phys/(RU*Tw)) // (Schwartzentruber)* S_o (Sticking Coefficient)
// */
//
///**
// * @todo Overall do
// * @todo Separate the different models (How?)
// * @todo Add the thermal velocity
// * @todo This is for atoms. Maybe dissociative physisorption.
// */
//  
//};
//class ThermalDesorption: public CatalysisRateLaw{
//  
///**
// * Guerta 
// * rate_coefficient = frequency_factor * exp(- activation_energy_for_th_desorp/ (RU * Tw))
// */
//
///**
// * @todo Same as above more or less
// */
//  
//};
//
//class Chemisorption: public CatalysisRateLaw{
//
///**
// * Guerta
// * rate_coefficient = steric_factor_for_chem_sites * fraction_of_surface_covered_with_chemisorption_sites *
// * thermal_velocity / (4 * total_number_of_chemisorption_sites) * exp(- activation_energy_for_chem/ (RU * Tw))
// */
//
///**
// * @todo Same as above more or less
// */
//  
//};
//
//class ERRecombination: public CatalysisRateLaw{
//  
//};
//class PhysisorptiontoChemisorption: public CatalysisRateLaw{
//  
//};
//class LHRecombination: public CatalysisRateLaw{
//  
//};
//

    } // namespace gsi
} // namespace Mutation

#endif // CATRATELAW_HPP