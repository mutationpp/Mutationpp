#ifndef CATRATELAW_H
#define CATRATELAW_H

#include "Utilities.h"

namespace Mutation{
    namespace gsi{
  
  
// Abstract class for the rate laws for a catalytic reaction
class CatalysisRateLaw{
    //CatalysisRateLaw(const Mutation::Utilities::IO::XmlElement& node);
};

/**
 * Gamma Model with constant recombination probability. 
 * @todo Add more details for that kind of models.
 */

class GammaModelConst: public CatalysisRateLaw{
  
public:
    GammaModelConst(const Mutation::Utilities::IO::XmlElement& node);
    ~GammaModelConst();
  
private:
    double m_gamma;
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

/*
class Physisorption: public CatalysisRateLaw{
  
};
class ThermalDesorption: public CatalysisRateLaw{
  
};
class Chemisorption: public CatalysisRateLaw{
  
};
class ERRecombination: public CatalysisRateLaw{
  
};
class PhysisorptiontoChemisorption: public CatalysisRateLaw{
  
};
class LHRecombination: public CatalysisRateLaw{
  
};
*/

    } // namespace gsi
} // namespace Mutation

#endif // CATRATELAW_HPP