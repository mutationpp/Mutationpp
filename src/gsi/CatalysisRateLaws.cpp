#include <cassert>

#include "CatalysisRateLaws.h"
#include "Constants.h"

namespace Mutation{
    namespace gsi{

//==============================================================================
//==============================================================================
//==============================================================================
      
GammaModelConst::GammaModelConst(const Mutation::Utilities::IO::XmlElement& node)
{
    assert( node.tag() == "gamma_const" );
    
    node.getAttribute("gamma", m_gamma, "In the gamma model the gamma parameter should be defined");
    
    std::cout << "Debug Point: CatalysisRateLaws gamma = " << m_gamma << std::endl;
    
    /**
     * @todo Add Controls here? Or in the rate manager
     * @todo What should I add here?
     */
    
}

//==============================================================================

//double GammaModelConst::getgammaCoefficient()
//{
//    return m_gamma;
//}

//==============================================================================


//double* GammaModelConst::reactionRate() 
//{
//  /**
//   * @todo This is for the non-general case of no slip and Maxwell. Add the other flux 
//   * functions and make them available in an interface
//   */
//  
//    double output = m_gamma * sqrt(KB / (2 * PI));
//}

//==============================================================================

/*
GammaModelT::GammaModelT(const Mutation::Utilities::IO::XmlElement& node){
    assert( node.tag() == "gamma_T" );
}

//==============================================================================

GammaModelTP::GammaModelTP(const Mutation::Utilities::IO::XmlElement& node){
    assert( node.tag() == "gamma_const" );
    
}
    */
//==============================================================================
      
    } // namespace gsi
} // namespace Mutation