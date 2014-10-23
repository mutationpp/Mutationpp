#include <cassert>

#include "CatalysisRateLaws.h"

namespace Mutation{
    namespace gsi{
      
//==============================================================================

GammaModelConst::GammaModelConst(const Mutation::Utilities::IO::XmlElement& node){
    assert( node.tag() == "gamma_const" );
    
    node.getAttribute("gamma", m_gamma, "In the gamma model the gamma parameter should be defined");
    
}
      
//==============================================================================
/*
GammaModelConst::GammaModelT(const Mutation::Utilities::IO::XmlElement& node){
    assert( node.tag() == "gamma_T" );
}

//==============================================================================

GammaModelConst::GammaModelTP(const Mutation::Utilities::IO::XmlElement& node){
    assert( node.tag() == "gamma_const" );
    
}
    */
//==============================================================================
      
    } // namespace gsi
} // namespace Mutation