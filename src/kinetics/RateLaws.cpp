#include <cassert>
#include <iostream>

#include "RateLaws.h"
#include "Constants.h"

using namespace std;
using namespace Mutation::Utilities;
using namespace Mutation::Utilities::IO;

namespace Mutation {
    namespace Kinetics {

void Arrhenius::setUnits(const XmlElement& node)
{
    assert( node.tag() == "arrhenius_units" );
    std::string a, e;
    node.getAttribute("A", a);
    node.getAttribute("E", e);
    sm_aunits = Units::split(a);
    sm_eunits = Units::split(e);
}

std::vector<Units> _default_aunits() {
    std::vector<Units> units;
    units.push_back("mol");
    units.push_back("m");
    units.push_back("s");
    units.push_back("K");
    return units;
}

std::vector<Units> _default_eunits() {
    std::vector<Units> units;
    units.push_back("J");
    units.push_back("mol");
    units.push_back("K");
    return units;
}

std::vector<Units> Arrhenius::sm_aunits = std::vector<Units>();
std::vector<Units> Arrhenius::sm_eunits = std::vector<Units>();

Arrhenius::Arrhenius(const XmlElement& node, const int order)
{
    assert( node.tag() == "arrhenius" );
    
    if (sm_aunits.empty())
        sm_aunits = _default_aunits();
    if (sm_eunits.empty())
        sm_eunits = _default_eunits();
    
    // Load the temperature exponent (defaults to 0)
	node.getAttribute("n", m_n, 0.0);
    
    // Load the pre-exponential factor (must be given)
    node.getAttribute("A", m_lnA, 
        "Arrhenius rate law must define coefficient A!");
    
    // Convert to correct units based on the order of the reaction and
    // store the log value
    Units A_units = 
        (((sm_aunits[1]^3) / sm_aunits[0])^(order-1)) /
        (sm_aunits[2] * (sm_aunits[3]^m_n));
    m_lnA = std::log(A_units.convertToBase(m_lnA));
    
    // Load the characteristic temperature
    if (node.hasAttribute("Ea")) {
        node.getAttribute("Ea", m_temp);
        // Convert to J/mol and divide by Ru to get characteristic temp
        m_temp = (sm_eunits[0]/sm_eunits[1]).convertToBase(m_temp) / RU;
    } else if (node.hasAttribute("T")) {
        node.getAttribute("T", m_temp);
        // Convert to K
        m_temp = sm_eunits[2].convertToBase(m_temp);
    } else {
        node.parseError("Arrhenius rate law must define coefficient Ea or T!");
    }
}

    } // namespace Kinetics
} // namespace Mutation

