
#include "MillikanWhite.h"
#include "Utilities.h"
#include "ParticleRRHO.h"

using namespace Mutation::Utilities;
using namespace Mutation::Utilities::IO;
using namespace Mutation::Thermodynamics;

namespace Mutation {
    namespace Transfer {

//==============================================================================

MillikanWhiteVibrator::MillikanWhiteVibrator(
    const XmlElement& node, const class Thermodynamics& thermo)
{
    assert(node.tag() == "vibrator");
        
    // Get the name of this species
    std::string name;
    node.getAttribute("species", name, "must provide species name!");
    m_index = thermo.speciesIndex(name);
    
    // Get the limiting cross-section if available
    node.getAttribute("omegav", m_omegav, 3.0E-21);
    
    const Species& vibrator = thermo.species(name);
    
    // Get characteristic vibrational temperature
    
    double theta=0.0;    //Probably the characteristic vibrational energy is not needed.
    /*//if (vibrator.hasRRHOParameters())
    //    theta = vibrator.getRRHOParameters()->vibrationalEnergy(0);
    //else {
        cout << "Cannot get characteristic vibrational temperature for "
             << "species " << vibrator.name() << " because there is no "
             << "RRHO data present in species.xml!" << endl;
    //    exit(1);
    //}*/
    
    // Loop over each heavy species in thermo
    int offset = (thermo.hasElectrons() ? 1 : 0);
    XmlElement::const_iterator partner_iter;
    
    double a, b, mu;
    for (int i = 0; i < thermo.nHeavy(); ++i) {
        // Get collision partner
        const Species& partner = thermo.species(i+offset);
        
        // Compute reduced mass of this pair
        mu = (vibrator.molecularWeight() * partner.molecularWeight()) /
             (vibrator.molecularWeight() + partner.molecularWeight()) * 1.0E3;
            
        // Use a and b data from data file or use the defaults if the pair
        // is not present in the file
        if ((partner_iter = node.findTagWithAttribute(
            "partner", "species", partner.name())) != node.end()) {
            // Get a, b from parnter node
            partner_iter->getAttribute("a", a, "must provide constant a!");
            partner_iter->getAttribute("b", b, "must provide constant b!");
            
            // add Millikan-White data for collision pair
            m_partners.push_back(MillikanWhitePartner(a, b, mu));
        } else {
            // Add Millikan-White data using defaults
            m_partners.push_back(MillikanWhitePartner(mu, theta));
        }
    }
}

//==============================================================================

MillikanWhiteVibrator::MillikanWhiteVibrator(
    const std::string& name, const class Thermodynamics& thermo)
    : m_omegav(3.0E-21), m_index(thermo.speciesIndex(name))
{
    const Species& vibrator = thermo.species(name);
    
    // Get characteristic vibrational temperature
    double theta, mu;
   // if (vibrator.hasRRHOParameters())
   //     theta = vibrator.getRRHOParameters()->vibrationalEnergy(0);
   // else {
        cout << "Cannot get characteristic vibrational temperature for "
             << "species " << vibrator.name() << " because there is no "
             << "RRHO data present in species.xml!" << endl;
//        exit(1);
    //}
    
    // Loop over each heavy species in thermo
    int offset = (thermo.hasElectrons() ? 1 : 0);
    
    for (int i = 0; i < thermo.nHeavy(); ++i) {
        // Get collision partner
        const Species& partner = thermo.species(i+offset);
        
        // Compute reduced mass of this pair
        mu = (vibrator.molecularWeight() * partner.molecularWeight()) /
             (vibrator.molecularWeight() + partner.molecularWeight()) * 1.0E3;
        
        // Add Millikan-White data using defaults
        m_partners.push_back(MillikanWhitePartner(mu, theta));
    }
}

//==============================================================================

MillikanWhite::MillikanWhite(const class Thermodynamics& thermo)
{
    std::string filename =
        getEnvironmentVariable("MPP_DATA_DIRECTORY") + "/transfer/VT.xml";
    
    // Open the VT file as an XML document
    XmlDocument doc(filename);
    XmlElement root = doc.root();
    
    // Search for the Millikan and White element
    XmlElement::const_iterator iter = root.begin();
    for ( ; iter != root.end(); ++iter) {
        if (iter->tag() == "Millikan-White") {
            root = *iter;
            break;
        }
    }

    // Loop over all molecules and load the Millikan-White data associated with
    // them
    int offset = (thermo.hasElectrons() ? 1 : 0);
    for (int i = 0; i < thermo.nHeavy(); ++i) {
        const Species& species = thermo.species(i+offset);
        
        // If this molecule can vibrate, add it to the list
        if (species.type() == MOLECULE) {
            if ((iter = root.findTagWithAttribute(
                "vibrator", "species", species.name())) != root.end())
                m_vibrators.push_back(
                    MillikanWhiteVibrator(*iter, thermo));
            else
                m_vibrators.push_back(
                    MillikanWhiteVibrator(species.name(), thermo));
        }
    }
}

//==============================================================================

    } // namespace Transfer
} // namespace Mutation