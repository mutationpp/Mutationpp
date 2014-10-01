
#include "Species.h"
#include "Utilities.h"

using namespace std;
using namespace Mutation::Utilities;

namespace Mutation {
    namespace Thermodynamics {

//==============================================================================

Element::Element(const IO::XmlElement &xml_element)
{
    xml_element.getAttribute("name", m_name, 
        "Element must have a name attribute!");
    
    xml_element.getAttribute("charge", m_charge, 0);
    
    IO::XmlElement::const_iterator child_iter = xml_element.begin();
    IO::XmlElement::const_iterator child_end  = xml_element.end();
    
    for ( ; child_iter != child_end; ++child_iter) {
        if (child_iter->tag() == "mw") {
            std::string str;
            child_iter->getAttribute("units", str, "kg/mol");
            m_atomic_mass = Units(str).convertToBase(
                atof(child_iter->text().c_str()));
        }
    }
}

//==============================================================================

Species::Species(const Species& to_copy)
    : m_name(to_copy.m_name),
      m_ground_state_name(to_copy.m_ground_state_name),
      m_mw(to_copy.m_mw),
      m_charge(to_copy.m_charge),
      m_phase(to_copy.m_phase),
      m_type(to_copy.m_type),
      m_level(to_copy.m_level),
      m_stoichiometry(to_copy.m_stoichiometry)
{ }

//==============================================================================

Species::Species(const Species& to_copy, const size_t level)
    : m_name(to_copy.m_name),
      m_ground_state_name(to_copy.m_name),
      m_mw(to_copy.m_mw),
      m_charge(to_copy.m_charge),
      m_phase(to_copy.m_phase),
      m_type(to_copy.m_type),
      m_level(level),
      m_stoichiometry(to_copy.m_stoichiometry)
{ 
    stringstream ss;
    ss << "(" << m_level << ")";
    m_name += ss.str();
}

//==============================================================================

Species::Species(
    const std::string& name, const PhaseType phase,
    const StoichList& stoichiometry, const std::vector<Element>& elements)
    : m_name(name),
      m_ground_state_name(name),
      m_mw(0.0),
      m_charge(0),
      m_phase(phase),
      m_type(ATOM),
      m_level(0),
      m_stoichiometry(stoichiometry)
{
    StoichList::const_iterator iter = m_stoichiometry.begin();
    std::vector<Element>::const_iterator el;
    for ( ; iter != m_stoichiometry.end(); ++iter) {
        for (el = elements.begin(); el != elements.end(); ++el)
            if (el->name() == iter->first) break;
    
        if (el == elements.end()) {
            std::cout << "Error trying to create species " << m_name
                      << " when element " << iter->first
                      << " does not exist in the database!" << std::endl;
            exit(1);
        }
        
        m_mw += iter->second * el->atomicMass();
        m_charge += iter->second * el->charge();
    }
    
    // Determine the species type
    int atoms = nAtoms();
    m_type = (atoms == 0 ? ELECTRON : (atoms == 1 ? ATOM : MOLECULE));
}

//==============================================================================

void Species::checkStoichiometryNameMatching(
    const string& name, map<string, int>& stoich, 
    const vector<Element>& elements)
{
    // First generate a guess of the stoichiometry based on the name of the 
    // species.  We will assume that element names are no more than 2 characters
    // long and that the species name is in "standard" form which is AaBbCc+
    // where the upper case letters represent element names, lower case letters
    // are the optional stoichiometric numbers for the element, and the '+'
    // represents either a '+' or '-' character which represents minus or plus
    // one electron respectively.  If the name does not fit this form, then we
    // forego the check.
    
    // Step 1: Check to see if the name is in standard form
    bool standard = isalpha(name[0]);
    int i = 1;
    
    if (name.length() > 1) {
        for ( ; i < name.length()-1; ++i)
            standard = standard && isalnum(name[i]);
        standard = 
            standard && (isalnum(name[i]) || name[i] == '+' || name[i] == '-');
    }
        
    if (!standard) 
        return;
    
    // Step 2: Now we know that the name is in standard form, parse it to get
    // the stoichiometry
    int length = name.length();
    map<string, int> guess;
    
    // Handle the special case of the electron
    if (name == "e-") {
        guess["e-"] = 1;
        length = 0;
    } else if (name[length-1] == '+') {
        guess["e-"] = -1;
        length--;
    } else if (name[length-1] == '-') {
        guess["e-"] = 1;
        length--;
    }
    
    // At this point we know that the name must be AaBbCc, use finite-state
    // machine to parse    
    enum {
        first_char,
        second_char,
        number
    } state = first_char;

    string element;
    string digits;
    bool found;
    
    i = 0;
    while (i < length) {
        // Decide action
        switch (state)
        {
            case first_char:
                element = name[i];
                if (i == length - 1) {
                    found = false;
                    for (int j = 0; j < elements.size(); ++j)
                        if ((found = (elements[j].name() == element))) break;
                    if (found)
                        guess[element] = guess[element] + 1;
                    else
                        return;
                } else {
                    state = second_char;
                }
                break;
            
            case second_char:
                // Is the second character an alphabetic letter?
                if (isalpha(name[i])) {
                    element += name[i];
                    // If so check and see if the combined first and second
                    // characters form an element name
                    found = false;
                    for (int j = 0; j < elements.size(); ++j)
                        if ((found = (elements[j].name() == element))) break;
                    if (found) {
                        // If it does, then the next character is either part of 
                        // the stoichiometry coefficient or beginning of a new 
                        // element name
                        if (i == length - 1)
                            guess[element] = guess[element] + 1;
                        else
                            state = number;
                    } else {
                        // Otherwise we assume that the first character is an
                        // element and the coefficient is 1 (if not, then the
                        // name is not in standard form so give up)
                        element = element[0];
                        found = false;
                        for (int j = 0; j < elements.size(); ++j)
                            if ((found = (elements[j].name() == element))) break;
                        if (found) {
                            // Add this element with a stoichiometry coefficient
                            // defaulting to 1
                            guess[element] = guess[element] + 1;
                            element = name[i];
                            // If we are on the last character then we should
                            // also add this element to the list
                            if (i == length - 1) {
                                found = false;
                                for (int j = 0; j < elements.size(); ++j)
                                    if ((found = (elements[j].name() == element)))
                                        break;
                                if (found)
                                    guess[element] = guess[element] + 1;
                                else
                                    return;
                            } else {                         
                                state = second_char;
                            }                        
                        } else {
                            return;
                        }
                    }
                // Otherwise it has to be a number
                } else {
                    digits = name[i];
                    
                    if (i == length - 1) {
                        found = false;
                        for (int j = 0; j < elements.size(); ++j)
                            if ((found = (elements[j].name() == element))) break;
                        if (found)
                            guess[element] = guess[element] + 
                                atoi(digits.c_str());
                        else
                            return;
                    } else {
                        state = number;
                    }
                }
                break;
            
            case number:
                if (isalpha(name[i])) {
                    found = false;
                    for (int j = 0; j < elements.size(); ++j)
                        if ((found = (elements[j].name() == element))) break;
                    
                    if (found)
                        guess[element] = guess[element] + 
                            std::max(atoi(digits.c_str()), 1);
                    else
                        return;
                    
                    element = name[i];
                    if (i == length - 1) {
                        found = false;
                        for (int j = 0; j < elements.size(); ++j)
                            if ((found = (elements[j].name() == element))) break;
                        if (found)
                            guess[element] = guess[element] + 1;
                    } else {                        
                        state = second_char;
                    }
                } else {
                    digits += name[i];
                    
                    if (i == length - 1) {
                        found = false;
                        for (int j = 0; j < elements.size(); ++j)
                            if ((found = (elements[j].name() == element))) break;
                        if (found)
                            guess[element] = guess[element] + 
                                atoi(digits.c_str());
                        else
                            return;                            
                    }
                }
                break;
        }
        
        // Update the position
        i++;
    }
    
    // Now compare the guessed stoichiometry with the given stoichiometry
    bool match ;
    map<string, int>::iterator iter;
    
    if (guess.size() > stoich.size()) {
        for (iter = guess.begin(); iter != guess.end(); ++iter)
            if (!(match = (stoich[iter->first] == iter->second))) break;
    } else {
        for (iter = stoich.begin(); iter != stoich.end(); ++iter)
            if (!(match = (guess[iter->first] == iter->second))) break;
    }
    
    if (!match) {
        cerr << "Warning: species \"" << name << "\" stoichiometry does "
             << "not match species name." << endl;
    }
}

//==============================================================================

void swap(Species& s1, Species& s2)
{
    std::swap(s1.m_name, s2.m_name);
    std::swap(s1.m_ground_state_name, s2.m_ground_state_name);
    std::swap(s1.m_mw, s2.m_mw);
    std::swap(s1.m_charge, s2.m_charge);
    std::swap(s1.m_phase, s2.m_phase);
    std::swap(s1.m_type, s2.m_type);
    std::swap(s1.m_level, s2.m_level);
    std::swap(s1.m_stoichiometry, s2.m_stoichiometry);
}

//==============================================================================

    } // namespace Thermodynamics
} // namespace Mutation


