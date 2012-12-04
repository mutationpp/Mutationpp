#include <cctype>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "Species.h"
#include "XMLite.h"
#include "utilities.h"
#include "ParticleRRHO.h"

using namespace std;

Element::Element(XmlElement &xml_element)
{
    xml_element.getAttribute("name", m_name, 
        "Element must have a name attribute!");
    
    xml_element.getAttribute("charge", m_charge, 0);
    
    XmlElement::Iterator child_iter = xml_element.begin();
    XmlElement::Iterator child_end  = xml_element.end();
    
    for ( ; child_iter != child_end; ++child_iter) {
        if (child_iter->tag() == "mw")
            m_atomic_mass = atof(child_iter->text().c_str());
    }
}

Species::Species(const Species& to_copy)
    : m_name(to_copy.m_name),
      m_nasa7_name(to_copy.m_nasa7_name),
      m_nasa9_name(to_copy.m_nasa9_name),
      mp_rrho_model(
          (to_copy.mp_rrho_model == NULL ? NULL : 
              new ParticleRRHO(*(to_copy.mp_rrho_model)))),
      m_mw(to_copy.m_mw),
      m_charge(to_copy.m_charge),
      m_phase(to_copy.m_phase),
      m_stoichiometry(to_copy.m_stoichiometry)
{ }

Species::Species(
    XmlElement &xml_element, const vector<Element> &elements, 
    set<int> &used_elements)
    : mp_rrho_model(NULL)
{
    // Load attribute information
    xml_element.getAttribute("name", m_name);
    m_nasa9_name = m_name;
    m_nasa7_name = m_name;
    
    string phase;
    xml_element.getAttribute("phase", phase, string("gas"));
    phase = utils::StringUtils::toLowerCase(phase);
    
    if (phase == "gas")
        m_phase = GAS;
    else if (phase == "liquid")
        m_phase = LIQUID;
    else if (phase == "solid")
        m_phase = SOLID;
    else {
        cerr << "Invalid phase description for species \"" << m_name
             << "\", can be \"gas\", \"liquid\", or \"solid\"." << endl;
        exit(1);
    }
    
    // Load information stored in child elements
    XmlElement::Iterator child_iter = xml_element.begin();
    XmlElement::Iterator child_end  = xml_element.end();
    
    for ( ; child_iter != child_end; ++child_iter) {
        if (child_iter->tag() == "stoichiometry")
            loadStoichiometry(child_iter->text(), elements);
        else if (child_iter->tag() == "thermodynamics") {
            std::string thermo_type;
            child_iter->getAttribute("type", thermo_type);
            if (thermo_type == "NASA-7")
                child_iter->getAttribute("db_name", m_nasa7_name, m_name);
            else if (thermo_type == "NASA-9")
                child_iter->getAttribute("db_name", m_nasa9_name, m_name);
            else if (thermo_type == "RRHO")
                mp_rrho_model = new ParticleRRHO(*child_iter);
        }
    }
    
    /*cout << m_name << endl;
    ParticleType ptype = type();
    cout << (ptype == ELECTRON ? "ELECTRON" : 
        (ptype == ATOM ? "ATOM" : "MOLECULE")) << endl;
    switch (ptype) {
        
        case MOLECULE:
            cout << "theta_r: " << mp_rrho_model->rotationalTemperature() << endl;
            cout << "linearity: " << mp_rrho_model->linearity() << endl;
            cout << "omega: " << mp_rrho_model->stericFactor() << endl;
            cout << "nvib: " << mp_rrho_model->nVibrationalLevels() << endl;
            for (int i = 0; i < mp_rrho_model->nVibrationalLevels(); ++i)
                cout << mp_rrho_model->vibrationalEnergy(i) << endl;
        case ATOM:
            cout << "nelec: " << mp_rrho_model->nElectronicLevels() << endl;
            for (int i = 0; i < mp_rrho_model->nElectronicLevels(); ++i)
                cout << mp_rrho_model->electronicEnergy(i).first << " "
                     << mp_rrho_model->electronicEnergy(i).second << endl;
        case ELECTRON:
            cout << "hf: " << mp_rrho_model->formationEnthalpy() << endl;
    }*/
    
    // Now use the elemental composition and elements list to determine the
    // molecular weight and charge of this species (ensures mass and charge
    // conservation)
    vector<Element>::const_iterator iter = elements.begin();
    set<string> element_names;
    
    m_mw = 0.0;
    m_charge = 0;
    
    for (int i = 0; iter < elements.end(); ++iter, ++i) {
        if (m_stoichiometry.count(iter->name()) > 0) {
            // Number of atoms of this element belonging to this species
            int atoms = nAtoms(iter->name());
            // Update molecular weight
            m_mw += atoms * iter->atomicMass();
            // Update species charge
            m_charge += atoms * iter->charge();
            // Let the user know that this element is in this species
            used_elements.insert(i);
        }
    }
}

Species::~Species()
{
    if (mp_rrho_model != NULL)
        delete mp_rrho_model;
}

void Species::loadStoichiometry(
    const string& stoichiometry, const vector<Element> &elements)
{
    // Split up the string in the format "A:a, B:b, ..." into a vector of
    // strings matching ["A", "a", "B", "b", ... ]
    vector<string> stoichiometry_tokens;
    utils::StringUtils::tokenize(
        utils::StringUtils::removeWhiteSpace(stoichiometry), 
        stoichiometry_tokens, ":,");
    
    // Check that the vector has an even number of tokens (otherwise there must
    // be an error in the syntax)
    if (stoichiometry_tokens.size() % 2 != 0) {
        cerr << "Error in species \"" << name() << "\" stoichiometry "
             << "definition, invalid syntax!" << endl;
        for (int i = 0; i < stoichiometry_tokens.size(); ++i)
            cerr << "\"" << stoichiometry_tokens[i] << "\"" << endl;
        exit(1); 
    }
    
    // Determine stoichiometry from the stoichiometry vector
    string element_name;
    int atoms;
    bool found;
    
    for (int i = 0; i < stoichiometry_tokens.size(); i+=2) {
        element_name = stoichiometry_tokens[i];
        atoms = atoi(stoichiometry_tokens[i+1].c_str());
    
        // Check that each element name in the stoichiometry list matches an
        // element in the element list that is loaded
        found = false;
        for (int j = 0; j < elements.size(); ++j)
            if (found = (elements[j].name() == element_name)) break;
        if (!found) {
            cerr << "Error in species \"" << name() << "\" stoichiometry "
                 << "definition, element \"" << element_name
                 << "\" is not defined!" << endl;
            exit(1);
        }
        
        // Now add the element's stoichiometry to the map
        m_stoichiometry[stoichiometry_tokens[i]] = 
            atoi(stoichiometry_tokens[i+1].c_str());
    }
    
    // Finally, do a check to determine whether or not the species name matches
    // the given stoichiometry if the name is in "standard" form
    checkStoichiometryNameMatching(m_name, m_stoichiometry, elements);
}

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
    
    int atoms;
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
                        if (found = (elements[j].name() == element)) break;
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
                        if (found = (elements[j].name() == element)) break;
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
                            if (found = (elements[j].name() == element)) break;
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
                                    if (found = (elements[j].name() == element)) 
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
                            if (found = (elements[j].name() == element)) break;
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
                        if (found = (elements[j].name() == element)) break;
                    if (found)
                        guess[element] = guess[element] + 
                            atoi(digits.c_str());
                    else
                        return;
                    
                    element = name[i];
                    if (i == length - 1) {
                        found = false;
                        for (int j = 0; j < elements.size(); ++j)
                            if (found = (elements[j].name() == element)) break;
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
                            if (found = (elements[j].name() == element)) break;
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


void swap(Species& s1, Species& s2)
{
    std::swap(s1.m_name, s2.m_name);
    std::swap(s1.m_nasa7_name, s2.m_nasa7_name);
    std::swap(s1.m_nasa9_name, s2.m_nasa9_name);
    std::swap(s1.mp_rrho_model, s2.mp_rrho_model);
    std::swap(s1.m_mw, s2.m_mw);
    std::swap(s1.m_charge, s2.m_charge);
    std::swap(s1.m_phase, s2.m_phase);
    std::swap(s1.m_stoichiometry, s2.m_stoichiometry);
}


