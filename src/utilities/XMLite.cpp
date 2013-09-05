/**
 * @file XMLite.cpp
 *
 * Implements the XmlDocument and XmlElement classes.
 *
 * @see XmlElement
 * @see XmlDocument
 * @see XMLite.h
 *
 * @author J.B. Scoggins (jbscoggi@gmail.com)
 * @date   January 24, 2012
 */

#include "XMLite.h"

#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

namespace Mutation {
    namespace Utilities {
        namespace IO {

//==============================================================================
//  XmlDocument member definitions
//==============================================================================
XmlDocument::XmlDocument(const std::string &filename)
    : m_filename(filename)
{
    ifstream xml_file(m_filename.c_str(), ios::in);
    
    if (!xml_file.is_open()) {
        cerr << "Could not open file " << m_filename << " for reading!" << endl;
        exit(1);
    }
    
    int line = 1;
    XmlElement element(NULL, this);
    while (element.parse(xml_file, line)) {
        m_elements.push_back(element);
        element = XmlElement(NULL, this);
    }
    
    xml_file.close();
}

//==============================================================================
//  XmlElement member definitions
//==============================================================================
void XmlElement::parseError(const std::string& error_msg) const
{
    _parseError(document(), line(), error_msg);
}

void XmlElement::_parseError(
    const XmlDocument *const p_document, const long int line, 
    const std::string& error)
{
    cerr << "XML error: "<< p_document->file() << ": line " << line 
         << ": " << error << endl;
    exit(1);
}

template < >
void XmlElement::getAttribute(
    const std::string &name, std::string &value)
{
    value = m_attributes[name];
}

template < >
void XmlElement::getAttribute(const std::string &name, int &value)
{
    value = atoi(m_attributes[name].c_str());
}

template < >
void XmlElement::getAttribute(const std::string &name, float &value)
{
    value = (float)atof(m_attributes[name].c_str());
}

template < >
void XmlElement::getAttribute(const std::string &name, double &value)
{
    value = atof(m_attributes[name].c_str());
}

template < >
void XmlElement::getAttribute(const std::string &name, bool &value)
{
    value = (m_attributes[name] == "true");
}

bool XmlElement::parse(
    istream &is, int &line, string name, ParseState state)
{
    // State-machine logic for reading in XML data
    char c;
    m_tag = name;
    name = "";
    string value = "";
    
    // Record the start line for error statements
    m_line_number = line;
    ParseState previous = state;
    bool done = false;
  
    while (!done && is.get(c)) {
        if (c == '\n')
            line++;
        
        switch (state) {
            case initial:
                if ( !(c == ' ' || c == '\t' || c == '\n' || c == '\r') ) {
                    m_line_number = line;
                    if (c == '<')
                        state = el_name;
                    else
                        _parseError(mp_document, line, 
                            string("encountered character other than ") + 
                            string("start-tag '<' at the top level of the") +
                            string(" document"));
                }
                break;
            case el_name:
                if (name == "!--") {
                    previous = initial;
                    state = comment;
                    name.clear();
                    break;
                }
                if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
                    if (name.empty()) {
                        _parseError(mp_document, line, 
                            string("element names must begin directly after") +
                            string(" the start-tag"));
                    } else {
                        state = attributes;
                        m_tag = name;
                        name.clear();
                   }
                } else if (c == '>') {
                    state = el_value;
                    m_tag = name;
                    name.clear();
                } else {
                    name += c;
                }
                break;
            case attributes:
                if (!(c == ' ' || c == '\t' || c == '=' || c == '\n' || c == '\r')) {
                    if (c == '/') {
                        done = true;
                        is.get();
                    } else if (c == '>') {
                        state = el_value;
                    } else if (c == '\"') {
                        state = att_value;
                    } else {
                        state = att_name;
                        is.unget();
                    }
                }
                break;
            case att_name:
                if (!(c == '=' || c == ' ' || c == '\t'))
                    name += c;
                else
                    state = attributes;
                break;
            case att_value:
                if (c != '\"') {
                    value += c;
                } else {
                    state = attributes;
                    m_attributes[name] = value;
                    name.clear();
                    value.clear();
                }
                break;
            case el_value:
                if (c == '<') {
                    state = child_name;
                    m_text = value;
                    value.clear();
                } else {
                    value += c;
                }
                break;
            case child_name:
                if (name == "!--") {
                    previous = el_value;
                    state = comment;
                    name.clear();
                    break;
                }
                if (c == ' ') {
                    XmlElement element(this, mp_document);
                    element.parse(is, line, name, attributes);
                    m_children.push_back(element);
                    name = "";
                    state = el_value;
                } else if (c == '>') {
                    if (name[0] == '/') {
                        if (name == "/"+m_tag)
                            done = true;
                        else
                            _parseError(mp_document, line,
                                "expecting end-tag </"+m_tag+"> but instead "+
                                "found <"+name+">");
                    } else {
                        XmlElement element(this, mp_document);
                        element.parse(is, line, name, el_value);
                        m_children.push_back(element);
                        name = "";
                        state = el_value;
                    }
                } else {
                    name += c;
                }
                break;
            case comment:
                if (c == '-')
                    state = comment_end_1;
                break;
            case comment_end_1:
                if (c == '-')
                    state = comment_end_2;
                else
                    state = comment;
                break;
            case comment_end_2:
                if (c == '>')
                    state = previous;
                else
                    _parseError(mp_document, line, 
                        "cannot use \"--\" in comment except for end-tag");
                break;
        } // state switch
    } // while not done with state-machine
    
    // Make sure we ended where we thought we would
    // note we could be more descriptive with the error statement by seeing
    // what state actually is
    if (!done) {
        if (m_tag.empty())
            return false;
        else
            _parseError(mp_document, m_line_number, 
                "reached end of file before element " + m_tag + " ended");
    }
    
    // If this element has children then force no value
    if (m_children.size() > 0)
        m_text.clear();
        
    return true;
}

        } // namespace IO
    } // namespace Utilities
} // namespace Mutation
