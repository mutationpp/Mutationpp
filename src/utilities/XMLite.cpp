/**
 * @file XMLite.cpp
 *
 * Implements the XmlDocument and XmlElement classes.
 *
 * @see XmlElement
 * @see XmlDocument
 */

/*
 * Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Mutation++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Mutation++.  If not, see
 * <http://www.gnu.org/licenses/>.
 */


#include "Errors.h"
#include "XMLite.h"
#include "StringUtils.h"

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

    if (!xml_file.is_open())
        throw FileNotFoundError(filename);
    
    int line = 1;
    m_elements.push_back(XmlElement(NULL, this));
    while (m_elements.back().parse(xml_file, line)) {
        m_elements.push_back(XmlElement(NULL, this));
    }
    m_elements.pop_back();
    
    xml_file.close();
}

//==============================================================================
//  XmlElement member definitions
//==============================================================================
void XmlElement::parseError(const std::string& error_msg) const
{
    _parseError(document(), line(), error_msg);
}

//==============================================================================

void XmlElement::_parseError(
    const XmlDocument *const p_document, const long int line, 
    const std::string& error)
{
    if (p_document == NULL)
        throw Error("XML error") << error;
    else
        throw FileParseError(p_document->file(), line) << error;
}

//==============================================================================

template < >
std::string XmlElement::getAttribute(
    const std::string& name, std::string& value) const
{
    std::map<std::string, std::string>::const_iterator iter =
        m_attributes.find(name);
    if (iter != m_attributes.end())
        value = iter->second;
    else
        value = "";
    return value;
}

//==============================================================================

template < >
int XmlElement::getAttribute(const std::string& name, int& value) const
{
    std::map<std::string, std::string>::const_iterator iter =
        m_attributes.find(name);
    if (iter != m_attributes.end())
        value = atoi(iter->second.c_str());
    else
        value = 0;
    return value;
}

//==============================================================================

template < >
float XmlElement::getAttribute(const std::string &name, float &value) const
{
    std::map<std::string, std::string>::const_iterator iter =
        m_attributes.find(name);
    if (iter != m_attributes.end())
        value = (float)atof(iter->second.c_str());
    else
        value = 0.0f;
    return value;
}

//==============================================================================

template < >
double XmlElement::getAttribute(const std::string &name, double &value) const
{
    std::map<std::string, std::string>::const_iterator iter =
        m_attributes.find(name);
    if (iter != m_attributes.end())
        value = atof(iter->second.c_str());
    else
        value = 0.0;
    return value;
}

//==============================================================================

template < >
bool XmlElement::getAttribute(const std::string &name, bool &value) const
{
    std::map<std::string, std::string>::const_iterator iter =
        m_attributes.find(name);
    if (iter != m_attributes.end()) {
        string lower = String::toLowerCase(iter->second);
        value = (lower == "true" || lower == "yes");
    } else
        value = false;
    return value;
}

//==============================================================================

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
                if (m_tag == "!--") {
                    previous = initial;
                    state = comment;
                    m_tag.clear();
                    break;
                }
                if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
                    if (m_tag.empty()) {
                        _parseError(mp_document, line, 
                            string("element names must begin directly after") +
                            string(" the start-tag"));
                    } else {
                        state = attributes;
                   }
                } else if (c == '>') {
                    state = el_value;
                } else {
                    m_tag += c;
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
                    //element.parse(is, line, name, attributes);
                    m_children.push_back(element);
                    m_children.back().parse(is, line, name, attributes);
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
                        //element.parse(is, line, name, el_value);
                        m_children.push_back(element);
                        m_children.back().parse(is, line, name, el_value);
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

    // Add special end element
    m_children.push_back(XmlEndElement());
        
    return true;
}

//==============================================================================

        } // namespace IO
    } // namespace Utilities
} // namespace Mutation
