/**
 * @file XMLite.h
 *
 * @brief Simple XML library which enables loading XML data from a file into
 * a simple data structure which can be mined.  The library consists of
 * minimal error handling but parses simple XML structures just fine.
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

#ifndef XMLITE_H
#define XMLITE_H

#include <istream>
#include <map>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>

namespace Mutation {
    namespace Utilities {
        namespace IO {

// Forward declaration of XmlDocument
class XmlDocument;

/**
 * Represents an XML element node which can contain (name,value) pair attributes
 * as well as either element children or a text value.
 *
 * @todo Should make the mp_document and mp_parent members be constant pointers
 * and then implement the = operator and clone stuff (forget the name of this
 * construct).
 */
class XmlElement
{
private:

    enum ParseState {
        initial,
        el_name,
        el_value,
        child_name,
        attributes,
        att_name,
        att_value,
        comment,
        comment_end_1,
        comment_end_2,
    };
    
public:

    typedef std::vector<XmlElement>::iterator iterator;
    typedef std::vector<XmlElement>::const_iterator const_iterator;

    /**
     * Constructs a new XmlElement given its parent and document pointers.
     */
    XmlElement(XmlElement *p_parent, XmlDocument *p_document)
        : mp_parent(p_parent), mp_document(p_document), m_line_number(0)
    { }
    
    /**
     * Constructs an element from an input string.
     */
    XmlElement(const std::string& str)
        : mp_parent(NULL), mp_document(NULL), m_line_number(0)
    {
        int line = 0;
        std::istringstream is(str);
        parse(is, line);
    }

    XmlElement() :
        mp_parent(NULL), 
        mp_document(NULL),
        m_attributes(),
        m_children(),
        m_tag(),
        m_text(),
        m_line_number(0)
    { }

    /**
     * Returns a pointer to the parent XmlElement object.
     */
    XmlElement* parent() const {
        return const_cast<XmlElement *const>(mp_parent);
    }
    
    /**
     * Returns a pointer to the XmlDocument object to which this element 
     * belongs.
     */
    XmlDocument* document() const {
        return const_cast<XmlDocument *const>(mp_document);
    }
    
    /**
     * Returns the tag associated with this element.
     */
    const std::string &tag() const {
        return m_tag;
    }
    
    /**
     * Returns the text associated with this element.
     */
    const std::string &text() const {
        return m_text;
    }
    
    /** 
     * Returns the line number in the file where this element begins.
     */
    long int line() const {
        return m_line_number;
    }
    
    /**
     * Returns true if the given attribute name exists in the attributes
     * belonging to this element, false otherwise.
     */
    bool hasAttribute(const std::string &name) const {
        return (m_attributes.count(name) > 0);
    }
    
    /**
     * Fills the value parameter with the value of the corresponding attribute
     * with the given name.  If the attribute is not present in the element, the
     * given default value is used.
     */
    template <typename T>
    void getAttribute(const std::string& name, T& value, const T& dflt) const
    {
        if (hasAttribute(name))
            getAttribute(name, value);
        else
            value = dflt;
    }
	
    /**
     * Fills the value parameter with the value of the corresponding attribute
     * with the given name.  If the attribute is not present in the element,
     * the error message is printed and execution is terminated.
     */
	template <typename T>
	void getAttribute(
        const std::string& name, T& value, const char *const error_msg) const
	{
        if (hasAttribute(name))
            getAttribute(name, value);
        else
            parseError(error_msg);
	}
    
    /**
     * Fills the value parameter with the value of the corresponding attribute
     * with the given name.
     */
    template <typename T>
    T getAttribute(const std::string& name, T& value = T()) const;
    
    /**
     * Returns an iterator pointing to the first child element in this
     * XmlElement.
     */
    const_iterator begin() const {
        return m_children.begin();
    }
    
    /**
     * Returns an iterator pointing to the end of the child element list.
     */
    const_iterator end() const {
        return (m_children.end()-1);
    }
    
    /**
     * Returns an iterator pointing to the first child element which has the
     * given tag.
     */
    const_iterator findTag(const std::string& tag) const
    {
        return findTag(tag, begin());
    }
    
    /**
     * Returns an iterator pointing to the first child element which has the
     * given tag starting from the given iterator.
     */
    const_iterator findTag(const std::string& tag, const_iterator iter) const
    {
        while (iter != end()) {
            if (iter->tag() == tag)
                break;
            ++iter;
        }
        
        return iter;
    }
    
    /**
     * Returns an iterator pointing to the first child element which has the 
     * given tag and given attribute/value pair.  If no such child element 
     * exists, the iterator equals end().
     */
    template <typename T>
    const_iterator findTagWithAttribute(
        const std::string& tag, const std::string& attribute, const T value) const
    {
        return findTagWithAttribute(tag, attribute, value, begin());
    }
    
    /**
     * Returns an iterator pointing to the first child element which has the 
     * given tag and given attribute/value pair.  If no such child element 
     * exists, the iterator equals end().
     */
    const_iterator findTagWithAttribute(
        const std::string& tag, const std::string& attribute, const char* value) const
    {
        return findTagWithAttribute(tag, attribute, std::string(value), begin());
    }
    
    /**
     * Returns an iterator pointing to the first child element (starting from
     * iter) which has the given tag and given attribute/value pair.  If no such
     * child element exists, the iterator equals end().
     */
    template <typename T>
    const_iterator findTagWithAttribute(
        const std::string& tag, const std::string& attribute, const T value,
        const_iterator iter) const
    {
        iter = findTag(tag, iter);
        while (iter != end()) {
            T val;
            iter->getAttribute(attribute, val);
                
            if (val == value)
                break;
            
            iter = findTag(tag, ++iter);
        }
        
        return iter;
    }
    
    
    /**
     * Returns an iterator pointing to the first child element (starting from
     * iter) which has the given tag and given attribute/value pair.  If no such
     * child element exists, the iterator equals end().
     */
    const_iterator findTagWithAttribute(
        const std::string& tag, const std::string& attribute, const char* value,
        const_iterator iter) const
    {
        return findTagWithAttribute(tag, attribute, std::string(value), iter);
    }
    
    template <typename T>
    bool getChildElementObject(const std::string& tag, T& object) const
    {
        getChildElementObject(tag, object, begin());
    }

    template <typename T>
    bool getChildElementObject(const std::string& tag, T& object,
         const_iterator iter) const
    {
        iter = findTag(tag, iter);
        if (iter != end()) {
            object = T(*iter);
            return true;
        }

        return false;
    }

    /** 
     * Prints the given error message and exits.
     */
    void parseError(const char *const error_msg) const {
        parseError(std::string(error_msg));
    }
    
    /** 
     * Prints the given error message and exits.
     */
    void parseError(const std::string& error_msg) const;
    
    /**
     * Same as
     * \code
     * if (!test) parseError(error_msg);
     * \endcode
     */
    void parseCheck(bool test, const char *const error_msg) const {
        if (!test) parseError(error_msg);
    }

    /**
     * Same as
     * \code
     * if (!test) parseError(error_msg);
     * \endcode
     */
    void parseCheck(bool test, const std::string& error_msg) const {
        if (!test) parseError(error_msg);
    }

    friend class XmlDocument;

private:    
    
    static void _parseError(
        const XmlDocument *const p_document, const long int line, 
        const std::string& error);
    
private:

    bool parse(
        std::istream &is, int &line, std::string name = "",
        ParseState state = initial);

    XmlElement*                        mp_parent;
    XmlDocument*                       mp_document;
    std::map<std::string, std::string> m_attributes;
    std::vector<XmlElement>            m_children;
    std::string                        m_tag;
    std::string                        m_text;
    long int                           m_line_number;
    
};

/**
 * @brief Special XmlElement which represents end of children list.
 * 
 * This class is necessary in order to avoid rare bug in which the end() pointer
 * for XmlElement accidently points to another real element (due to how the
 * memory is allcoated).
 */
class XmlEndElement : public XmlElement
{
public:
    XmlEndElement()
        : XmlElement()
    { }
};

/**
 * Represents an XML document which stores an array of each top-level element
 * contained in a file.
 */
class XmlDocument
{
public:

    typedef std::vector<XmlElement>::iterator iterator;
    typedef std::vector<XmlElement>::const_iterator const_iterator;

    XmlDocument(const std::string &filename);
    
    XmlElement &root() {
        return m_elements[0];
    }
    
    const std::string& file() const {
        return m_filename;
    }
    
    const_iterator begin() const {
        return m_elements.begin();
    }
    
    const_iterator end() const {
        return m_elements.end();
    }
    
    friend void _parseError(const int, const std::string &);

private:

    const std::string       m_filename;
    std::vector<XmlElement> m_elements;
};

        } // namespace IO
    } // namespace Utilities
} // namespace Mutation

#endif // XMLITE_H

