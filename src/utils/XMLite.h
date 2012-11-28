/**
 * @file Very simple XML library which enables loading XML data from a file into
 * a simple data structure which can be mined.  The library consists of
 * minimal error handling but parses simple XML structures just fine.
 *
 * @author J.B. Scoggins (jbscoggi@gmail.com)
 * @date   January 24, 2012
 */
#ifndef XMLITE_H
#define XMLITE_H

#include "IteratorWrapper.h"

#include <istream>
#include <map>
#include <string>
#include <sstream>
#include <vector>

// Forward declaration of XmlDocument
class XmlDocument;

/**
 * Represents an XML element node which can contain (name,value) pair attributes
 * as well as either element children or a value.
 *
 * @todo Should make the mp_document and mp_parent members be constatnt pointers
 * and then implement the = operator and clone stuff (forget the name of this
 * construct).
 */
class XmlElement
{
private:

    enum ParseState {
        initial,
        prename,
        el_name,
        el_value,
        child_name,
        attributes,
        att_name,
        att_value,
        comment,
        comment_end_1,
        comment_end_2,
        done
    };
    
public:

    XmlElement(XmlElement *p_parent, XmlDocument *p_document)
        : mp_parent(p_parent), mp_document(p_document), m_line_number(0)
    { }
    
    XmlElement* parent() const {
        return const_cast<XmlElement *const>(mp_parent);
    }
    
    XmlDocument* document() const {
        return const_cast<XmlDocument *const>(mp_document);
    }
    
    const std::string &tag() const {
        return m_tag;
    }
    
    const std::string &text() const {
        return m_text;
    }
    
    long int line() const {
        return m_line_number;
    }
    
    bool hasAttribute(const std::string &name) const {
        return (m_attributes.count(name) > 0);
    }
    
    template <typename T>
    void getAttribute(const std::string& name, T& value, const T& dflt)
    {
        if (hasAttribute(name))
            getAttribute(name, value);
        else
            value = dflt;
    }
	
	template <typename T>
	void getAttribute(const std::string& name, T& value, 
        const char *const error_msg)
	{
        if (hasAttribute(name))
            getAttribute(name, value);
        else
            parseError(error_msg);
	}
        
    template <typename T>
    void getAttribute(const std::string& name, T& value);
    
    /**
     * Children iterator class.
     */
    class Iterator 
        : public IteratorWrapper<XmlElement, std::vector<XmlElement>::iterator>
    {
        friend class XmlElement;
    };
    
    Iterator &begin()
    {
        m_iterator.m_iter = m_children.begin();
        return m_iterator;
    }
    
    Iterator &end()
    {
        m_iterator.m_iter = m_children.end();
        return m_iterator;
    }
    
    void parseError(const char *const error_msg) {
        parseError(std::string(error_msg));
    }
    
    void parseError(const std::string& error_msg);
    
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
    Iterator                           m_iterator;
    long int                           m_line_number;
    
};

/**
 * Represents an XML document which stores an array of each top-level element
 * contained in a file.
 */
class XmlDocument
{
public:    

    XmlDocument(const std::string &filename);
    
    XmlElement &root() {
        return m_elements[0];
    }
    
    const std::string& file() const {
        return m_filename;
    }

    class Iterator 
        : public IteratorWrapper<XmlElement, std::vector<XmlElement>::iterator>
    {
        friend class XmlDocument;
    };
    
    Iterator &begin()
    {
        m_iterator.m_iter = m_elements.begin();
        return m_iterator;
    }
    
    Iterator &end()
    {
        m_iterator.m_iter = m_elements.end();
        return m_iterator;
    }
    
    friend void _parseError(const int, const std::string &);

private:

    const std::string       m_filename;
    std::vector<XmlElement> m_elements;
    Iterator                m_iterator;
};

#endif // XMLITE_H
