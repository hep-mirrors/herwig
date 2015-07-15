// -*- C++ -*-
//
// Element.hpp is a part of myXML
// Copyright (C) 2012-2013 Simon Platzer
//
// myXML is licenced under version 2 of the GPL, see COPYING for details.
//

#ifndef MYXML_Element_hpp_included
#define MYXML_Element_hpp_included

#include <string>
#include <map>
#include <list>
#include <sstream>
#include <iomanip>

namespace XML {

  /**
   * \brief ElementTypes contains the element type enumeration
   * \author Simon Platzer
   */
  namespace ElementTypes {

    /**
     * The element type enumeration
     */
    enum EnumerateElementTypes {

      Unknown = -1,
      /** Unknown element type */
      Root = 0,
      /** The entire document */
      EmptyElement = 1,
      /** An empty element */
      Element = 2,
      /** An element */
      ProcessingInstruction = 3,
      /** A processing instruction */
      CharacterData = 4,
      /** Character data */
      ParsedCharacterData = 5,
      /** Parsed character data */
      Comment = 6
      /** A comment */

    };

  }

  /**
   * \brief Element represents a (tree of) XML elements
   * \author Simon Platzer
   */
  class Element {

  public:

    /**
     * The default constructor
     */
    Element()
      : theType(ElementTypes::Unknown) {}

    /**
     * The standard constructor
     */
    explicit Element(int newType,
		     const std::string& newNameOrContent = "")
      : theType(newType), theNameOrContent(newNameOrContent) {}


    /**
     * The copy constructor
     */
    Element(const Element& other);

    /**
     * Assignment
     */
    Element& operator=(const Element& other);
  
    /** 
     * Comparison operator
     *
     */
    friend bool operator==(const XML::Element &one, const XML::Element &two);
    friend bool operator!=(const XML::Element &one, const XML::Element &two);

    /**
     * Combine operator
     * 
     * Operator checks if type and name is equal otherwises throws exception
     */
    friend XML::Element operator+(const XML::Element& one, const XML::Element& two);
    
  public:

    /**
     * Return the type of this element
     */
    int type() const { return theType; }

    /**
     * Return the name of this element
     */
    const std::string& name() const;

  public:

    /**
     * Return the content of this element
     */
    const std::string& content() const;

    /**
     * Access the content of this element
     */
    std::string& content();

    /**
     * Append to the content of this element
     */
    template<class T>
    Element& operator<<(const T& t) {
      assertContent();
      std::ostringstream out;
      out << std::setprecision(16) << t;
      theNameOrContent += out.str();
      return *this;
    }

  public:

    /**
     * Return true, if this element has attributes
     */
    bool hasAttributes() const {
      return 
	type() == ElementTypes::Element ||
	type() == ElementTypes::EmptyElement;
    }

    /**
     * Return true, if this element contains an attribute of the given name
     */
    bool hasAttribute(const std::string&) const;

    /**
     * Return the attributes
     */
    const std::map<std::string,std::string>& attributes() const;

    /**
     * Return the attribute of the given name
     */
    const std::string& attribute(const std::string&) const;

    /**
     * Access the attribute of the given name
     */
    std::string& attribute(const std::string&);

    /**
     * Represent an attribute
     */
    struct Attribute {

      /**
       * The attribute name
       */
      std::string name;

      /**
       * The attribute value
       */
      std::string value;

      /**
       * Construct an attribute
       */
      template<class T>
      Attribute(const std::string& newName,
		const T& newValue)
	: name(newName) {
	std::ostringstream valueOut;
	valueOut << std::setprecision(16) << newValue;
	value = valueOut.str();
	

	  
      }
	
    /**
      * Comparison operators for attributes
      */
    inline bool operator==(const Attribute& other) {
      return ( name == other.name &&
      value == other.value);
    }
    inline bool operator!=(const Attribute& other) {
      return !(*this == other);
    }
    };

    /**
     * Append an attribute to this element
     */
    Element& operator<<(const Attribute&);

    /**
     * Append an attribute to this element
     */
    template<class T>
    void appendAttribute(const std::string& name, const T& t) {
      *this << Attribute(name,t);
    }

    /**
     * Get a value from an attribute
     */
    template<class T>
    void getFromAttribute(const std::string& name, T& t) const {
      std::istringstream in(attribute(name));
      in >> std::setprecision(16) >> t;
    }

  public:

    /**
     * Return true, if this element has children
     */
    bool hasChildren() const {
      return 
	type() == ElementTypes::Root ||
	type() == ElementTypes::Element;
    }

    /**
     * Return the list of children elements
     */
    const std::list<Element>& children() const;

    /**
     * Access the list of children elements
     */
    std::list<Element>& children();

    /**
     * Append a child element to this element
     */
    Element& append(const Element&);

    /**
     * Prepend a child element to this element
     */
    Element& prepend(const Element&);

    /**
     * Append a child element to this element
     */
    Element& operator<<(const Element&);

    /**
     * Insert an element before the given position
     */
    void insert(std::list<Element>::iterator, const Element&);

    /**
     * Erase an element at the given position
     */
    void erase(std::list<Element>::iterator);

    /**
     * Find the first element of the given type and name
     */
    std::list<Element>::const_iterator 
    findFirst(int type, const std::string& name) const;

    /**
     * Find all elements of the given type and name
     */
    std::pair<std::multimap<std::pair<int,std::string>,std::list<Element>::iterator>::const_iterator,
	      std::multimap<std::pair<int,std::string>,std::list<Element>::iterator>::const_iterator>
    findAll(int type, const std::string& name) const;

  private:

    /**
     * The type of this element
     */
    int theType;

    /**
     * The name or character content of the element
     */
    std::string theNameOrContent;

    /**
     * The attributes of this element
     */
    std::map<std::string,std::string> theAttributes;

    /**
     * The children elements of this element
     */
    std::list<Element> theChildren;

    /**
     * Index children elements by type and name
     */
    std::multimap<std::pair<int,std::string>,std::list<Element>::iterator> theIndex;

  private:

    /**
     * Index elements by type and name
     */
    void index();

    /**
     * Assert that this element got a character content
     */ 
    void assertContent() const;

    /**
     * Assert that this element is named
     */
    void assertNamed() const;

    /**
     * Assert that this element contains attributes
     */
    void assertAttributes() const;

    /**
     * Assert that this element contains children elements
     */
    void assertChildren() const;

  };

}

#endif // MYXML_Element_hpp_included
