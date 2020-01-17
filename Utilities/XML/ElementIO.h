// -*- C++ -*-
//
// ElementIO.hpp is a part of myXML
// Copyright (C) 2012-2019 Simon Platzer, The Herwig Collaboration
//
// myXML is licenced under version 3 of the GPL, see COPYING for details.
//

#ifndef MYXML_ElementIO_hpp_included
#define MYXML_ElementIO_hpp_included

#include <iostream>

#include "Element.h"

namespace XML {

  /**
   * \brief ElementIO handles in/output of XML elements
   * \author Simon Platzer
   */
  class ElementIO {

  public:

    /**
     * Write a tree to an ostream
     */
    template<class OStream>
    static void put(Element, OStream&);

    /**
     * Get the next element from a stream
     */
    static Element get(std::istream&);

    /**
     * Parse the entire stream
     */
    static Element getAll(std::istream&);

  private:

    /**
     * Strip whitespaces from a string
     */
    static void strip(std::string&, const std::string& skip = "\n\t ");

    /**
     * Read from istream until occurence of the given pattern
     */
    static std::istream& getline(std::istream&, std::string&, const std::string&);

    /**
     * Skip charecters encountered on the given istream
     */
    static void skip(std::istream&, const std::string& skip_chars = "\n\t ");

    /**
     * Helper struct representing a single tag or parsed content
     */
    struct Tag {

      /**
       * The tag type enumeration
       */
      enum EnumerateTagTypes {

	Unknown = -1,
	/** Tag type unknown */
	EmptyElement = 1,
	/** An empty element tag */
	ElementBegin = 2,
	/** An element begin tag */
	ProcessingInstruction = 3,
	/** A processing instruction tag */
	CharacterData = 4,
	/** Character data tag */
	ParsedCharacterData = 5,
	/** Parsed character data */
	Comment = 6,
	/** A comment tag */
	ElementEnd = 20
	/** An element end tag */
      };

      /**
       * The type of the tag
       */
      int type;

      /**
       * The content or name of the tag
       */
      std::string content;

      /**
       * A lis of attributes, if present
       */
      std::map<std::string,std::string> attributes;

      /**
       * Produce an element
       */
      Element produce() const;

    };

    /**
     * Get a single tag or parsed content from an istream 
     */
    static void getTag(Tag&, std::istream&);

    /**
     * Produce element tree from parsed stack of tags
     */
    static Element produce(std::list<Tag>&);

    /**
     * Produce element tree from parsed stack of tags
     */
    static void produce(std::list<Tag>&, Element&);

    /**
     * Throw an unexpected end of file exception
     */
    static void unexpectedEOF(const std::string& what = "");

    /**
     * Throw a parse error exception
     */
    static void parseError(const std::string& what = "");

  };

}

#include "ElementIO.tcc"

#endif // MYXML_ElementIO_hpp_included
