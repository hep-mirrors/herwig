// -*- C++ -*-
//
// ElementIO.tpp is a part of myXML
// Copyright (C) 2012-2017 Simon Platzer, The Herwig Collaboration
//
// myXML is licenced under version 3 of the GPL, see COPYING for details.
//

namespace XML {

  template<class OStream>
  void ElementIO::put(Element e, OStream& os) {

    if ( e.type() == ElementTypes::Unknown )
      return;

    if ( e.type() != ElementTypes::Root ) {

      if ( e.type() == ElementTypes::Element ||
	   e.type() == ElementTypes::EmptyElement ) {
	os << "<" << e.name();
	for ( std::map<std::string,std::string>::const_iterator a =
		e.attributes().begin(); a != e.attributes().end(); ++a ) {
	  os << " " << a->first << "=";
	  std::string delim = "\"";
	  if (a->second.find_first_of(delim) != std::string::npos)
	    delim = "'";
	  os << delim << a->second << delim;
	}
	if ( e.type() == ElementTypes::Element ) {
	  os << ">\n";
	}
	if ( e.type() == ElementTypes::EmptyElement ) {
	  os << "/>\n";
	}
      }

      if ( e.type() == ElementTypes::ProcessingInstruction ) {
	os << "<?" << e.content() << "?>\n";
      }

      if ( e.type() == ElementTypes::CharacterData ) {
	os << "<![CDATA[" << e.content() << "]]>\n";
      }

      if ( e.type() == ElementTypes::ParsedCharacterData ) {
	os << e.content() << "\n";
      }

      if ( e.type() == ElementTypes::Comment ) {
	os << "<!--" << e.content() << "-->\n";
      }

    }

    if ( e.hasChildren() ) {
      for ( std::list<Element>::const_iterator c = e.children().begin();
	    c != e.children().end(); ++c )
	put(*c,os);
    }

    if ( e.type() == ElementTypes::Element ) {
      os << "</" << e.name() << ">\n";
    }

  }

}
