// -*- C++ -*-
//
// Element.cpp is a part of myXML
// Copyright (C) 2012-2013 Simon Platzer
//
// myXML is licenced under version 2 of the GPL, see COPYING for details.
//

#include "Element.h"

#include <exception>
#include <stdexcept>

using namespace XML;
using namespace std;

Element::Element(const Element& other) 
  : theType(other.theType), 
    theNameOrContent(other.theNameOrContent),
    theAttributes(other.theAttributes),
    theChildren(other.theChildren) {
  index();
}

Element& Element::operator=(const Element& other) {
  theType = other.theType;
  theNameOrContent = other.theNameOrContent;
  theAttributes = other.theAttributes;
  theChildren = other.theChildren;
  index();
  return *this;
}

namespace XML {
  bool operator==(const XML::Element &one, const XML::Element &two) {
    return(one.theType == two.theType &&
    one.theNameOrContent == two.theNameOrContent &&
    one.theAttributes == two.theAttributes &&
    one.theChildren == two.theChildren
    );
  }

  bool operator!=(const XML::Element &one, const XML::Element &two) {
    return !(one == two);
  }
 
  XML::Element operator+(const XML::Element& one, const XML::Element& two) {
    if ( one.type() != two.type())
      throw logic_error("[XML::Element] Trying to combine elements with different types.");
    if ( one.name() != two.name())
      throw logic_error("[XML::Element] Trying to combine elements with different names.");
    
    XML::Element returnElement = XML::Element(one);
    if(two.hasAttributes()) {
      std::map<std::string,std::string> attributesOfSecond = two.attributes();
      for(std::map<std::string,std::string>::const_iterator attributesofSecondIter=attributesOfSecond.begin(); attributesofSecondIter!=attributesOfSecond.end(); ++attributesofSecondIter) {
	returnElement.appendAttribute(attributesofSecondIter->first, attributesofSecondIter->second);
      }
    }
    if(two.hasChildren()) {
     std::list<Element> childrenOfSecond = two.children();
     for(std::list<Element>::const_iterator childrenOfSecondIter = childrenOfSecond.begin(); childrenOfSecondIter != childrenOfSecond.end(); ++childrenOfSecondIter)
       returnElement.append(*childrenOfSecondIter);
    }
    return returnElement; 
  }
}



const string& Element::name() const {
  assertNamed();
  return theNameOrContent;
}

const string& Element::content() const {
  assertContent();
  return theNameOrContent;
}

string& Element::content() {
  assertContent();
  return theNameOrContent;
}

bool Element::hasAttribute(const string& id) const {
  assertAttributes();
  return theAttributes.find(id) != theAttributes.end();
}

const map<string,string>& Element::attributes() const {
  assertAttributes();
  return theAttributes;
}

const string& Element::attribute(const string& id) const {
  assertAttributes();
  map<string,string>::const_iterator ait = theAttributes.find(id);
  if ( ait == theAttributes.end() )
    throw runtime_error("[XML::Element] no such attribute found.");
  return ait->second;
}

string& Element::attribute(const string& id) {
  assertAttributes();
  string& ret = theAttributes[id];
  return ret;
}

Element& Element::operator<<(const Element::Attribute& a) {
  attribute(a.name) = a.value;
  return *this;
}

const list<Element>& Element::children() const {
  assertChildren();
  return theChildren;
}

list<Element>& Element::children() {
  assertChildren();
  return theChildren;
}

Element& Element::append(const Element& e) {
  assertChildren();
  theChildren.push_back(e);
  index();
  return theChildren.back();
}

Element& Element::prepend(const Element& e) {
  assertChildren();
  theChildren.push_front(e);
  index();
  return theChildren.front();
}

void Element::insert(list<Element>::iterator pos, 
		     const Element& e) {
  assertChildren();
  theChildren.insert(pos,e);
  index();
}

void Element::erase(list<Element>::iterator pos) {
  assertChildren();
  theChildren.erase(pos);
  index();
}

Element& Element::operator<<(const Element& e) {
  append(e);
  return *this;
}

void Element::index() {
  theIndex.clear();
  for ( list<Element>::iterator i = theChildren.begin();
	i != theChildren.end(); ++i ) {
    if ( i->type() == ElementTypes::Element ||
	 i->type() == ElementTypes::EmptyElement )
      theIndex.insert(make_pair(pair<int,string>(i->type(),i->name()),i));
    else
      theIndex.insert(make_pair(pair<int,string>(i->type(),string("")),i));
  }
}

list<Element>::const_iterator 
Element::findFirst(int type, const string& name) const {
  assertChildren();
  typedef
    multimap<pair<int,string>,list<Element>::iterator>::const_iterator
    IndexIterator;
  IndexIterator i = theIndex.find(make_pair(type,name));
  if ( i == theIndex.end() )
    return theChildren.end();
  return i->second;
}

pair<multimap<pair<int,string>,list<Element>::iterator>::const_iterator,
     multimap<pair<int,string>,list<Element>::iterator>::const_iterator>
Element::findAll(int type, const string& name) const {
  assertChildren();
  return theIndex.equal_range(make_pair(type,name));
}

void Element::assertContent() const {
  if ( type() == ElementTypes::ProcessingInstruction ||
       type() == ElementTypes::CharacterData ||
       type() == ElementTypes::ParsedCharacterData ||
       type() == ElementTypes::Comment )
    return;
  throw logic_error("[XML::Element] element has no plain character content.");
}

void Element::assertNamed() const {
  if ( type() == ElementTypes::EmptyElement ||
       type() == ElementTypes::Element )
    return;
  throw logic_error("[XML::Element] element has no name.");
}

void Element::assertAttributes() const {
  if ( type() == ElementTypes::EmptyElement ||
       type() == ElementTypes::Element )
    return;
  throw logic_error("[XML::Element] element has no attributes.");
}

void Element::assertChildren() const {
  if ( type() == ElementTypes::Element ||
       type() == ElementTypes::Root )
    return;
  throw logic_error("[XML::Element] element has no children elements.");
}

