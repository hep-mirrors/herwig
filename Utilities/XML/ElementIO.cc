// -*- C++ -*-
//
// ElementIO.cpp is a part of myXML
// Copyright (C) 2012-2019 Simon Platzer, The Herwig Collaboration
//
// myXML is licenced under version 3 of the GPL, see COPYING for details.
//

#include "ElementIO.h"

#include <exception>
#include <stdexcept>

using namespace XML;
using namespace std;

void ElementIO::strip(string& str, 
		      const string& skip) {
  string::size_type i = str.find_first_not_of(skip);
  if (i != string::npos) str = str.substr(i);
  i = str.find_last_not_of(skip);
  str = str.substr(0, i + 1);
}

istream& ElementIO::getline(istream& is, string& res, const string& pattern) {
    
  res.clear();

  if (is.eof() || is.fail() || pattern == "") {
    return is;
  }

  string buffer;

  if (pattern.size() == 1) {
    std::getline(is,res,*pattern.begin());
    return is;
  }

  char first = *pattern.begin();
  char last = *(--pattern.end());
  string remainder_pattern = pattern.substr(1,pattern.size()-2);
  string between;

  while (!is.eof() && !is.fail()) {
	
    std::getline(is,buffer,first);

    res += buffer;

    if (is.eof() || is.fail()) {
      return is;
    }

    std::getline(is,between,last);

    if (is.eof() || is.fail()) {
      res += string(1,first) + between;
      return is;
    }

    // pattern needs to be at the end
    string::size_type pos = between.find(remainder_pattern,between.size()-remainder_pattern.size());

    if (pos != string::npos) {

      if (between.size() != remainder_pattern.size())
	res += string(1,first) + between.substr(0,pos-1);
      return is;
    }

    res += string(1,first) + between + string(1,last);

  }

  return is;

}

void ElementIO::skip(istream& is, const string& skip_chars) {
  char check;
  is.get(check);

  while (!is.eof() && !is.fail()) {
    if (skip_chars.find_first_of(check) == string::npos)
      break;
    is.get(check);
  }

  is.putback(check);
}

void ElementIO::getTag(ElementIO::Tag& t, istream& is) {

  t.type = Tag::Unknown;
  t.content.clear();
  t.attributes.clear();

  skip(is);

  if (is.eof() || is.fail())
    unexpectedEOF();

  char next;

  is.get(next);

  if (is.eof() || is.fail())
    unexpectedEOF();

  if (next != '<') {
    is.putback(next);
    t.type = Tag::ParsedCharacterData;
    getline(is,t.content,"<");
    strip(t.content);
    if (is.eof() || is.fail())
      unexpectedEOF();
    is.putback('<');
    return;
  }

  is.get(next);

  if (next == '?') {
    t.type = Tag::ProcessingInstruction;
    getline(is,t.content,"?>");
    strip(t.content);
    if (is.eof() || is.fail())
      unexpectedEOF();
    return;
  }

  if (next == '!') {

    is.get(next);

    if (is.eof() || is.fail())
      unexpectedEOF();

    if (next == '-') {

      is.get(next);

      if (is.eof() || is.fail())
	unexpectedEOF();

      if (next == '-') {
	t.type = Tag::Comment;
	getline(is,t.content,"-->");
	strip(t.content);
	if (is.eof() || is.fail())
	  unexpectedEOF();
	return;
      }

      parseError();

    }

    if (next == '[') {

      char check_buf[7];
      is.get(check_buf,7);

      if (is.eof() || is.fail())
	unexpectedEOF();

      if (string(check_buf) == "CDATA[") {
	t.type = Tag::CharacterData;
	getline(is,t.content,"]]>");
	strip(t.content);
	if (is.eof() || is.fail())
	  unexpectedEOF();
	return;
      }

    }

    throw runtime_error("[XML::ElementIO] unsupported XML content.");

  }

  if (next == '/') {
    t.type = Tag::ElementEnd;
    getline(is,t.content,">");
    strip(t.content);
    if (is.eof() || is.fail())
      unexpectedEOF();
    return;
  }

  is.putback(next);

  // start or emptyelem otherwise

  string buffer;

  getline(is,buffer,">");

  if (is.eof() || is.fail())
    unexpectedEOF();

  if (!buffer.size())
    parseError();

  if (*(--buffer.end()) == '/') {
    t.type = Tag::EmptyElement;
    buffer.erase(--buffer.end());
  } else {
    t.type = Tag::ElementBegin;
  }

  string::size_type pos = buffer.find_first_of("\n\t ");

  t.content = buffer.substr(0,pos);

  if (pos != string::npos) {
    buffer = buffer.substr(pos);
  } else {
    buffer.clear();
  }
  strip(buffer);

  // now parse attributes

  string name;
  string delim;
  string value;

  while (buffer != "") {

    name = buffer.substr(0,buffer.find_first_of("="));
    strip(name);

    buffer = buffer.substr(buffer.find_first_of("=")+1);
    strip(buffer);
    
    if (buffer == "")
      parseError();

    delim = buffer.substr(0,1);
    buffer.erase(buffer.begin());

    value = buffer.substr(0,buffer.find_first_of(delim));
    buffer = buffer.substr(buffer.find_first_of(delim)+1);
    strip(buffer);
    strip(value);

    t.attributes.insert(make_pair(name,value));

  }

}

Element ElementIO::Tag::produce() const {
  Element res(type,content);
  for ( map<string,string>::const_iterator a = attributes.begin();
	a != attributes.end(); ++a )
    res << Element::Attribute(a->first,a->second);
  return res;
}

Element ElementIO::produce(list<Tag>& tagStack) {

  if ( tagStack.front().type == Tag::ElementEnd )
    parseError("element started with closing tag " + tagStack.front().content);

  Element res = tagStack.front().produce();
  tagStack.pop_front();

  if ( tagStack.empty() )
    return res;

  if ( res.type() != ElementTypes::Element )
    parseError();

  produce(tagStack,res);

  return res;

}

void ElementIO::produce(list<Tag>& tagStack, Element& parent) {

  if ( tagStack.empty() )
    parseError();

  while ( tagStack.front().type != Tag::ElementEnd ) {

    Element nextElement = tagStack.front().produce();
    tagStack.pop_front();
    Element& next = parent.append(nextElement);
    if ( next.type() == ElementTypes::Element )
      produce(tagStack,next);

  }

  tagStack.pop_front();

}

Element ElementIO::get(std::istream& is) {

  list<Tag> tagStack;
  list<string> openClose;

  Tag next;

  do {

    getTag(next,is);
    tagStack.push_back(next);

    if ( next.type == Tag::ElementBegin ) {
      openClose.push_back(next.content);
    } else if ( next.type == Tag::ElementEnd ) {
      if ( openClose.empty() )
	parseError();
      if ( openClose.back() != next.content )
	parseError("element " + openClose.back() + " closed by " + next.content);
      openClose.pop_back();
    }

  } while ( !openClose.empty() );

  return produce(tagStack);

}

Element ElementIO::getAll(std::istream& is) {

  Element res(ElementTypes::Root);

  while ( !is.eof() && !is.fail() ) {
    Element next = get(is);
    res.append(next);
  }

  return res;

}

void ElementIO::unexpectedEOF(const string& what) {
  string msg("[XML::ElementIO] unexpected end of file");
  if ( what != "" )
    msg += ": " + what;
  throw runtime_error(msg.c_str());
}

void ElementIO::parseError(const string& what) {
  string msg("[XML::Element] parse error");
  if ( what != "" )
    msg += ": " + what;
  throw runtime_error(msg.c_str());
}

