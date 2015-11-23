// -*- C++ -*-
//
// xmlTests.cc is a part of myXML
// Copyright (C) 2015 Simon Plaetzer, Marco A. Harrendorf
//
// myXML is licenced under version 2 of the GPL, see COPYING for details.
//
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#define BOOST_TEST_MODULE xmlTest

#include "Herwig/Utilities/XML/Element.h"
#include "Herwig/Utilities/XML/ElementIO.h"

/*
 * Fixture which defines the common variables for testing Element class
 */
struct Fix1 {
    Fix1() : elementUnknown(XML::ElementTypes::Unknown,""), 
	     elementRoot(XML::ElementTypes::Root,"Root"),
	     elementEmpty(XML::ElementTypes::EmptyElement,""),
	     elementElement(XML::ElementTypes::Element, "Element"),
	     elementProcessingInstruction(XML::ElementTypes::ProcessingInstruction, "ProcessingInstruction"),
	     elementCharacterData(XML::ElementTypes::CharacterData, "CharacterData"),
	     elementParsedCharacterData(XML::ElementTypes::ParsedCharacterData, "ParsedCharacterData"),
	     elementComment(XML::ElementTypes::Comment, "Comment"),
	     elementGrids(XML::ElementTypes::Element,"Grids")	  
    {BOOST_TEST_MESSAGE( "setup fixture for xmlTest" ); }
    
    ~Fix1()  { BOOST_TEST_MESSAGE( "teardown fixture for xmlTest" ); }

    XML::Element elementUnknown;
    XML::Element elementRoot;
    XML::Element elementEmpty;
    XML::Element elementElement;
    XML::Element elementProcessingInstruction;
    XML::Element elementCharacterData;
    XML::Element elementParsedCharacterData;
    XML::Element elementComment;
    XML::Element elementGrids;
};

//____________________________________________________________________________//

/*
 * Start of boost unit tests
 * @todo Implement tests for ElementIO
 */
BOOST_FIXTURE_TEST_SUITE(xmlTestElement, Fix1 )

/*
 * Boost unit tests for Element.h
 * @todo Implement further tests for child operations
 */
BOOST_AUTO_TEST_CASE(elementTypes)
{
    BOOST_CHECK_EQUAL(elementUnknown.type(), -1);
    BOOST_CHECK_EQUAL(elementRoot.type(), 0);
    BOOST_CHECK_EQUAL(elementEmpty.type(), 1);
    BOOST_CHECK_EQUAL(elementElement.type(), 2);
    BOOST_CHECK_EQUAL(elementGrids.type(), 2);
    BOOST_CHECK_EQUAL(elementProcessingInstruction.type(), 3);
    BOOST_CHECK_EQUAL(elementCharacterData.type(), 4);
    BOOST_CHECK_EQUAL(elementParsedCharacterData.type(), 5);
    BOOST_CHECK_EQUAL(elementComment.type(), 6);
}

BOOST_AUTO_TEST_CASE(constructors)
{
    BOOST_CHECK(elementUnknown == XML::Element());
    BOOST_CHECK(elementComment == XML::Element(XML::ElementTypes::Comment, "Comment"));

    
    BOOST_CHECK(elementUnknown != elementComment);
    
    XML::Element elementCopy = elementComment;
    BOOST_CHECK(elementCopy == elementComment);
    
    XML::Element elementAssignment = XML::Element();
    elementAssignment = elementGrids;
    BOOST_CHECK(elementAssignment == elementGrids);
}

BOOST_AUTO_TEST_CASE(combineOperator)
{
    BOOST_CHECK(elementElement + elementElement == elementElement);
    BOOST_CHECK_THROW((elementElement + elementComment), std::logic_error);
    XML::Element elementName1 = XML::Element(XML::ElementTypes::Element, "Name1") ;
    XML::Element elementName2 = XML::Element(XML::ElementTypes::Element, "Name2") ;
    BOOST_CHECK_THROW(elementElement + elementGrids, std::logic_error);
    
    XML::Element elementWithBothChildren = XML::Element(elementElement);
    elementWithBothChildren.append(elementComment);
    elementWithBothChildren.append(elementParsedCharacterData);
    XML::Element elementWithChild1 = XML::Element(elementElement);
    elementWithChild1.append(elementComment);
    XML::Element elementWithChild2 = XML::Element(elementElement);
    elementWithChild2.append(elementParsedCharacterData);
    BOOST_CHECK(elementWithChild1 + elementWithChild2 == elementWithBothChildren);   
    
    XML::Element elementWithBothAttributes = XML::Element(elementElement);
    elementWithBothAttributes.appendAttribute("Entry1", "Value1");
    elementWithBothAttributes.appendAttribute("Entry2", "Value2");
    XML::Element elementWithAttribute1 = XML::Element(elementElement);
    elementWithAttribute1.appendAttribute("Entry1", "Value1");
    XML::Element elementWithAttribute2 = XML::Element(elementElement);
    elementWithAttribute2.appendAttribute("Entry2", "Value2");
    BOOST_CHECK(elementWithAttribute1 + elementWithAttribute2 == elementWithBothAttributes);  
}

BOOST_AUTO_TEST_CASE(type)
{
    BOOST_CHECK(elementUnknown.type() != elementGrids.type());
    BOOST_CHECK_EQUAL(elementGrids.type(), XML::ElementTypes::Element);
}
 
BOOST_AUTO_TEST_CASE(name)
{
    BOOST_CHECK(elementElement.name() != elementEmpty.name());
    BOOST_CHECK_EQUAL(elementGrids.name(), "Grids");
} 
  
BOOST_AUTO_TEST_CASE(content)
{
    BOOST_CHECK_EQUAL(elementProcessingInstruction.content(), "ProcessingInstruction");
    BOOST_CHECK_EQUAL(elementCharacterData.content(), "CharacterData");
    BOOST_CHECK_EQUAL(elementParsedCharacterData.content(), "ParsedCharacterData");
    BOOST_CHECK_EQUAL(elementComment.content(), "Comment");    
    XML::Element elementAppendContent = XML::Element(elementComment);
    elementAppendContent << std::string("Test");
    BOOST_CHECK_EQUAL(elementAppendContent.content(), "CommentTest");
}
    
BOOST_AUTO_TEST_CASE(attributes)
{    
    BOOST_CHECK(elementEmpty.hasAttributes());
    BOOST_CHECK(elementElement.hasAttributes());
    BOOST_CHECK(elementGrids.hasAttributes());

    
    XML::Element elementEmptyWithAttributes = XML::Element(elementEmpty);
    elementEmptyWithAttributes.appendAttribute("Entry1", "Value1");
    BOOST_CHECK(elementEmptyWithAttributes.hasAttribute("Entry1"));
    
    XML::Element elementElementWithAttributes = XML::Element(elementElement);
    XML::Element::Attribute elementAttribute = XML::Element::Attribute("Entry2","Value2"); 
    elementElementWithAttributes << elementAttribute;
    elementElementWithAttributes.appendAttribute("Entry3","Value3");
    BOOST_CHECK(elementElementWithAttributes.hasAttribute("Entry2"));
    BOOST_CHECK(elementElementWithAttributes.hasAttribute("Entry3"));

    
    std::map<std::string,std::string> comparisonAttributes;
    comparisonAttributes.insert(std::pair<std::string,std::string>("Entry2", "Value2"));
    comparisonAttributes.insert(std::pair<std::string,std::string>("Entry3", "Value3"));
    std::map<std::string,std::string> elementElementWithAttributesMap = elementElementWithAttributes.attributes();
    BOOST_CHECK(comparisonAttributes == elementElementWithAttributesMap);
  
    BOOST_CHECK_EQUAL(elementElementWithAttributes.attribute("Entry3"), "Value3");
 
    std::string elementAttributeOutput;
    elementElementWithAttributes.getFromAttribute("Entry3", elementAttributeOutput);
    BOOST_CHECK_EQUAL(elementAttributeOutput, "Value3");
    
    XML::Element::Attribute elementAttribute4 = XML::Element::Attribute("Entry4","Value4");
    XML::Element::Attribute elementAttribute5 = XML::Element::Attribute("Entry5","Value5"); 
    
    BOOST_CHECK(elementAttribute4 == XML::Element::Attribute("Entry4","Value4"));
    BOOST_CHECK(elementAttribute4 != XML::Element::Attribute("Entry4",""));
    BOOST_CHECK(elementAttribute4 != XML::Element::Attribute("","Value4"));
    BOOST_CHECK(elementAttribute4 != elementAttribute5);   
} 
    
BOOST_AUTO_TEST_CASE(children)
{
    BOOST_CHECK(elementRoot.hasChildren());
    BOOST_CHECK(elementElement.hasChildren());
    BOOST_CHECK(elementGrids.hasChildren());
    
    std::list<XML::Element> listWithoutElement;
    BOOST_CHECK(elementRoot.children() == listWithoutElement);
    BOOST_CHECK(elementElement.children() == listWithoutElement);    
    

    std::list<XML::Element> listWithElement;
    listWithElement.push_back(elementComment);
    listWithElement.push_front(elementParsedCharacterData);
    XML::Element elementWithChildren = XML::Element(elementElement);
    elementWithChildren.append(elementComment);
    elementWithChildren.prepend(elementParsedCharacterData);
    BOOST_CHECK(elementWithChildren.children() == listWithElement);
    
    elementWithChildren << elementGrids;
    listWithElement.push_back(elementGrids);
    BOOST_CHECK(elementWithChildren.children() == listWithElement);
}
    
BOOST_AUTO_TEST_SUITE_END()
