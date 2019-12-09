// -*- C++ -*-
//
// CellGrid.cpp is a part of ExSample
// Copyright (C) 2012-2019 Simon Platzer, The Herwig Collaboration
//
// ExSample is licenced under version 3 of the GPL, see COPYING for details.
//

#include "CellGrid.h"

#include <exception>
#include <stdexcept>
#include <sstream>
#include <iomanip>

using namespace ExSample;
using namespace std;

double CellGrid::volume(const vector<double>& lowerLeft,
			const vector<double>& upperRight) const {
  assert(lowerLeft.size() == upperRight.size());
  double res = 1.0;
  vector<double>::const_iterator upper = upperRight.begin();
  vector<double>::const_iterator lower = lowerLeft.begin();
  for ( ; lower != lowerLeft.end(); ++lower, ++upper ) {
    assert(*lower <= *upper);
    res *= *upper - *lower;
  }
  return res;
}

CellGrid::CellGrid(const vector<double>& newLowerLeft,
		   const vector<double>& newUpperRight,
		   double newWeight)
  : theLowerLeft(newLowerLeft), theUpperRight(newUpperRight),
    theVolumeOrIntegral(volume(newLowerLeft,newUpperRight)),
    theWeight(newWeight) {
  theUpperBoundInclusive.resize(lowerLeft().size(),true);
}

CellGrid::~CellGrid() {
  if ( !isLeaf() ) {
    delete theChildren[0];
    delete theChildren[1];
    theChildren.clear();
  }
}

void CellGrid::boundaries(const std::vector<double>& newLowerLeft,
			  const std::vector<double>& newUpperRight) {
  if ( !lowerLeft().empty() )
    throw runtime_error("[ExSample::CellGrid] Cannot set the boundaries of non-empty grids.");
  theLowerLeft = newLowerLeft;
  theUpperRight = newUpperRight;
  theUpperBoundInclusive.resize(lowerLeft().size(),true);
  theVolumeOrIntegral = volume(newLowerLeft,newUpperRight);
}

CellGrid* CellGrid::makeInstance() const { 
  return new CellGrid();
}

CellGrid* CellGrid::makeInstance(const vector<double>& newLowerLeft,
				 const vector<double>& newUpperRight,
				 double newWeight) const { 
  return new CellGrid(newLowerLeft,newUpperRight,newWeight);
}

void CellGrid::weight(double w) { 
  if ( !isLeaf() )
    throw runtime_error("[ExSample::CellGrid] Cannot set the weight of branching nodes.");
  theWeight = w;
}

void CellGrid::updateIntegral() {
  if ( isLeaf() )
    return;
  firstChild().updateIntegral();
  secondChild().updateIntegral();
  theVolumeOrIntegral = firstChild().integral() + secondChild().integral();
  theWeight = 0.0;
}

pair<size_t,double> CellGrid::splitPoint() const {
  if ( isLeaf() )
    throw runtime_error("[ExSample::CellGrid] Leaf nodes have no splits.");
  pair<size_t,double> res;
  for ( res.first = 0; res.first != firstChild().upperRight().size(); ++res.first ) {
    if ( firstChild().upperRight()[res.first] != secondChild().upperRight()[res.first] ) {
      res.second = firstChild().upperRight()[res.first];
      break;
    }
  }
  assert(res.first < firstChild().upperRight().size());
  return res;
}

const CellGrid& CellGrid::firstChild() const {
  if ( isLeaf() )
    throw runtime_error("[ExSample::CellGrid] Cannot access children of leaf nodes.");
  return *theChildren[0];
}

CellGrid& CellGrid::firstChild() {
  if ( isLeaf() )
    throw runtime_error("[ExSample::CellGrid] Cannot access children of leaf nodes.");
  return *theChildren[0];
}

const CellGrid& CellGrid::secondChild() const {
  if ( isLeaf() )
    throw runtime_error("[ExSample::CellGrid] Cannot access children of leaf nodes.");
  return *theChildren[1];
}

CellGrid& CellGrid::secondChild() {
  if ( isLeaf() )
    throw runtime_error("[ExSample::CellGrid] Cannot access children of leaf nodes.");
  return *theChildren[1];
}

void CellGrid::split(size_t newSplitDimension, double newSplitCoordinate) {
  if ( !isLeaf() )
    throw runtime_error("[ExSample::CellGrid] Cannot split an already branching node.");
  if ( newSplitDimension > lowerLeft().size() )
    throw runtime_error("[ExSample::CellGrid] Cannot split along non-existing dimension.");
  assert(lowerLeft()[newSplitDimension] <= newSplitCoordinate && 
	 newSplitCoordinate <= upperRight()[newSplitDimension]);
  vector<double> firstUpperRight = upperRight();
  firstUpperRight[newSplitDimension] = newSplitCoordinate;
  vector<double> secondLowerLeft = lowerLeft();
  secondLowerLeft[newSplitDimension] = newSplitCoordinate;
  theChildren.resize(2);
  theChildren[0] = makeInstance(lowerLeft(),firstUpperRight);
  theChildren[1] = makeInstance(secondLowerLeft,upperRight());
  firstChild().theUpperBoundInclusive = upperBoundInclusive();
  firstChild().theUpperBoundInclusive[newSplitDimension] = false;
  secondChild().theUpperBoundInclusive = upperBoundInclusive();
}

void CellGrid::splitCoordinates(size_t dimension, set<double>& coordinates) const {
  if ( dimension > lowerLeft().size() )
    throw runtime_error("[ExSample::CellGrid] Cannot get splits for non-existing dimension.");
  if ( isLeaf() ) {
    coordinates.insert(lowerLeft()[dimension]);
    coordinates.insert(upperRight()[dimension]);
    return;
  }
  firstChild().splitCoordinates(dimension,coordinates);
  secondChild().splitCoordinates(dimension,coordinates);
}

bool CellGrid::contains(const vector<double>& point,
			const vector<bool>& parameterFlags) const {
  assert(point.size()==parameterFlags.size());
  assert(point.size()==lowerLeft().size());
  for ( size_t k = 0; k < point.size(); ++k ) {
    if ( !parameterFlags[k] )
      continue;
    if ( !upperBoundInclusive()[k] ) {
      if ( point[k] < lowerLeft()[k] || point[k] >= upperRight()[k] )
	return false;
    } else {
      if ( point[k] < lowerLeft()[k] || point[k] > upperRight()[k] )
	return false;
    }
  }
  return true;
}

double CellGrid::nonParametricVolume(const vector<double>& point,
				     const vector<bool>& parameterFlags) const {
  assert(point.size()==parameterFlags.size());
  assert(point.size()==lowerLeft().size());
  double v = 1.0;
  for ( size_t k = 0; k < point.size(); ++k ) {
    if ( parameterFlags[k] )
      continue;
    v *= upperRight()[k] - lowerLeft()[k];
  }
  return v;
}

void CellGrid::updateIntegral(const vector<double>& point,
			      const vector<bool>& parameterFlags,
			      vector<bool>::iterator hashPosition) {
  if ( !contains(point,parameterFlags) ) {
    theVolumeOrIntegral = 0.0;
    *hashPosition = false;
    return;
  }
  if ( isLeaf() ) {
    theVolumeOrIntegral = nonParametricVolume(point,parameterFlags);
    *hashPosition = true;
    return;
  }
  firstChild().updateIntegral(point,parameterFlags,hashPosition+1);
  secondChild().updateIntegral(point,parameterFlags,hashPosition+firstChild().size()+1);
  theVolumeOrIntegral = firstChild().integral() + secondChild().integral();
  theWeight = 0.0;
  *hashPosition = true;
}

void CellGrid::doMinimumSelection(double r,
				  double ref) {
  if ( isLeaf() ) {
    theWeight = ref/theVolumeOrIntegral;
    return;
  }
  double refFirst = 
    integral() != 0.0 ? firstChild().integral()*ref/integral() : 0.5*ref;
  double refSecond = 
    integral() != 0.0 ? secondChild().integral()*ref/integral() : 0.5*ref;
  if ( refFirst/ref < r &&
       refSecond/ref < r ) {
    refFirst = 0.5*ref;
    refSecond = 0.5*ref;
  } else if ( refFirst/ref < r &&
	      refSecond/ref >= r ) {
    refFirst = r*ref;
    refSecond = (1.-r)*ref;
  } else if ( refFirst/ref >= r &&
	      refSecond/ref < r ) {
    refFirst = (1.-r)*ref;
    refSecond = r*ref;
  }
  theVolumeOrIntegral = refFirst + refSecond;
  firstChild().doMinimumSelection(r,refFirst);
  secondChild().doMinimumSelection(r,refSecond);
}

void CellGrid::minimumSelection(double p) {
  updateIntegral();
  double ref = integral();
  if ( ref == 0.0 )
    ref = 1.0;
  doMinimumSelection(p,ref);
}

double CellGrid::volume() const { 
  if ( !isLeaf() )
    throw runtime_error("[ExSample::CellGrid] No volume is stored for branching nodes.");
  return theVolumeOrIntegral;
}

double CellGrid::integral() const { 
  if ( !isLeaf() )
    return theVolumeOrIntegral;
  return theWeight*theVolumeOrIntegral;
}

bool CellGrid::active() const { 
  return theVolumeOrIntegral != 0.0;
}

size_t CellGrid::depth() const {
  if ( !isLeaf() )
    return max(firstChild().depth(),secondChild().depth()) + 1;
  return 0;
}

size_t CellGrid::size() const {
  if ( !isLeaf() )
    return firstChild().size() + secondChild().size() + 1;
  return 1;
}

double CellGrid::projectInterval(const pair<double,double>& interval,
				 size_t dimension) const {
  if ( dimension > lowerLeft().size() )
    throw runtime_error("[ExSample::CellGrid] Cannot project to non-existing dimension.");
  if ( (interval.first <= lowerLeft()[dimension] &&
	interval.second <= lowerLeft()[dimension]) ||
       (interval.first >= upperRight()[dimension] &&
	interval.second >= upperRight()[dimension]) )
    return 0.0;
  if ( interval.first >= lowerLeft()[dimension] &&
       interval.first <= upperRight()[dimension] &&
       interval.second >= lowerLeft()[dimension] &&
       interval.second <= upperRight()[dimension] ) {
    if ( !isLeaf() ) {
      double res = 0.0;
      if ( firstChild().active() )
	res += firstChild().projectInterval(interval,dimension);
      if ( secondChild().active() )
	res += secondChild().projectInterval(interval,dimension);
      return res;
    }
    double res = integral();
    res /= upperRight()[dimension]-lowerLeft()[dimension];
    return res;
  } else {
    throw runtime_error("[ExSample::CellGrid] Integration interval needs to fully be contained in the grid.");
  }
}

map<pair<double,double>,double> 
CellGrid::project(pair<double,double> interval,
		  size_t dimension) const {
  set<double> splits;
  splitCoordinates(dimension,splits);
  assert(splits.size() > 1);
  while ( *splits.begin() < interval.first ) {
    splits.erase(splits.begin());
    if ( splits.empty() )
      break;
  }
  while ( *(--splits.end()) > interval.second ) {
    splits.erase(--splits.end());
    if ( splits.empty() )
      break;
  }
  splits.insert(interval.first);
  splits.insert(interval.second);
  map<pair<double,double>,double> res;
  while ( splits.size() > 1 ) {
    interval.first = *splits.begin();
    interval.second = *(++splits.begin());
    res[interval] = projectInterval(interval,dimension);
    splits.erase(splits.begin());
  }
  return res;
}

void CellGrid::fromXML(const XML::Element& grid) {

  size_t dimension = 0;
  bool branching = false;

  grid.getFromAttribute("dimension",dimension);
  grid.getFromAttribute("branching",branching);

  if ( branching ) {
    grid.getFromAttribute("integral",theVolumeOrIntegral);
  } else {
    grid.getFromAttribute("volume",theVolumeOrIntegral);
    grid.getFromAttribute("weight",theWeight);
  }

  theLowerLeft.resize(dimension);
  theUpperRight.resize(dimension);
  theUpperBoundInclusive.resize(dimension);

  list<XML::Element>::const_iterator cit;

  cit = grid.findFirst(XML::ElementTypes::Element,"Boundaries");
  if ( cit == grid.children().end() )
    throw runtime_error("[ExSample::CellGrid] Expected a Boundaries element.");

  const XML::Element& boundaries = *cit;
  cit = boundaries.findFirst(XML::ElementTypes::ParsedCharacterData,"");
  if ( cit == boundaries.children().end() )
    throw runtime_error("[ExSample::CellGrid] Expected boundary data.");
  istringstream bdata(cit->content());

  for ( size_t k = 0; k < lowerLeft().size(); ++k ) {
    bdata >> theLowerLeft[k] >> theUpperRight[k];
    bool x; bdata >> x;
    theUpperBoundInclusive[k] = x;
  }

  if ( branching ) {

    theChildren.resize(2,0);
    theChildren[0] = makeInstance();
    theChildren[1] = makeInstance();

    cit = grid.findFirst(XML::ElementTypes::Element,"FirstChild");
    if ( cit == grid.children().end() )
      throw runtime_error("[ExSample::CellGrid] Expected a FirstChild element.");
    const XML::Element& first = *cit;
    cit = first.findFirst(XML::ElementTypes::Element,"CellGrid");
    if ( cit == first.children().end() )
      throw runtime_error("[ExSample::CellGrid] Expected a CellGrid element.");
    firstChild().fromXML(*cit);

    cit = grid.findFirst(XML::ElementTypes::Element,"SecondChild");
    if ( cit == grid.children().end() )
      throw runtime_error("[ExSample::CellGrid] Expected a SecondChild element.");
    const XML::Element& second = *cit;
    cit = second.findFirst(XML::ElementTypes::Element,"CellGrid");
    if ( cit == second.children().end() )
      throw runtime_error("[ExSample::CellGrid] Expected a CellGrid element.");
    secondChild().fromXML(*cit);

  }

}

XML::Element CellGrid::toXML() const {

  XML::Element grid(XML::ElementTypes::Element,"CellGrid");
  grid.appendAttribute("dimension",lowerLeft().size());
  grid.appendAttribute("branching",!isLeaf());
  if ( !isLeaf() ) {
    ostringstream i; i << setprecision(17) << integral();
    grid.appendAttribute("integral",i.str());
  } else {
    ostringstream v; v << setprecision(17) << volume();
    grid.appendAttribute("volume",v.str());
    ostringstream w; w << setprecision(17) << weight();
    grid.appendAttribute("weight",w.str());
  }

  XML::Element boundaries(XML::ElementTypes::Element,"Boundaries");

  ostringstream bdata;
  bdata << setprecision(17);
  for ( size_t k = 0; k < lowerLeft().size(); ++k )
    bdata << lowerLeft()[k] << " " << upperRight()[k] << " "
	  << upperBoundInclusive()[k] << " ";

  XML::Element belem(XML::ElementTypes::ParsedCharacterData,bdata.str());
  boundaries.append(belem);

  grid.append(boundaries);

  if ( !isLeaf() ) {

    XML::Element first(XML::ElementTypes::Element,"FirstChild");
    first.append(firstChild().toXML());
    grid.append(first);

    XML::Element second(XML::ElementTypes::Element,"SecondChild");
    second.append(secondChild().toXML());
    grid.append(second);

  }

  return grid;

}

#ifdef units

void CellGrid::dumpToC(ostream& os,
		       const string& name) const {
  os << "double " << name << "::evaluate(const vector<double>& p) const {\n";
  dumpPartToC(os,"  ");
  os << "}\n";
}

void CellGrid::dumpPartToC(ostream& os,
			   string prefix) const {
  if ( isLeaf() ) {
    os << prefix << "return " << weight() << ";\n";
    return;
  }  
  pair<size_t,double> sp = splitPoint();
  os << prefix << "if ( p[" << sp.first << "] < " << sp.second << " ) {\n";
  firstChild().dumpPartToC(os,prefix+"  ");
  os << prefix << "} else {\n";
  secondChild().dumpPartToC(os,prefix+"  ");
  os << prefix << "}\n";
}

#endif

