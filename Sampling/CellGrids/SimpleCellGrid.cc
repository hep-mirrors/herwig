// -*- C++ -*-
//
// SimpleCellGrid.cpp is a part of ExSample
// Copyright (C) 2012-2019 Simon Platzer, The Herwig Collaboration
//
// ExSample is licenced under version 3 of the GPL, see COPYING for details.
//

#include "SimpleCellGrid.h"

#include <exception>
#include <stdexcept>
#include <sstream>
#include <iomanip>

using namespace ExSample;
using namespace std;

SimpleCellGrid::SimpleCellGrid(const std::vector<double>& newLowerLeft,
			       const std::vector<double>& newUpperRight,
			       bool keepWeightInformation,
			       double newWeight)
  : CellGrid(newLowerLeft,newUpperRight,newWeight),
    theReferenceWeight(0.0) {
  if ( keepWeightInformation )
    weightInformation().resize(newLowerLeft.size());
}

void SimpleCellGrid::split(std::size_t newSplitDimension, double newSplitCoordinate) {
  CellGrid::split(newSplitDimension,newSplitCoordinate);
  weightInformation().clear();
}

CellGrid* SimpleCellGrid::makeInstance() const {
  return new SimpleCellGrid();
}

CellGrid* SimpleCellGrid::makeInstance(const std::vector<double>& newLowerLeft,
				       const std::vector<double>& newUpperRight,
				       double newWeight) const {
  return new SimpleCellGrid(newLowerLeft,newUpperRight,!weightInformation().empty(),newWeight);
}

void SimpleCellGrid::adapt(double gain, double epsilon,
			   set<SimpleCellGrid*>& newCells) {
  if ( !isLeaf() ) {
    firstChild().adapt(gain,epsilon,newCells);
    secondChild().adapt(gain,epsilon,newCells);
    return;
  }
  if ( weightInformation().empty() )
    throw runtime_error("[ExSample::SimpleCellGrid] Cannot adapt without weight information.");
  theReferenceWeight = 0.0;
  double avg =
    (weightInformation().front().first.sumOfWeights +
     weightInformation().front().second.sumOfWeights)/
    (weightInformation().front().first.nPoints +
     weightInformation().front().second.nPoints);
  double maxw =
    std::max(weightInformation().front().first.maxWeight,
	     weightInformation().front().second.maxWeight);
  if ( maxw != 0.0 )
    if ( avg/maxw > epsilon )
      return;
  double max = 0.0;
  size_t maxdim = 0;
  for ( size_t k = 0; k < lowerLeft().size(); ++k ) {
    double first = weightInformation()[k].first.averageWeight();
    double second = weightInformation()[k].second.averageWeight();
    if ( first + second == 0.0 ) {
      continue;
    }
    double diff = abs(first-second)/(first+second); 
    if ( diff >= max ) {
      max = diff;
      maxdim = k;
    }
  }
  if ( max >= gain ) {
    split(maxdim,0.5*(lowerLeft()[maxdim]+upperRight()[maxdim]));
    newCells.insert(&firstChild());
    newCells.insert(&secondChild());
  }
}

void SimpleCellGrid::splitter(size_t dim, int rat) {
    if ( !isLeaf() ) {
    firstChild().splitter( dim, rat);
    secondChild().splitter(dim, rat);
    return;
    }
    SimpleCellGrid* b=this;
    for (int i=rat;i>1;i--){
       b->split(dim,1./i*(b->lowerLeft()[dim]+b->upperRight()[dim]));
       b=&(b->secondChild());
    }
}




void SimpleCellGrid::setWeights() {
  if ( !isLeaf() ) {
    firstChild().setWeights();
    secondChild().setWeights();
    return;
  }
  if ( weightInformation().empty() )
    throw runtime_error("[ExSample::SimpleCellGrid] Cannot set weights without weight information.");
  double avg =
    (weightInformation().front().first.sumOfWeights +
     weightInformation().front().second.sumOfWeights)/
    (weightInformation().front().first.nPoints +
     weightInformation().front().second.nPoints);
  weight(avg);
}

void SimpleCellGrid::updateWeightInformation(const vector<double>& point,
					     double w) {
  if ( !isLeaf() )
    throw runtime_error("[ExSample::SimpleCellGrid] Cannot update weight information of a branching node.");
  for ( std::size_t d = 0; d < weightInformation().size(); ++d ) {
    if ( point[d] < 0.5*(lowerLeft()[d]+upperRight()[d]) )
      weightInformation()[d].first.book(w);
    else
      weightInformation()[d].second.book(w);
  }
}

void SimpleCellGrid::fromXML(const XML::Element& grid) {

  CellGrid::fromXML(grid);
  grid.getFromAttribute("referenceWeight",theReferenceWeight);
  if ( !grid.hasAttribute("keepWeightInformation") )
    return;
  bool keepWeightInformation = false;
  grid.getFromAttribute("keepWeightInformation",keepWeightInformation);
  if ( keepWeightInformation )
    weightInformation().resize(lowerLeft().size());

}

XML::Element SimpleCellGrid::toXML() const {

  XML::Element grid = CellGrid::toXML();
  grid.appendAttribute("keepWeightInformation",!weightInformation().empty());
  grid.appendAttribute("referenceWeight",theReferenceWeight);
  return grid;

}
