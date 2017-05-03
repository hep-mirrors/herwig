// -*- C++ -*-
//
// CellGrid.hpp is a part of ExSample
// Copyright (C) 2012-2017 Simon Platzer, The Herwig Collaboration
//
// ExSample is licenced under version 3 of the GPL, see COPYING for details.
//

#ifndef EXSAMPLE_CellGrid_hpp_included
#define EXSAMPLE_CellGrid_hpp_included

#include <vector>
#include <utility>
#include <set>
#include <map>
#include <cstdlib>
#include <cassert>
#include <cmath>

#define units

#ifdef units
#include <iostream>
#include <string>
#endif

#include "Herwig/Utilities/XML/Element.h"

namespace ExSample {

  /**
   * Simple helper
   */
  inline double sqr(double x) {
    return x*x;
  }

  /**
   * \brief A binary cell grid
   * \author Simon Platzer
   */
  class CellGrid {

  public:

    /**
     * Default constructor
     */
    CellGrid()
      : theVolumeOrIntegral(0.0), theWeight(0.0) {}

    /**
     * Construct given boundaries and a weight
     */
    CellGrid(const std::vector<double>& newLowerLeft,
	     const std::vector<double>& newUpperRight,
	     double newWeight = 0.0);

    /**
     * Produce a new instance of a cell grid
     */
    virtual CellGrid* makeInstance() const;

    /**
     * Produce a new instance of a cell grid
     */
    virtual CellGrid* makeInstance(const std::vector<double>& newLowerLeft,
				   const std::vector<double>& newUpperRight,
				   double newWeight = 0.0) const;

    /**
     * The destructor
     */
    virtual ~CellGrid();

  public:

    /**
     * Set the boundaries
     */
    void boundaries(const std::vector<double>& newLowerLeft,
		    const std::vector<double>& newUpperRight);

    /**
     * Return the lower left corner of the cell grid
     */
    const std::vector<double>& lowerLeft() const { return theLowerLeft; }

    /**
     * Return the upper right corner of the cell grid
     */
    const std::vector<double>& upperRight() const { return theUpperRight; }

    /**
     * Flag in which dimension this cell is upper bound inclusive
     */
    const std::vector<bool>& upperBoundInclusive() const { return theUpperBoundInclusive; }

    /**
     * Calculate a volume given upper and lower bound
     */
    double volume(const std::vector<double>& lowerLeft,
		  const std::vector<double>& upperRight) const;

    /**
     * Return true, if this is a leaf in the tree
     */
    bool isLeaf() const { return theChildren.empty(); }

    /**
     * Return the depth of this grid
     */
    std::size_t depth() const;

    /**
     * Return the number of nodes contained in this grid
     */
    std::size_t size() const;

    /**
     * Split this cell grid in the given dimension and coordinate, if
     * it is a leaf
     */
    virtual void split(std::size_t newSplitDimension, double newSplitCoordinate);

    /**
     * Return the dimension and coordinate along which the first split
     * of this grid occurs
     */
    std::pair<std::size_t,double> splitPoint() const;

    /**
     * Return the first child
     */
    const CellGrid& firstChild() const;

    /**
     * Access the first child
     */
    CellGrid& firstChild();

    /**
     * Return the second child
     */
    const CellGrid& secondChild() const;

    /**
     * Access the second child
     */
    CellGrid& secondChild();

    /**
     * Fill split coordinates along a given dimension
     */
    void splitCoordinates(std::size_t, std::set<double>&) const;

  public:

    /**
     * Return true, if this grid is active with respect to the last
     * parameter point passed.
     */
    bool active() const;

    /**
     * Return the volume relevant for the last parameter point set.
     */
    double volume() const;

    /**
     * Return the integral relevant for the last parameter point set.
     */
    double integral() const;

    /**
     * Return the weight
     */
    double weight() const { return theWeight; }

    /**
     * Set the weight
     */
    void weight(double w);

    /**
     * Update the integrals
     */
    void updateIntegral();

    /**
     * Ensure a minimum cell selection probability
     */
    void minimumSelection(double p = 0.1);

    /**
     * Return true, if this grid contains the given parameter point
     */
    bool contains(const std::vector<double>& point,
		  const std::vector<bool>& parameterFlags) const;

    /**
     * Get the volume relevant for the dimensions not considered
     * parameters
     */
    double nonParametricVolume(const std::vector<double>& point,
			       const std::vector<bool>& parameterFlags) const;

    /**
     * Set a parameter point and flag which dimensions are considered
     * parameters. Calculate a hash for the relevant subgrid.
     */
    void updateIntegral(const std::vector<double>& point,
			const std::vector<bool>& parameterFlags,
			std::vector<bool>::iterator hashPosition);

    /**
     * Return the projection to the given interval and dimension,
     * provided the interval is not overlapping with more than one
     * cell in the given dimension. Use splitCoordinates() to
     * determine appropriate intervals.
     */
    double projectInterval(const std::pair<double,double>& interval,
			   std::size_t dimension) const;

    /**
     * Return the projection to the given interval and dimension
     */
    std::map<std::pair<double,double>,double> 
    project(std::pair<double,double> interval,
	    std::size_t dimension) const;

  public:

    /**
     * Fill CellGrid data from an XML element
     */
    virtual void fromXML(const XML::Element&);

    /**
     * Return an XML element for the data of this CellGrid
     */
    virtual XML::Element toXML() const;

  private:

    /**
     * Ensure a minimum cell selection probability
     */
    void doMinimumSelection(double r,
			    double ref);

    /**
     * The lower left corner of the cell grid
     */
    std::vector<double> theLowerLeft;

    /**
     * The upper right corner of the cell grid
     */
    std::vector<double> theUpperRight;

    /**
     * Flag in which dimension this cell is upper bound inclusive
     */
    std::vector<bool> theUpperBoundInclusive;

    /**
     * The volume (leafs) or integral (nodes) relevant for the last
     * parameter point set.
     */
    double theVolumeOrIntegral;

    /**
     * The weight
     */
    double theWeight;

    /**
     * The children cell grids
     */
    std::vector<CellGrid*> theChildren;

#ifdef units

  public:

    /**
     * Generate a random grid of given maximum depth
     */
    template<class RndGenerator>
    void randomGrid(RndGenerator& rnd, 
		    std::size_t maxDepth,
		    bool forceSplit = true) {
      if ( maxDepth > 0 ) {
	if ( !forceSplit )
	  if ( rnd.rnd() < 0.5 )
	    return;
	std::size_t dimension = (std::size_t)(std::floor(rnd.rnd()*lowerLeft().size()));
	double point = 0.5*(lowerLeft()[dimension]+upperRight()[dimension]);
	split(dimension,point);
	firstChild().weight(rnd.rnd());
	secondChild().weight(rnd.rnd());
	firstChild().randomGrid(rnd,maxDepth-1,false);
	secondChild().randomGrid(rnd,maxDepth-1,false);
      }
    }

    /**
     * Write out C code corresponding to the cell grid function
     */
    void dumpToC(std::ostream& os,
		 const std::string& name) const;

    /**
     * Write out C code corresponding to the cell grid function
     */
    void dumpPartToC(std::ostream& os,
		     std::string prefix = "") const;

#endif

  };

}

#endif // EXSAMPLE_CellGrid_hpp_included

