// -*- C++ -*-
//
// SimpleCellGrid.hpp is a part of ExSample
// Copyright (C) 2012-2013 Simon Platzer
//
// ExSample is licenced under version 2 of the GPL, see COPYING for details.
//

#ifndef EXSAMPLE_SimpleCellGrid_hpp_included
#define EXSAMPLE_SimpleCellGrid_hpp_included

#include "CellGrid.h"

namespace ExSample {

  /**
   * \brief A simple cell grid providing basic adaption and sampling
   * \author Simon Platzer
   */
  class SimpleCellGrid
    : public CellGrid {

  public:

    /**
     * Default constructor
     */
    SimpleCellGrid()
      : CellGrid() {}

    /**
     * Construct given boundaries and a weight
     */
    SimpleCellGrid(const std::vector<double>& newLowerLeft,
		   const std::vector<double>& newUpperRight,
		   bool keepWeightInformation = true,
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
     * Split this cell grid in the given dimension and coordinate, if
     * it is a leaf
     */
    virtual void split(std::size_t newSplitDimension, double newSplitCoordinate);

  public:

    /**
     * Return the first child
     */
    const SimpleCellGrid& firstChild() const {
      return dynamic_cast<const SimpleCellGrid&>(CellGrid::firstChild());
    }

    /**
     * Access the first child
     */
    SimpleCellGrid& firstChild() {
      return dynamic_cast<SimpleCellGrid&>(CellGrid::firstChild());
    }

    /**
     * Return the second child
     */
    const SimpleCellGrid& secondChild() const {
      return dynamic_cast<const SimpleCellGrid&>(CellGrid::secondChild());
    }

    /**
     * Access the second child
     */
    SimpleCellGrid& secondChild() {
      return dynamic_cast<SimpleCellGrid&>(CellGrid::secondChild());
    }

  public:

    /**
     * A simple counter to store information used for adaption
     */
    struct Counter {

      /**
       * Default constructor
       */
      Counter()
	: nPoints(0.0), sumOfWeights(0.0), 
	  sumOfSquaredWeights(0.0) {}

      /**
       * The number of points
       */
      double nPoints;

      /**
       * The sum of weights
       */
      double sumOfWeights;

      /**
       * The sum of squared weights
       */
      double sumOfSquaredWeights;

      /**
       * Book a point
       */
      void book(double weight) {
	nPoints += 1.0;
	sumOfWeights += abs(weight);
	sumOfSquaredWeights += sqr(weight);
      }

      /**
       * Return the average weight
       */
      double averageWeight() const { return nPoints != 0.0 ? sumOfWeights/nPoints : 0.0; }

      /**
       * Return the variance of the weights
       */
      double varianceOfAverage() const {
	return 
	  nPoints > 1.0 ?
	  abs(sumOfSquaredWeights/nPoints - sqr(sumOfWeights/nPoints))/(nPoints-1) : 0.0;
      }

    };

    /**
     * Return weight information for adaption steps
     */
    const std::vector<std::pair<Counter,Counter> >& weightInformation() const { return theWeightInformation; }

    /**
     * Access weight information for adaption steps
     */
    std::vector<std::pair<Counter,Counter> >& weightInformation() { return theWeightInformation; }

    /**
     * Update the weight information for the given point
     */
    virtual void updateWeightInformation(const std::vector<double>& p,
					 double w);

    /**
     * Perform a default adaption step, splitting along the dimension
     * which shows up the largest difference in average weights; if
     * this exceeds gain, perform the split.
     */
    virtual void adapt(double gain, std::set<SimpleCellGrid*>& newCells);

  public:

    /**
     * Sample a point flat in this cell
     */
    template<class RndGenerator>
    void sampleFlatPoint(std::vector<double>& p, 
			 RndGenerator& rnd) const {
      assert(p.size() == lowerLeft().size());
      for ( size_t k = 0; k < p.size(); ++k ) {
	p[k] = lowerLeft()[k] + rnd.rnd()*(upperRight()[k]-lowerLeft()[k]);
      }
    }

    /**
     * Sample a point flat in this cell, keeping parameters fixed
     */
    template<class RndGenerator>
    void sampleFlatPoint(std::vector<double>& p, 
			 const std::vector<bool>& parameterFlags,
			 RndGenerator& rnd) const {
      assert(p.size() == lowerLeft().size());
      for ( size_t k = 0; k < p.size(); ++k ) {
	if ( parameterFlags[k] )
	  continue;
	p[k] = lowerLeft()[k] + rnd.rnd()*(upperRight()[k]-lowerLeft()[k]);
      }
    }

    /**
     * Explore the cell grid, given a number of points to be sampled
     * in each cell; the weights of the cell will contain the maximum
     * weight encountered. If newCells is non-empty explore only these
     * cells, otherwise explore all cells.
     */
    template<class RndGenerator, class Function>
    void explore(std::size_t nPoints,
		 RndGenerator& rnd,
		 Function& f,
		 std::set<SimpleCellGrid*>& newCells) {
      if ( !isLeaf() ) {
	firstChild().explore(nPoints,rnd,f,newCells);
	secondChild().explore(nPoints,rnd,f,newCells);
	return;
      }
      if ( !newCells.empty() ) {
	if ( newCells.find(this) == newCells.end() )
	  return;
      }
      double avg = 0.0;
      std::vector<double> point(lowerLeft().size());
      for ( std::size_t k = 0; k < nPoints; ++k ) {
	sampleFlatPoint(point,rnd);
	double w = f.evaluate(point);
	avg += abs(w);
	updateWeightInformation(point,abs(w));
      }
      weight(avg/nPoints);
    }

    /**
     * Perform a number of exploration and adaption steps. Less steps
     * may be performed, if no split has been encountered while
     * adapting.
     */
    template<class RndGenerator, class Function>
    void exploreAndAdapt(std::size_t nSteps,
			 std::size_t nPoints,
			 double gain,
			 double minProb,
			 RndGenerator& rnd,
			 Function& f) {
      std::set<SimpleCellGrid*> newCells;
      for ( std::size_t step = 0; step < nSteps; ++step ) {
	explore(nPoints,rnd,f,newCells);
	newCells.clear();
	adapt(gain,newCells);
	if ( newCells.empty() )
	  break;
      }
      updateIntegral();
      minimumSelection(minProb);
    }

    /**
     * Select a cell
     */
    template<class RndGenerator>
    SimpleCellGrid* selectCell(RndGenerator& rnd) {
      if ( isLeaf() )
	return this;
      if ( firstChild().active() &&
	   secondChild().active() ) {
	double p = firstChild().integral()/integral();
	if ( rnd.rnd() <= p )
	  return firstChild().selectCell(rnd);
	else
	  return secondChild().selectCell(rnd);
      }
      if ( firstChild().active() &&
	   !secondChild().active() )
	return firstChild().selectCell(rnd);
      else
	return secondChild().selectCell(rnd);
    }

    /**
     * Sample a point and return its weight
     */
    template<class RndGenerator, class Function>
    std::pair<double,double> generate(RndGenerator& rnd,
				      Function& f,
				      std::vector<double>& p) {
      SimpleCellGrid* selected = selectCell(rnd);
      selected->sampleFlatPoint(p,rnd);
      double w = f.evaluate(p);
      selected->updateWeightInformation(p,w);
      return std::make_pair(w,selected->weight());
    }

    /**
     * Sample a point and return its weight
     */
    template<class RndGenerator, class Function>
    std::pair<double,double> generate(RndGenerator& rnd,
				      Function& f,
				      std::vector<double>& p,
				      const std::vector<bool>& parameterFlags) {
      SimpleCellGrid* selected = selectCell(rnd);
      selected->sampleFlatPoint(p,parameterFlags,rnd);
      double w = f.evaluate(p);
      selected->updateWeightInformation(p,w);
      return std::make_pair(w,selected->weight());
    }

  public:

    /**
     * Fill CellGrid data from an XML element
     */
    virtual void fromXML(XML::Element);

    /**
     * Return an XML element for the data of this CellGrid
     */
    virtual XML::Element toXML() const;

  private:

    /**
     * Weight information for adaption steps
     */
    std::vector<std::pair<Counter,Counter> > theWeightInformation;

  };

}

#endif // EXSAMPLE_SimpleCellGrid_hpp_included

