// -*- C++ -*-
//
// CellGridSampler2.hpp is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_CellGridSampler2_H
#define Herwig_CellGridSampler2_H
//
// This is the declaration of the CellGridSampler2 class.
//

#include "Herwig++/Sampling/BinSampler.h"

#include "SimpleCellGrid.h"

#include "Herwig++/Utilities/XML/Element.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <cstdlib>
namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Johannes Bellm,Simon Platzer,Michael Rauch
 *
 * \brief CellGridSampler2 samples XCombs bins using CellGrids and Monaco
 *
 * @see \ref CellGridSampler2Interfaces "The interfaces"
 * defined for CellGridSampler2.
 */
class CellGridSampler2: 
    public Herwig::BinSampler, ExSample::SimpleCellGrid {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  CellGridSampler2();

  /**
   * The destructor.
   */
  virtual ~CellGridSampler2();
  //@}

public:

  /**
   * Clone this object.
   */
  Ptr<CellGridSampler2>::ptr cloneMe() const {
    return dynamic_ptr_cast<Ptr<CellGridSampler2>::ptr>(clone());
  }

public:

  /**
   * Generate the next point; store the point in lastPoint() and its
   * weight using select(); if noMaxInfo is true, do not throw
   * NewMaximum or UpdateCrossSections exceptions.
   */
  virtual double generate();

  /**
   * Initialize this bin sampler. This default version calls runIteration.
   */
  virtual void initialize(bool progress);

  /**
   * Finalize this sampler.
   */
  virtual void finalize(bool);

  /**
   * Evaluate the cross section.
   */
  double evaluate(const vector<double>& p);

  /**
   * Adapt
   */
  virtual void adapt();
  
  /**
   * Return the dimension -subtracted the monacoDimensions.
   */
  virtual int dimension() const { return Herwig::BinSampler::dimension()-monacoDimensions(); }
  
  int monacoDimensions() const {return theMonacoDimensions.size();}
  

    
  const vector<int>& pre_adaption_splits() const { return the_pre_adaption_splits; }
  
  
  
  
  
  
  
  
  
  
    /**
   * Fill Monaco grid data from an XML element
   */
  virtual void Monaco_fromXML(const XML::Element&);

  /**
   * Return an XML element for the data of the Monaco grid
   */
  virtual XML::Element Monaco_toXML() const;
  
  
  
    /**
   * Adapt this sampler after an iteration has been run
   */
  virtual void Monaco_adapt();
  
  
  
   
   string  doMonacoDimensions(string in){
    int n;
    std::stringstream ss(in);
    ss >> n;
    theMonacoDimensions.insert(make_pair(n,1.0));
    return in;
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  CellGridSampler2 & operator=(const CellGridSampler2 &);

  

  /**
   * The number of points used to explore a cell
   */
  size_t theExplorationPoints;


  /**
   * The number of exploration steps
   */
  size_t theExplorationSteps;


  /**
   * The adaption threshold.
   */
  double theGain;


  /**
   * The adaption threshold.
   */
  double theEpsilon;


  /**
   * The minimum probability for cell selection.
   */
  double theMinimumSelection;


  /**
   * The splittings for each dimension befor adaption.
   */
  vector<int>  the_pre_adaption_splits;


  /**
   * The number of splits to put into parton luminiosity degrees of
   * freedom.
   */
  int theLuminositySplits;


  /**
   * The number of splits to put into channel degrees of freedom.
   */
  int theChannelSplits;




  
  map<int,double> theMonacoDimensions;

  
  vector<double> lastpoint_cellgrid;

  vector<double> lastpoint_monaco;
  
  
  private:

   /**
    * Rate of grid modification (0 for no modification)
    */
  double theAlpha;

   /**
    * Number of grid divisions per dimension
    */
  size_t theGridDivisions;

   /**
    * Grid boundaries 
    * (first index: dimension of random numbers,
    * second index: dimension of partitions per random number) 
    */
  boost::numeric::ublas::matrix<double> theGrid;

   /**
    * Collected value per grid bin
    */
  boost::numeric::ublas::matrix<double> theGridData;

   private:

   /**
    * Number of points collected in iteration so far
    */
  size_t theIterationPoints;
  
  
  
  
  
  
  
};

}

#endif /* Herwig_CellGridSampler2_H */
