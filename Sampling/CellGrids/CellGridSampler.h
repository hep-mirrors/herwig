// -*- C++ -*-
//
// CellGridSampler.hpp is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_CellGridSampler_H
#define Herwig_CellGridSampler_H
//
// This is the declaration of the CellGridSampler class.
//

#include "Herwig/Sampling/BinSampler.h"

#include "SimpleCellGrid.h"


namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief CellGridSampler samples XCombs bins using CellGrids
 *
 * @see \ref CellGridSamplerInterfaces "The interfaces"
 * defined for CellGridSampler.
 */
class CellGridSampler: 
    public Herwig::BinSampler, ExSample::SimpleCellGrid {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  CellGridSampler();

  /**
   * The destructor.
   */
  virtual ~CellGridSampler();
  //@}

public:

  /**
   * Clone this object.
   */
  Ptr<CellGridSampler>::ptr cloneMe() const {
    return dynamic_ptr_cast<Ptr<CellGridSampler>::ptr>(clone());
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
   * Adapt
   */
  virtual void adapt();

  /**
   * Return true, if grid data exists for this sampler.
   */
  virtual bool existsGrid() const;

  /**
   * Save grid data
   */
  virtual void saveGrid() const;
  
  /**
   * The splittings for each dimension befor adaption.
   */
    
  const vector<int>& pre_adaption_splits() const { return the_pre_adaption_splits; }

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
  CellGridSampler & operator=(const CellGridSampler &);

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

  /**
   * Perform splits for all channels
   */
  bool theAllChannelSplits;

  /**
   * Perform unweighting in cells
   */
  bool theUnweightCells;

};

}

#endif /* Herwig_CellGridSampler_H */
