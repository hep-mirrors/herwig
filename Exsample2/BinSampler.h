// -*- C++ -*-
//
// BinSampler.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_BinSampler_H
#define Herwig_BinSampler_H
//
// This is the declaration of the BinSampler class.
//

#include "ThePEG/Handlers/StandardEventHandler.h"
#include "ThePEG/Utilities/Exception.h"
#include "ThePEG/Repository/UseRandom.h"

#include "MultiIterationStatistics.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief BinSampler samples XCombs bins. This default implementation
 * performs flat MC integration.
 *
 * @see \ref BinSamplerInterfaces "The interfaces"
 * defined for BinSampler.
 */
class BinSampler: public Herwig::MultiIterationStatistics {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  BinSampler();

  /**
   * The destructor.
   */
  virtual ~BinSampler();
  //@}

public:

  /**
   * Clone this object.
   */
  Ptr<BinSampler>::ptr cloneMe() const {
    return dynamic_ptr_cast<Ptr<BinSampler>::ptr>(clone());
  }

public:

  /**
   * Set the event handler
   */
  void eventHandler(tStdEHPtr eh) { theEventHandler = eh; }

  /**
   * Return the event handler
   */
  tStdEHPtr eventHandler() const { return theEventHandler; }

  /**
   * Return the bin
   */
  int bin() const { return theBin; }

  /**
   * Set the bin
   */
  void bin(int b) { theBin = b; }

  /**
   * Return a string describing the process handled by this sampler.
   */
  string process() const;

  /**
   * Return the last generated point.
   */
  const vector<double>& lastPoint() const { return theLastPoint; }

  /**
   * Access the last generated point.
   */
  vector<double>& lastPoint() { return theLastPoint; }

  /**
   * Return true, if this bin sampler produces unweighted events.
   */
  virtual bool isUnweighting() const { return false; }

  /**
   * Return true, if this sampler is in a compensating mode.
   */
  virtual bool compensating() const { return false; }

  /**
   * Exception to be thrown, if a new maximum weight has been encountered
   * which should be taken care of by the unweighting performed in GeneralSampler
   */
  struct NewMaximum {
    double oldMaxWeight;
    double newMaxWeight;
    NewMaximum(double oldm, double newm)
      : oldMaxWeight(oldm), newMaxWeight(newm) {}
  };

  /**
   * Exception to be thrown if cross section information should be updated.
   */
  struct UpdateCrossSections {};

  /**
   * Generate the next point; store the point in lastPoint() and its
   * weight using select(); if noMaxInfo is true, do not throw
   * NewMaximum or UpdateCrossSections exceptions.
   */
  virtual void generate(bool noMaxInfo = false);

  /**
   * Run a single iteration of n points, optionally printing a
   * progress bar to cout. Calls generate n times.
   */
  void runIteration(unsigned long n, bool progress);

  /**
   * Initialize this bin sampler. This default version calls runIteration.
   */
  virtual void initialize(bool progress);

  /**
   * Return true, if this sampler has already been initialized.
   */
  bool initialized() const { return theInitialized; }

  /**
   * Indicate that this sampler has already been initialized.
   */
  void isInitialized() { theInitialized = true; }

  /**
   * Finalize this sampler.
   */
  virtual void finalize(bool) {}

public:

  /**
   * Return the dimension.
   */
  int dimension() const { return theEventHandler->nDim(bin()); }

  /**
   * Return the number of points to be used for initial integration.
   */
  unsigned long initialPoints() const { return (unsigned long)(theEnhanceInitialPoints*theInitialPoints); }

  /**
   * Set an enhancement factor for the number of initial points.
   */
  void enhanceInitialPoints(double f) { theEnhanceInitialPoints = f; }

  /**
   * Set an oversampling factor for this sampler.
   */
  void oversamplingFactor(double f) { theOversamplingFactor = f; }

  /**
   * Return enhancement factor for the number of initial points.
   */
  double enhanceInitialPoints() const { return theEnhanceInitialPoints; }

  /**
   * Return the oversampling factor for this sampler.
   */
  double oversamplingFactor() const { return theOversamplingFactor; }

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
   * The number of points to use for initial integration.
   */
  unsigned long theInitialPoints;

  /**
   * An enhancement factor for the number of initial points.
   */
  double theEnhanceInitialPoints;

  /**
   * An oversampling factor for this sampler.
   */
  double theOversamplingFactor;

  /**
   * The bin to be sampled.
   */
  int theBin;

  /**
   * Wether or not this sampler has already been initialized.
   */
  bool theInitialized;

  /**
   * The last generated point.
   */
  vector<double> theLastPoint;

protected:

  /**
   * The event handler to be used.
   */
  tStdEHPtr theEventHandler;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  BinSampler & operator=(const BinSampler &);

};

}

#endif /* Herwig_BinSampler_H */
