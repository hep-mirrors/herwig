// -*- C++ -*-
//
// ProjectingSampler.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_ProjectingSampler_H
#define Herwig_ProjectingSampler_H
//
// This is the declaration of the ProjectingSampler class.
//

#include "BinSampler.h"
#include "BinnedStatistics.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief ProjectingSampler does adaption from projections of the integrand.
 *
 * @see \ref ProjectingSamplerInterfaces "The interfaces"
 * defined for ProjectingSampler.
 */
class ProjectingSampler: public Herwig::BinSampler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ProjectingSampler();

  /**
   * The destructor.
   */
  virtual ~ProjectingSampler();
  //@}

public:

  /**
   * Generate the next point; store the point in lastPoint() and its
   * weight using select(); if noMaxInfo is true, do not throw
   * NewMaximum or UpdateCrossSections exceptions.
   */
  virtual void generate(bool noMaxInfo = false);

  /**
   * Initialize this bin sampler. This default version calls runIteration.
   */
  virtual void initialize(bool progress);

  /**
   * Finish an iteration, performing the adaption.
   */
  void adapt();

public:

  /**
   * Select an event
   */
  virtual void select(double weight);

  /**
   * Accept an event.
   */
  virtual void accept();

  /**
   * Reject an event.
   */
  virtual void reject();

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
  ProjectingSampler & operator=(const ProjectingSampler &);

  /**
   * True, if we're running the first iteration.
   */
  bool theFirstIteration;

  /**
   * The number of iterations to be considered for initialization.
   */
  unsigned long theNIterations;

  /**
   * Factor to enhance the number of points for the next iteration.
   */
  double theEnhancementFactor;

  /**
   * The initial number of bins to use.
   */
  unsigned int theNBins;

  /**
   * The adaption threshold.
   */
  double theEpsilon;

  /**
   * The number of points used for the last iteration.
   */
  unsigned long theLastNPoints;

  /**
   * The projections to use.
   */
  vector<BinnedStatistics> theProjections;

  /**
   * The last integrand value.
   */
  double theLastValue;

  /**
   * The weight threshold which governs the minimum bin weight.
   */
  double theWeightThreshold;

};

}

#endif /* Herwig_ProjectingSampler_H */
