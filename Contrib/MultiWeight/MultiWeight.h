// -*- C++ -*-
//
// MultiWeight.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MultiWeight_H
#define HERWIG_MultiWeight_H
//
// This is the declaration of the MultiWeight class.
//

#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig/Utilities/Histogram.h"



#include <cmath>
#include <fstream>

namespace Herwig {

using namespace ThePEG;

/**
 * The MultiWeight class is designed to perform some simple analysis of
 * gauge boson, W and Z, distributions in hadron-hadron collisions. The main 
 * distriubtions are the transverse momentum and rapidities of the gauge bosons
 * which are of physical interest, and the azimuthal angle distribution for
 * testing purposes.
 *
 * @see \ref MultiWeightInterfaces "The interfaces"
 * defined for MultiWeight.
 */
class MultiWeight: public AnalysisHandler {

public:

  /**
   * The default constructor.
   */
  inline MultiWeight();
  
  /** @name Virtual functions required by the AnalysisHandler class. */
  //@{
  /**
   * Analyze a given Event. Note that a fully generated event
   * may be presented several times, if it has been manipulated in
   * between. The default version of this function will call transform
   * to make a lorentz transformation of the whole event, then extract
   * all final state particles and call analyze(tPVector) of this
   * analysis object and those of all associated analysis objects. The
   * default version will not, however, do anything on events which
   * have not been fully generated, or have been manipulated in any
   * way.
   * @param event pointer to the Event to be analyzed.
   * @param ieve the event number.
   * @param loop the number of times this event has been presented.
   * If negative the event is now fully generated.
   * @param state a number different from zero if the event has been
   * manipulated in some way since it was last presented.
   */
  virtual void analyze(tEventPtr event, long ieve, int loop, int state);


  //@}


public:

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

protected:
  inline virtual void doinitrun();
  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static NoPIOClassDescription<MultiWeight> initMultiWeight;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MultiWeight & operator=(const MultiWeight &);

private:

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MultiWeight. */
template <>
struct BaseClassTrait<Herwig::MultiWeight,1> {
  /** Typedef of the first base class of MultiWeight. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MultiWeight class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MultiWeight>
  : public ClassTraitsBase<Herwig::MultiWeight> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MultiWeight"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the MultiWeight class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "MultiWeight.so"; }
};

/** @endcond */

}


#endif /* HERWIG_SimpleLHCAnalysis_H */
