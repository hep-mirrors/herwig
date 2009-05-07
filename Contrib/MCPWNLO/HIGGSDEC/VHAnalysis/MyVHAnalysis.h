// -*- C++ -*-
//
// MyVHAnalysis.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MyVHAnalysis_H
#define HERWIG_MyVHAnalysis_H
//
// This is the declaration of the MyVHAnalysis class.
//

#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig++/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MyVHAnalysis class is designed to perform some simple analysis of
 * gauge boson, W and Z, distributions in hadron-hadron collisions. The main 
 * distriubtions are the transverse momentum and rapidities of the gauge bosons
 * which are of physical interest, and the azimuthal angle distribution for
 * testing purposes.
 *
 * @see \ref MyVHAnalysisInterfaces "The interfaces"
 * defined for MyVHAnalysis.
 */
class MyVHAnalysis: public AnalysisHandler {

public:

  /**
   * The default constructor.
   */
  MyVHAnalysis();

  /** @name Virtual Functions required by the AnalysisHandler class. */
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
  inline virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

protected:

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
  static NoPIOClassDescription<MyVHAnalysis> initMyVHAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MyVHAnalysis & operator=(const MyVHAnalysis &);

private:

  /**
   *   Mass of the b bbar pair
   */
  Histogram _Mbb;

  /**
   *   Energy of the b bbar pair
   */
  Histogram _Ebb;

  /**
   *   pT of the b bbar pair wrt the beam axis
   */
  Histogram _PTbmH;

  /**
   *  pT of the b bbar pair wrt the vector boson
   */
  Histogram _PTVH;

  /**
   * Longitudinal momentum of the b bbar pair
   */
  Histogram _Lbb;

  /**
   *  Rapidity of the b bbar pair
   */
  Histogram _rapbb;

  
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MyVHAnalysis. */
template <>
struct BaseClassTrait<Herwig::MyVHAnalysis,1> {
  /** Typedef of the first base class of MyVHAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MyVHAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MyVHAnalysis>
  : public ClassTraitsBase<Herwig::MyVHAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MyVHAnalysis"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the MyVHAnalysis class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "MyVHAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MyVHAnalysis_H */
