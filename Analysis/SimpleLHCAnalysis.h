// -*- C++ -*-
//
// SimpleLHCAnalysis.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SimpleLHCAnalysis_H
#define HERWIG_SimpleLHCAnalysis_H
//
// This is the declaration of the SimpleLHCAnalysis class.
//

#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Analysis
 * The SimpleLHCAnalysis class is designed to perform some simple analysis of
 * gauge boson, W and Z, distributions in hadron-hadron collisions. The main 
 * distriubtions are the transverse momentum and rapidities of the gauge bosons
 * which are of physical interest, and the azimuthal angle distribution for
 * testing purposes.
 *
 * @see \ref SimpleLHCAnalysisInterfaces "The interfaces"
 * defined for SimpleLHCAnalysis.
 */
class SimpleLHCAnalysis: public AnalysisHandler {

public:

  /**
   * The default constructor.
   */
  SimpleLHCAnalysis();

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
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
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
  static NoPIOClassDescription<SimpleLHCAnalysis> initSimpleLHCAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SimpleLHCAnalysis & operator=(const SimpleLHCAnalysis &) = delete;

private:

  /**
   *   \f$p_T\f$ of the Z boson
   */
  vector<Histogram> _ptZ;

  /**
   *   \f$p_T\f$ of the \f$W^+\f$ boson
   */
  vector<Histogram> _ptWp;

  /**
   *   \f$p_T\f$ of the \f$W^-\f$ boson
   */
  vector<Histogram> _ptWm;

  /**
   * Mass of the Z boson
   */
  Histogram _mZ;

  /**
   * Mass of the \f$W^+\f$ boson
   */
  Histogram _mWp;

  /**
   * Mass of the \f$W^-\f$ boson
   */
  Histogram _mWm;

  /**
   *  Rapidity of Z
   */
  Histogram _rapZ;

  /**
   *  Rapidity of \f$W^+\f$ boson
   */
  Histogram _rapWp;

  /**
   *  Rapidity of \f$W^-\f$ boson
   */
  Histogram _rapWm;

  /**
   *  Azimuth of Z
   */
  Histogram _phiZ;

  /**
   *  Azimuth of \f$W^+\f$ boson
   */
  Histogram _phiWp;

  /**
   *  Azimuth of \f$W^-\f$ boson
   */
  Histogram _phiWm;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SimpleLHCAnalysis. */
template <>
struct BaseClassTrait<Herwig::SimpleLHCAnalysis,1> {
  /** Typedef of the first base class of SimpleLHCAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SimpleLHCAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SimpleLHCAnalysis>
  : public ClassTraitsBase<Herwig::SimpleLHCAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SimpleLHCAnalysis"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the SimpleLHCAnalysis class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_SimpleLHCAnalysis_H */
