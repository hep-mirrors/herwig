// -*- C++ -*-
//
// CLEOCharmAnalysis.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_CLEOCharmAnalysis_H
#define HERWIG_CLEOCharmAnalysis_H
//
// This is the declaration of the CLEOCharmAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Analysis
 * The CLEOCharmAnalysis class compares the results of Herwig at 10.52 GeV
 * with data on Charm hadron spectra from the CLEO experiment
 *
 * @see \ref CLEOCharmAnalysisInterfaces "The interfaces"
 * defined for CLEOCharmAnalysis.
 */
class CLEOCharmAnalysis: public AnalysisHandler {

public:

  /**
   * The default constructor.
   */
  CLEOCharmAnalysis() : _s() {}

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

  /**
   * Analyze the given particle.
   * @param particle pointer to the particle to be analyzed.
   * @param weight The weight for the event
   */
  virtual void analyze(tPPtr particle, double weight);

  using AnalysisHandler::analyze;
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
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

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
  static NoPIOClassDescription<CLEOCharmAnalysis> initCLEOCharmAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  CLEOCharmAnalysis & operator=(const CLEOCharmAnalysis &);

private:

  /**
   *  Histograms
   */
  //@{
  /**
   *  Histogram for \f$D^+\f$
   */
  HistogramPtr _histDplus;

  /**
   *  Histogram for \f$D^0\f$
   */
  HistogramPtr _histD0;

  /**
   *  Histogram for \f$D^{*+}\f$
   */
  HistogramPtr _histDstarplus;

  /**
   *  Histogram for \f$D^{*0}\f$
   */
  HistogramPtr _histDstar0;
  //@}

  /**
   *  CMF energy squared
   */
  Energy2 _s;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of CLEOCharmAnalysis. */
template <>
struct BaseClassTrait<Herwig::CLEOCharmAnalysis,1> {
  /** Typedef of the first base class of CLEOCharmAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the CLEOCharmAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::CLEOCharmAnalysis>
  : public ClassTraitsBase<Herwig::CLEOCharmAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::CLEOCharmAnalysis"; }
  /**
   * The name of a file containing the dynamic library where the class
   * CLEOCharmAnalysis is implemented. It may also include several, space-separated,
   * libraries if the class CLEOCharmAnalysis depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwAnalysis.so HwLEPAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_CLEOCharmAnalysis_H */
