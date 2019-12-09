// -*- C++ -*-
//
// HiggsJetAnalysis.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_HiggsJetAnalysis_H
#define HERWIG_HiggsJetAnalysis_H
//
// This is the declaration of the HiggsJetAnalysis class.
//

#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Analysis
 * HiggsJetAnalysis assumes that there is one Higgs in the final state
 * and books some observables computed from its four momentum.  It
 * shouldn't do anything in case there is no Higgs in the event.
 *
 * @see \ref HiggsJetAnalysisInterfaces "The interfaces"
 * defined for HiggsJetAnalysis.
 */
class HiggsJetAnalysis: public AnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  HiggsJetAnalysis();
  //@}

public:

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

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HiggsJetAnalysis & operator=(const HiggsJetAnalysis &) = delete;

private:
  /**
   *   \f$p_T\f$ of the h boson
   */
  Histogram _pth;
  Histogram _pthZoom;

  /**
   *  Rapidity of h
   */
  Histogram _raph;

  /**
   *  Azimuth of h
   */
  Histogram _phih;

};

}

#endif /* HERWIG_HiggsJetAnalysis_H */
