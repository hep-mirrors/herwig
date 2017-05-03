// -*- C++ -*-
//
// GammaGammaAnalysis.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_GammaGammaAnalysis_H
#define HERWIG_GammaGammaAnalysis_H
//
// This is the declaration of the GammaGammaAnalysis class.
//

#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Analysis
 * GammaGammaAnalysis is for the analysis of events with a pair of hard
 * photons produced.  These are selected as the two highest pt photons
 * in the final state of the event.  A topdrawer file with histograms
 * is written to the working directory.  
 *
 * @see \ref GammaGammaAnalysisInterfaces "The interfaces"
 * defined for GammaGammaAnalysis.
 */
class GammaGammaAnalysis: public AnalysisHandler {

public:

  /**
   * The default constructor.
   */
  GammaGammaAnalysis();

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
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class without persistent data.
   */
  static NoPIOClassDescription<GammaGammaAnalysis> initGammaGammaAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GammaGammaAnalysis & operator=(const GammaGammaAnalysis &);

private:

  /**
   *   \f$p_T\f$ of the harder photon
   */
  Histogram _ptharder;

  /**
   *   \f$p_T\f$ of the softer photon
   */
  Histogram _ptsofter;

  /**
   *   \f$p_T\f$ of the photon pair
   */
  Histogram _ptpair;

  /**
   *   Energy of the harder photon
   */
  Histogram _Eharder;

  /**
   *   Energy of the softer photon
   */
  Histogram _Esofter;

  /**
   *   Energy of the photon pair
   */
  Histogram _Epair;

  /**
   *  Rapidity of the harder photon
   */
  Histogram _rapharder;

  /**
   *  Rapidity of the softer photon
   */
  Histogram _rapsofter;

  /**
   *  Rapidity of the photon pair
   */
  Histogram _rappair;

  /**
   *  Azimuth of the harder photon
   */
  Histogram _phiharder;

  /**
   *  Azimuth of the softer photon
   */
  Histogram _phisofter;

  /**
   *  Azimuth of the photon pair
   */
  Histogram _deltaphi;
  
  /**
   *  Invariant mass of the pair
   */
  Histogram _mpair;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of GammaGammaAnalysis. */
template <>
struct BaseClassTrait<Herwig::GammaGammaAnalysis,1> {
  /** Typedef of the first base class of GammaGammaAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the GammaGammaAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::GammaGammaAnalysis>
  : public ClassTraitsBase<Herwig::GammaGammaAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::GammaGammaAnalysis"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the GammaGammaAnalysis class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_GammaGammaAnalysis_H */
