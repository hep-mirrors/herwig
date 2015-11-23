// -*- C++ -*-
//
// TTbarAnalysis.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_TTbarAnalysis_H
#define HERWIG_TTbarAnalysis_H
//
// This is the declaration of the TTbarAnalysis class.
//
       
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Analysis
 * The TTbarAnalysis class tries to find a top antitop pair in the final
 * state and books a number of histograms.  It only makes sense if
 * hadronization and decays are switched off.  However, if there is is
 * no top quark pair in the final state then a warning is printed and
 * nothing is booked.  Some of the histograms will be sensitive to the
 * initial state shower.
 *
 * @see \ref TTbarAnalysisInterfaces "The interfaces"
 * defined for TTbarAnalysis.
 */
class TTbarAnalysis: public AnalysisHandler {

public:

  /**
   * The default constructor.
   */
  TTbarAnalysis();

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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class without persistent data.
   */
  static NoPIOClassDescription<TTbarAnalysis> initTTbarAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  TTbarAnalysis & operator=(const TTbarAnalysis &);

private:
  /**
   *   \f$p_T\f$ of the tops
   */
  Histogram _pttop;
  Histogram _pttbar;
  Histogram _ptpair;

  /**
   *   \f$E_T\f$ of the tops
   */
  Histogram _ettop;
  Histogram _ettbar;
  Histogram _etpair;

  /**
   *   Energy of the tops
   */
  Histogram _etop;
  Histogram _etbar;
  Histogram _epair;

  /**
   *  Rapidity of the tops
   */
  Histogram _raptop;
  Histogram _raptbar;
  Histogram _rappair;

  /**
   *  Azimuth of the tops
   */
  Histogram _phitop;
  Histogram _phitbar;
  Histogram _deltaphi;
  
  /**
   *  Invariant mass of the pair
   */
  Histogram _mpair;

  /**
   *  scalar sums of Et, pt
   */
  Histogram _etsum;
  Histogram _ptsum;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of TTbarAnalysis. */
template <>
struct BaseClassTrait<Herwig::TTbarAnalysis,1> {
  /** Typedef of the first base class of TTbarAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the TTbarAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::TTbarAnalysis>
  : public ClassTraitsBase<Herwig::TTbarAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::TTbarAnalysis"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the TTbarAnalysis class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_TTbarAnalysis_H */
