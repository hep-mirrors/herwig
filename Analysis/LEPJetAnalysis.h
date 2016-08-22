// -*- C++ -*-
//
// LEPJetAnalysis.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_LEPJetAnalysis_H
#define HERWIG_LEPJetAnalysis_H
//
// This is the declaration of the LEPJetAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Analysis
 * The LEPJetAnalysis class compares the results of Heriwg++ with LEP data for
 * various jet distributions.
 *
 * @see \ref LEPJetAnalysisInterfaces "The interfaces"
 * defined for LEPJetAnalysis.
 */
class LEPJetAnalysis: public AnalysisHandler {

public:

  /// Default constructor
  LEPJetAnalysis() : _nevent() {}

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
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<LEPJetAnalysis> initLEPJetAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LEPJetAnalysis & operator=(const LEPJetAnalysis &);

private:

  /**
   *  Histograms ofr the \f$y\f$ distributions
   */
  //@{
  /**
   *  \f$y_{23}\f$
   */
  HistogramPtr _y23;

  /**
   *  \f$y_{34}\f$
   */
  HistogramPtr _y34;

  /**
   *  \f$y_{45}\f$
   */
  HistogramPtr _y45;

  /**
   *  \f$y_{56}\f$
   */
  HistogramPtr _y56;
  //@}

  /**
   *  Bins for the y fractions
   */
  vector<double> _yc_frac;

  /**
   *  Points for the y fractions
   */
  //@{
  /**
   *   1 jet
   */
  vector<int> _frac1;

  /**
   *   2 jet
   */
  vector<int> _frac2;

  /**
   *   3 jet
   */
  vector<int> _frac3;

  /**
   *   4 jet
   */
  vector<int> _frac4;

  /**
   *   5 jet
   */
  vector<int> _frac5;

  /**
   *   6 jet
   */
  vector<int> _frac6;
  //@}

  /**
   *  Number of events analysed
   */
  unsigned int _nevent;

  /**
   *  N jet distribution
   */
  vector<Statistic> _njet;

  /**
   *  For different jet rates
   */
  //@{
  /**
   *  Differential two jet rate
   */
  vector<double> _d2dbins;

  /**
   *  Differential three jet rate
   */
  vector<double> _d3dbins;

  /**
   *  Differential four jet rate
   */
  vector<double> _d4dbins;

  /**
   *  differential 2->2
   */
  vector<int> _d2dN2;

  /**
   *  differential 3->2
   */
  vector<int> _d3dN2;

  /**
   *  differential 3->2
   */
  vector<int> _d3dN3;

  /**
   *  differential 4->2
   */
  vector<int> _d4dN2;

  /**
   *  differential 4->3
   */
  vector<int> _d4dN3;

  /**
   *  differential 4->4
   */
  vector<int> _d4dN4;

  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of LEPJetAnalysis. */
template <>
struct BaseClassTrait<Herwig::LEPJetAnalysis,1> {
  /** Typedef of the first base class of LEPJetAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the LEPJetAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::LEPJetAnalysis>
  : public ClassTraitsBase<Herwig::LEPJetAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::LEPJetAnalysis"; }
  /**
   * The name of a file containing the dynamic library where the class
   * LEPJetAnalysis is implemented. It may also include several, space-separated,
   * libraries if the class LEPJetAnalysis depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwLEPJetAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_LEPJetAnalysis_H */
