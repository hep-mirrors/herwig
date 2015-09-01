// -*- C++ -*-
//
// BFragmentationAnalysisHandler.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_BFragmentationAnalysisHandler_H
#define HERWIG_BFragmentationAnalysisHandler_H
//
// This is the declaration of the BFragmentationAnalysisHandler class.
//

#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig/Utilities/Histogram.h"
#include "ThePEG/EventRecord/Event.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Analysis
 * The BFragmentationAnalysisHandler class is designed to compare
 * the fragmentation function for weakly decaying B hadrons with data
 * from SLD and ALEPH.
 *
 * @see \ref BFragmentationAnalysisHandlerInterfaces "The interfaces"
 * defined for BFragmentationAnalysisHandler.
 */
class BFragmentationAnalysisHandler: public AnalysisHandler {

public:

  /**
   * The default constructor.
   */
  BFragmentationAnalysisHandler() : _emax() {}

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
   *  Identifies which step(2) final state particles originate
   *  from the b/bbar...
   */
  void analyze_bquarks(ParticleSet, double);

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
  static NoPIOClassDescription<BFragmentationAnalysisHandler> initBFragmentationAnalysisHandler;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  BFragmentationAnalysisHandler & operator=(const BFragmentationAnalysisHandler &);

private:
  
  /**
   *  Histogram for the SLD binning
   */
  HistogramPtr _fragBxE;

  /**
   *  Histogram for the ALEPH binning
   */
  HistogramPtr _fragBxEa;

  /**
   *  Histograms for quark energy fraction
   */
  HistogramPtr _fragbquarkxE;

  /**
   * Histograms for b jet mass
   */
  HistogramPtr _fragbquarkjetmass;

  /**
   *  Centre-of-mass energy of the collision
   */
  Energy _emax;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of BFragmentationAnalysisHandler. */
template <>
struct BaseClassTrait<Herwig::BFragmentationAnalysisHandler,1> {
  /** Typedef of the first base class of BFragmentationAnalysisHandler. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the BFragmentationAnalysisHandler class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::BFragmentationAnalysisHandler>
  : public ClassTraitsBase<Herwig::BFragmentationAnalysisHandler> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::BFragmentationAnalysisHandler"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the BFragmentationAnalysisHandler class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwAnalysis.so HwLEPAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_BFragmentationAnalysisHandler_H */
