// -*- C++ -*-
//
// BELLECharmAnalysis.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_BELLECharmAnalysis_H
#define HERWIG_BELLECharmAnalysis_H
//
// This is the declaration of the BELLECharmAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Analysis
 * The BELLECharmAnalysis class is designed to compared the results of
 * Herwig at 10.52 GeV with data from the BELLE experiment.
 *
 * @see \ref BELLECharmAnalysisInterfaces "The interfaces"
 * defined for BELLECharmAnalysis.
 */
class BELLECharmAnalysis: public AnalysisHandler {

public:

  /**
   * The default constructor.
   */
  BELLECharmAnalysis() : _s(), _onshell(false), _ratioDstar(),
			 _ratioDs(), _ratioLambda(), _weight()
  {}

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
   * Analyze the given vector of particles. The default version calls
   * analyze(tPPtr) for each of the particles.
   * @param particles the vector of pointers to particles to be analyzed
   */
  virtual void analyze(const tPVector & particles);

  /**
   * Analyze the given particle.
   * @param particle pointer to the particle to be analyzed.
   */
  virtual void analyze(tPPtr particle);
  //@}

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
  static ClassDescription<BELLECharmAnalysis> initBELLECharmAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  BELLECharmAnalysis & operator=(const BELLECharmAnalysis &);

private:

  /**
   *  CMS energy squared
   */
  Energy2 _s;

  /**
   *  Whether this is on-shell or continuum
   */
  bool _onshell;

  /**
   *  Histogram for \f$D^{*+}\f$
   */
  HistogramPtr _histDstarplus;

  /**
   *  Histogram for \f$D^{*0}\f$
   */
  HistogramPtr _histDstar0;

  /**
   *  Histogram for \f$D^0\f$
   */
  HistogramPtr _histD0;

  /**
   *  Histogram for \f$D^+\f$
   */
  HistogramPtr _histDplus;

  /**
   *  Histogram for \f$D^+_s\f$
   */
  HistogramPtr _histDs;

  /**
   *  Histogram for \f$\Lambda_c^+\f$
   */
  HistogramPtr _histLambda;

  /**
   *  Ratios
   */
  //@{
  /**
   *  Ratio for \f$D^{*0}+D^{*+}\f$
   */
  double _ratioDstar;

  /**
   *  Ratio for \f$D_s\f$
   */
  double _ratioDs;

  /**
   *  Ratio for \f$\Lambda_c^+\f$
   */
  double _ratioLambda;
  //@}

  /**
   *  Statistics for the ratios
   */
  //@{
  /**
   *  Statistics for \f$D^0+D^+\f$
   */
  Statistic _statD;

  /**
   *  Statistics for  \f$D^{*0}+D^{*+}\f$
   */
  Statistic _statDstar;
  
  /**
   *  Statistics for \f$D_s\f$
   */
  Statistic _statDs;

  /**
   *  Statistics for \f$\Lambda_c^+\f$
   */
  Statistic _statLambda;
  //@}

  /**
   *  The weight for an event
   */
  double _weight;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of BELLECharmAnalysis. */
template <>
struct BaseClassTrait<Herwig::BELLECharmAnalysis,1> {
  /** Typedef of the first base class of BELLECharmAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the BELLECharmAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::BELLECharmAnalysis>
  : public ClassTraitsBase<Herwig::BELLECharmAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::BELLECharmAnalysis"; }
  /**
   * The name of a file containing the dynamic library where the class
   * BELLECharmAnalysis is implemented. It may also include several, space-separated,
   * libraries if the class BELLECharmAnalysis depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwAnalysis.so HwLEPAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_BELLECharmAnalysis_H */
