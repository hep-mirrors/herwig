// -*- C++ -*-
#ifndef HERWIG_BeamCorrelationAnalysis_H
#define HERWIG_BeamCorrelationAnalysis_H
//
// This is the declaration of the BeamCorrelationAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The documentation of the BeamCorrelationAnalysis class
 * is designed to plot the angle in the lab frame of an outgoing particle
 * w.r.t to an incoming particle to test the spin correlations.
 *
 * @see \ref BeamCorrelationAnalysisInterfaces "The interfaces"
 * defined for BeamCorrelationAnalysis.
 */
class BeamCorrelationAnalysis: public AnalysisHandler {

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
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<BeamCorrelationAnalysis> initBeamCorrelationAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  BeamCorrelationAnalysis & operator=(const BeamCorrelationAnalysis &) = delete;

private:

  /**
   *  PDG code of the outgoing particle
   */
  vector<int> _outgoing;

  /**
   *  PDG code of the incoming particle
   */
  vector<int> _incoming;

  /**
   *  Histograms
   */
  vector<HistogramPtr> _angle;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of BeamCorrelationAnalysis. */
template <>
struct BaseClassTrait<Herwig::BeamCorrelationAnalysis,1> {
  /** Typedef of the first base class of BeamCorrelationAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the BeamCorrelationAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::BeamCorrelationAnalysis>
  : public ClassTraitsBase<Herwig::BeamCorrelationAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::BeamCorrelationAnalysis"; }
  /**
   * The name of a file containing the dynamic library where the class
   * BeamCorrelationAnalysis is implemented. It may also include several, space-separated,
   * libraries if the class BeamCorrelationAnalysis depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwDecayAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_BeamCorrelationAnalysis_H */
