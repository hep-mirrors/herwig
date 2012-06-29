// -*- C++ -*-
#ifndef HERWIG_Upsilon4SSpectrumAnalysis_H
#define HERWIG_Upsilon4SSpectrumAnalysis_H
//
// This is the declaration of the Upsilon4SSpectrumAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig++/Utilities/Histogram.h"

namespace Herwig {
using namespace ThePEG;

/**
 * Here is the documentation of the Upsilon4SSpectrumAnalysis class.
 *
 * @see \ref Upsilon4SSpectrumAnalysisInterfaces "The interfaces"
 * defined for Upsilon4SSpectrumAnalysis.
 */
class Upsilon4SSpectrumAnalysis: public AnalysisHandler {

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
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
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
  static NoPIOClassDescription<Upsilon4SSpectrumAnalysis> initUpsilon4SSpectrumAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Upsilon4SSpectrumAnalysis & operator=(const Upsilon4SSpectrumAnalysis &);

private:

  /**
   * The \f$\pi^\pm\f$ spectrum without \f$K^0_S\f$ and \f$\Lambda^0\f$ decay products
   */
  HistogramPtr _spectrumpipNo;

  /**
   * The \f$\pi^\pm\f$ spectrum with \f$K^0_S\f$ and \f$\Lambda^0\f$ decay products
   */
  HistogramPtr _spectrumpipAll;

  /**
   *  The \f$\pi^0\f$ spectrum
   */
  HistogramPtr _spectrumpi0;

  /**
   *  The \f$K^\pm\f$ spectrum
   */
  HistogramPtr _spectrumKpA,_spectrumKpB;

  /**
   * The proton spectrum without \f$\Lambda^0\f$ decay products
   */
  HistogramPtr _spectrumprotonNo;

  /**
   * The proton spectrum with \f$\Lambda^0\f$ decay products
   */
  HistogramPtr _spectrumprotonAll;

  /**
   *  \f$K_0^S\f$ spectrum
   */
  HistogramPtr _spectrumK0;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of Upsilon4SSpectrumAnalysis. */
template <>
struct BaseClassTrait<Herwig::Upsilon4SSpectrumAnalysis,1> {
  /** Typedef of the first base class of Upsilon4SSpectrumAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the Upsilon4SSpectrumAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::Upsilon4SSpectrumAnalysis>
  : public ClassTraitsBase<Herwig::Upsilon4SSpectrumAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::Upsilon4SSpectrumAnalysis"; }
  /**
   * The name of a file containing the dynamic library where the class
   * Upsilon4SSpectrumAnalysis is implemented. It may also include several, space-separated,
   * libraries if the class Upsilon4SSpectrumAnalysis depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwDecayAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_Upsilon4SSpectrumAnalysis_H */
