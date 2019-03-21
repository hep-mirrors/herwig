// -*- C++ -*-
#ifndef RADIATIVEZPRIME_SimpleZPrimeAnalysis_H
#define RADIATIVEZPRIME_SimpleZPrimeAnalysis_H
//
// This is the declaration of the SimpleZPrimeAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig/Utilities/Histogram.h"

namespace RadiativeZPrime {

using namespace ThePEG;
using Herwig::HistogramPtr;

/**
 * Here is the documentation of the SimpleZPrimeAnalysis class.
 *
 * @see \ref SimpleZPrimeAnalysisInterfaces "The interfaces"
 * defined for SimpleZPrimeAnalysis.
 */
class SimpleZPrimeAnalysis: public AnalysisHandler {

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
  static NoPIOClassDescription<SimpleZPrimeAnalysis> initSimpleZPrimeAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SimpleZPrimeAnalysis & operator=(const SimpleZPrimeAnalysis &) = delete;

private:

  /**
   *  Histogram for the Z \f$p_T\f$
   */
  HistogramPtr _ptZ;

  /**
   *  Histogram for the Z rapidity
   */
  HistogramPtr _rapZ;

  /**
   *  Histogram for the Z azimuth
   */
  HistogramPtr _phiZ;

  /**
   *  Histogram for the photon \f$p_T\f$
   */
  HistogramPtr _ptgamma;

  /**
   *  Histogram for the photon rapidity
   */
  HistogramPtr _rapgamma;

  /**
   *  Histogram for the photon azimuth
   */
  HistogramPtr _phigamma;

  /**
   *  Histogram for the \f$e^+\f$ \f$p_T\f$
   */
  HistogramPtr _ptep;

  /**
   *  Histogram for the \f$e^+\f$ rapidity
   */
  HistogramPtr _rapep;

  /**
   *  Histogram for the \f$e^+\f$ azimuth
   */
  HistogramPtr _phiep;

  /**
   *  Histogram for the \f$e^-\f$ \f$p_T\f$
   */
  HistogramPtr _ptem;

  /**
   *  Histogram for the \f$e^-\f$ rapidity
   */
  HistogramPtr _rapem;

  /**
   *  Histogram for the \f$e^-\f$ azimuth
   */
  HistogramPtr _phiem;

  /**
   *  Histogram for the mass of the photon and Z
   */
  HistogramPtr _masspair;

  /**
   *  Histgram for the mass of the Z
   */
  HistogramPtr _massZ;

  /**
   *  Histogram for the mass of the photon and \f$e^+\f$
   */
  HistogramPtr _massgammaep;

  /**
   *  Histogram for the mass of the photon and \f$e^-\f$
   */
  HistogramPtr _massgammaem;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SimpleZPrimeAnalysis. */
template <>
struct BaseClassTrait<RadiativeZPrime::SimpleZPrimeAnalysis,1> {
  /** Typedef of the first base class of SimpleZPrimeAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SimpleZPrimeAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<RadiativeZPrime::SimpleZPrimeAnalysis>
  : public ClassTraitsBase<RadiativeZPrime::SimpleZPrimeAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "RadiativeZPrime::SimpleZPrimeAnalysis"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SimpleZPrimeAnalysis is implemented. It may also include several, space-separated,
   * libraries if the class SimpleZPrimeAnalysis depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "RadiativeZPrime.so"; }
};

/** @endcond */

}

#endif /* RADIATIVEZPRIME_SimpleZPrimeAnalysis_H */
