// -*- C++ -*-
#ifndef HERWIG_TauCorrelationAnalysis_H
#define HERWIG_TauCorrelationAnalysis_H
//
// This is the declaration of the TauCorrelationAnalysis class.
//

#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The documentation of the TauCorrelationAnalysis class is designed to 
 * perform some analysis of the distributions of tau decay products in
 * Higgs decays.
 *
 * The analysis of the pion decays is based on hep-ph/0202007 and the 
 * rho decays is based on hep-ph/0204292.
 *
 * @see \ref TauCorrelationAnalysisInterfaces "The interfaces"
 * defined for TauCorrelationAnalysis.
 */
class TauCorrelationAnalysis: public AnalysisHandler {

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

protected:

  /**
   *  Methods to perform the analysis
   */
  //@{
  /**
   * Analyze the given particle for correlations with pi.
   * @param particle pointer to the particle to be analyzed.
   */
  void analyzePi(tPPtr particle, ParticleVector children);
  /**
   * Analyze the given particle for correlations with rho.
   * @param particle pointer to the particle to be analyzed.
   */
  void analyzeRho(tPPtr particle, ParticleVector children);
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<TauCorrelationAnalysis> initTauCorrelationAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  TauCorrelationAnalysis & operator=(const TauCorrelationAnalysis &);

private:

  /**
   *  Histogram of the \f$\phi\f$ angle in 
   *  \f$H\to\tau^+\tau^-\to\pi^+\bar{\nu}_tau\pi^-\nu_tau\f$
   */
  HistogramPtr _phi;

  /**
   * Histogram of the \f$\delta\f$ angle in 
   *  \f$H\to\tau^+\tau^-\to\pi^+\bar{\nu}_tau\pi^-\nu_tau\f$
   */
  HistogramPtr _delta;

  /**
   * Histogram of the \f$\phi\f$ angle in 
   *  \f$H\to\tau^+\tau^-\to\rho^+\bar{\nu}_tau\rho^-\nu_tau\f$
   */
  //@{
  /**
   *  First angle
   */
  HistogramPtr _rhoangle1;

  /**
   *  Second angle
   */
  HistogramPtr _rhoangle2;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of TauCorrelationAnalysis. */
template <>
struct BaseClassTrait<Herwig::TauCorrelationAnalysis,1> {
  /** Typedef of the first base class of TauCorrelationAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the TauCorrelationAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::TauCorrelationAnalysis>
  : public ClassTraitsBase<Herwig::TauCorrelationAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::TauCorrelationAnalysis"; }
  /**
   * The name of a file containing the dynamic library where the class
   * TauCorrelationAnalysis is implemented. It may also include several, space-separated,
   * libraries if the class TauCorrelationAnalysis depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwTauAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_TauCorrelationAnalysis_H */
