// -*- C++ -*-
#ifndef HERWIG_HiggsJetAnalysis_H
#define HERWIG_HiggsJetAnalysis_H
//
// This is the declaration of the HiggsJetAnalysis class.
//

#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "HiggsJetAnalysis.fh"
#include "Herwig++/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/**
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
  inline HiggsJetAnalysis();
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

  /**
   * Transform the event to the desired Lorentz frame and return the
   * corresponding LorentzRotation.
   * @param event a pointer to the Event to be transformed.
   * @return the LorentzRotation used in the transformation.
   */
  virtual LorentzRotation transform(tEventPtr event) const;

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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class without persistent data.
   */
  static NoPIOClassDescription<HiggsJetAnalysis> initHiggsJetAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HiggsJetAnalysis & operator=(const HiggsJetAnalysis &);

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

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of HiggsJetAnalysis. */
template <>
struct BaseClassTrait<Herwig::HiggsJetAnalysis,1> {
  /** Typedef of the first base class of HiggsJetAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the HiggsJetAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::HiggsJetAnalysis>
  : public ClassTraitsBase<Herwig::HiggsJetAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::HiggsJetAnalysis"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the HiggsJetAnalysis class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwAnalysis.so"; }
};

/** @endcond */

}

#include "HiggsJetAnalysis.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "HiggsJetAnalysis.tcc"
#endif

#endif /* HERWIG_HiggsJetAnalysis_H */
