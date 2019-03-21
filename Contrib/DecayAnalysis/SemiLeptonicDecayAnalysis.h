// -*- C++ -*-
#ifndef HERWIG_SemiLeptonicDecayAnalysis_H
#define HERWIG_SemiLeptonicDecayAnalysis_H
//
// This is the declaration of the SemiLeptonicDecayAnalysis class.
//

#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the SemiLeptonicDecayAnalysis class.
 *
 * @see \ref SemiLeptonicDecayAnalysisInterfaces "The interfaces"
 * defined for SemiLeptonicDecayAnalysis.
 */
class SemiLeptonicDecayAnalysis: public AnalysisHandler {

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
  static NoPIOClassDescription<SemiLeptonicDecayAnalysis> initSemiLeptonicDecayAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SemiLeptonicDecayAnalysis & operator=(const SemiLeptonicDecayAnalysis &) = delete;

private:

  /**
   *  PDG codes of the decaying mesons
   */ 
  vector<long> _incoming;

  /**
   *  PDG codes of the decay products
   */
  vector<long> _outgoing;  

  /**
   *  Identidies of the leptons
   */
  vector<long> _outgoingL;

  /**
   *  Histograms
   */
  //@{
  /**
   *  The lepton energy
   */
  vector<HistogramPtr> _energy;

  /**
   *  The \f$q\f$ value
   */
  vector<HistogramPtr> _scale;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SemiLeptonicDecayAnalysis. */
template <>
struct BaseClassTrait<Herwig::SemiLeptonicDecayAnalysis,1> {
  /** Typedef of the first base class of SemiLeptonicDecayAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SemiLeptonicDecayAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SemiLeptonicDecayAnalysis>
  : public ClassTraitsBase<Herwig::SemiLeptonicDecayAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SemiLeptonicDecayAnalysis"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SemiLeptonicDecayAnalysis is implemented. It may also include several, space-separated,
   * libraries if the class SemiLeptonicDecayAnalysis depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwDecayAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_SemiLeptonicDecayAnalysis_H */
