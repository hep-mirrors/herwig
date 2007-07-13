// -*- C++ -*-
#ifndef HERWIG_SpinCorrAnalysis_H
#define HERWIG_SpinCorrAnalysis_H
//
// This is the declaration of the SpinCorrAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig++/Utilities/Histogram.h"
#include "ThePEG/PDT/ParticleData.h"
#include "SpinCorrAnalysis.fh"
#include "ThePEG/Repository/CurrentGenerator.h"

namespace Herwig {
using namespace ThePEG;
/**
 * Here is the documentation of the SpinCorrAnalysis class.
 *
 * @see \ref SpinCorrAnalysisInterfaces "The interfaces"
 * defined for SpinCorrAnalysis.
 */
class SpinCorrAnalysis: public AnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline SpinCorrAnalysis();

  /**
   * The destructor.
   */
  virtual ~SpinCorrAnalysis();
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

  /**
   *  Find the decay particle
   */
  tPPtr findDecayingParticle(tPPtr) const;

  /**
   *  Find the child
   */
  tPPtr findChild(tPPtr, long) const;

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

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);
  
  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  inline virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SpinCorrAnalysis> initSpinCorrAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SpinCorrAnalysis & operator=(const SpinCorrAnalysis &);

private:
  
  /**
   * Mass dist. quark-lepton
   */
  HistogramPtr thelPlus;

  /**
   * Mass dist quark-lepton
   */
  HistogramPtr thelMinus;

  /**
   * Maximum invariant mass of UED Spectrum
   */
  Energy theMaxInvMass;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SpinCorrAnalysis. */
template <>
struct BaseClassTrait<Herwig::SpinCorrAnalysis,1> {
  /** Typedef of the first base class of SpinCorrAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SpinCorrAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SpinCorrAnalysis>
  : public ClassTraitsBase<Herwig::SpinCorrAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SpinCorrAnalysis"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SpinCorrAnalysis is implemented. It may also include several, space-separated,
   * libraries if the class SpinCorrAnalysis depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwBSMAnalysis.so"; }
};

/** @endcond */

}

#include "SpinCorrAnalysis.icc"

#endif /* THEPEG_SpinCorrAnalysis_H */
