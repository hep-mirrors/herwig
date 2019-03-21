// -*- C++ -*-
#ifndef HERWIG_SemiLeptonicDPiAnalysis_H
#define HERWIG_SemiLeptonicDPiAnalysis_H
//
// This is the declaration of the SemiLeptonicDPiAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig/Utilities/Histogram.h"
#include "ThePEG/EventRecord/Event.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the SemiLeptonicDPiAnalysis class.
 *
 * @see \ref SemiLeptonicDPiAnalysisInterfaces "The interfaces"
 * defined for SemiLeptonicDPiAnalysis.
 */
class SemiLeptonicDPiAnalysis: public AnalysisHandler {

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

protected:

  /**
   *  Find the stable decay products
   */
  void findChildren(tPPtr,ParticleVector &);

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SemiLeptonicDPiAnalysis> initSemiLeptonicDPiAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SemiLeptonicDPiAnalysis & operator=(const SemiLeptonicDPiAnalysis &) = delete;

private:

  /**
   *  The PDG code of the incoming particle
   */
  vector<int> _incoming;

  /**
   *  The PDG code of the outgoing D meson
   */
  vector<int> _outgoingD;

  /**
   *  The PDG code of the outgoing pion
   */
  vector<int> _outgoingP;

  /**
   *  The PDG code of the outgoing lepton
   */
  vector<long> _outgoingL;

  /**
   *  The energy of the leptons
   */
  vector<HistogramPtr> _energy;

  /**
   *  The mass of the lepton-neutrino pair
   */
  vector<HistogramPtr> _scale;

  /**
   *  the mass of the D pi system
   */
  vector<HistogramPtr> _mDpi;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SemiLeptonicDPiAnalysis. */
template <>
struct BaseClassTrait<Herwig::SemiLeptonicDPiAnalysis,1> {
  /** Typedef of the first base class of SemiLeptonicDPiAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SemiLeptonicDPiAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SemiLeptonicDPiAnalysis>
  : public ClassTraitsBase<Herwig::SemiLeptonicDPiAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SemiLeptonicDPiAnalysis"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SemiLeptonicDPiAnalysis is implemented. It may also include several, space-separated,
   * libraries if the class SemiLeptonicDPiAnalysis depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwDecayAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_SemiLeptonicDPiAnalysis_H */
