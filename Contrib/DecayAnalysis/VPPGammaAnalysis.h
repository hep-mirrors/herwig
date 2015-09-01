// -*- C++ -*-
#ifndef HERWIG_VPPGammaAnalysis_H
#define HERWIG_VPPGammaAnalysis_H
//
// This is the declaration of the VPPGammaAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig/Utilities/Histogram.h"
#include "ThePEG/PDT/EnumParticles.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the VPPGammaAnalysis class.
 *
 * @see \ref VPPGammaAnalysisInterfaces "The interfaces"
 * defined for VPPGammaAnalysis.
 */
class VPPGammaAnalysis: public AnalysisHandler {

public:

  /**
   * The default constructor.
   */
  VPPGammaAnalysis();

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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<VPPGammaAnalysis> initVPPGammaAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VPPGammaAnalysis & operator=(const VPPGammaAnalysis &);

private:

  /**
   *  id's of the vector mesons to consider
   */
  vector<long> _id;

  /**
   *  id's of the outgoing decay products
   */
  //@{
  /**
   *  First outgoing particle
   */
  vector<long> _outgoing1;

  /**
   *  Second outgoing particle
   */
  vector<long> _outgoing2;
  //@}

  /**
   *  histogram for the mass
   */
  vector<HistogramPtr> _masstotal;

  /**
   *  Histograms for the energies
   */
  /**
   *  Total photon energy
   */
  vector<HistogramPtr> _etotal;

  /**
   *  Energy of all the photons
   */
  vector<HistogramPtr> _eall;

  /**
   *  Single photon energy
   */
  vector<HistogramPtr> _esingle;
  //@}

  /**
   *  histograms for the multiplicities
   */
  vector<HistogramPtr> _nphoton;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of VPPGammaAnalysis. */
template <>
struct BaseClassTrait<Herwig::VPPGammaAnalysis,1> {
  /** Typedef of the first base class of VPPGammaAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the VPPGammaAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::VPPGammaAnalysis>
  : public ClassTraitsBase<Herwig::VPPGammaAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::VPPGammaAnalysis"; }
  /**
   * The name of a file containing the dynamic library where the class
   * VPPGammaAnalysis is implemented. It may also include several, space-separated,
   * libraries if the class VPPGammaAnalysis depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwDecayAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_VPPGammaAnalysis_H */
