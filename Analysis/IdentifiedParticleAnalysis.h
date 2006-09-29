// -*- C++ -*-
#ifndef HERWIG_IdentifiedParticleAnalysis_H
#define HERWIG_IdentifiedParticleAnalysis_H
//
// This is the declaration of the IdentifiedParticleAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "EventShapes.h"
#include "Herwig++/Utilities/Histogram.h"
#include "IdentifiedParticleAnalysis.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the IdentifiedParticleAnalysis class.
 *
 * @see \ref IdentifiedParticleAnalysisInterfaces "The interfaces"
 * defined for IdentifiedParticleAnalysis.
 */
class IdentifiedParticleAnalysis: public AnalysisHandler {

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

  /**
   *  Work out the flavour of the quarks produced
   */
  inline int getFlavour(const tPVector &);

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
  static ClassDescription<IdentifiedParticleAnalysis> initIdentifiedParticleAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  IdentifiedParticleAnalysis & operator=(const IdentifiedParticleAnalysis &);

private:

  /**
   *  Single particle spectra
   */
  //@{
  /**
   * Histogram for the \f$\xi\f$ distribution for all particles from all quarks
   */
  HistogramPtr _xpa;

  /**
   * Histogram for the \f$\xi\f$ distribution for all particles from light quarks
   */
  HistogramPtr _xpl;

  /**
   * Histogram for the \f$\xi\f$ distribution for all particles from charm quarks
   */
  HistogramPtr _xpc;

  /**
   * Histogram for the \f$\xi\f$ distribution for all particles from bottom quarks
   */
  HistogramPtr _xpb;

  /**
   * Histogram for the \f$\xi\f$ distribution for charged pions from all quarks
   */
  HistogramPtr _pipma;

  /**
   * Histogram for the \f$\xi\f$ distribution for charged pions from light quarks
   */
  HistogramPtr _pipml;

  /**
   * Histogram for the \f$\xi\f$ distribution for charged pions from charm quarks
   */
  HistogramPtr _pipmc;

  /**
   * Histogram for the \f$\xi\f$ distribution for charged pions from bottom quarks
   */
  HistogramPtr _pipmb;

  /**
   * Histogram for the \f$\xi\f$ distribution for charged pions from OPAL
   */
  HistogramPtr _pipm;

  /**
   * Histogram for the \f$\xi\f$ distribution for charged kaons from all quarks
   */
  HistogramPtr _kpma;

  /**
   * Histogram for the \f$\xi\f$ distribution for charged kaons from light quarks
   */
  HistogramPtr _kpml;

  /**
   * Histogram for the \f$\xi\f$ distribution for charged kaons from charm quarks
   */
  HistogramPtr _kpmc;

  /**
   * Histogram for the \f$\xi\f$ distribution for charged kaons from bottom quarks
   */
  HistogramPtr _kpmb;

  /**
   * Histogram for the \f$\xi\f$ distribution for charged kaons from OPAL
   */
  HistogramPtr _kpm;


  /**
   * Histogram for the \f$\xi\f$ distribution for protons from all quarks
   */
  HistogramPtr _ppma;

  /**
   * Histogram for the \f$\xi\f$ distribution for protons from light quarks
   */
  HistogramPtr _ppml;

  /**
   * Histogram for the \f$\xi\f$ distribution for protons from charm quarks
   */
  HistogramPtr _ppmc;

  /**
   * Histogram for the \f$\xi\f$ distribution for protons from bottom quarks
   */
  HistogramPtr _ppmb;

  /**
   * Histogram for the \f$\xi\f$ distribution for protons from OPAL
   */
  HistogramPtr _ppm;

  /**
   * Histogram for the \f$x\f$ distribution for light quark events (lin)
   */ 
  HistogramPtr _udsxp;

  /**
   * Histogram for the \f$\xi\f$ distribution for light quark events (lin)
   */ 
  HistogramPtr _udsxip;

  /**
   *  Histogram for the \f$\xi\f$ distribution for \f$\Lambda\f$ 
   */
  HistogramPtr _lpm;

  /**
   *  Pointer to the event shapes object
   */
  EventShapesPtr _shapes;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of IdentifiedParticleAnalysis. */
template <>
struct BaseClassTrait<Herwig::IdentifiedParticleAnalysis,1> {
  /** Typedef of the first base class of IdentifiedParticleAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the IdentifiedParticleAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::IdentifiedParticleAnalysis>
  : public ClassTraitsBase<Herwig::IdentifiedParticleAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::IdentifiedParticleAnalysis"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the IdentifiedParticleAnalysis class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwAnalysis.so HwLEPAnalysis.so"; }
};

/** @endcond */

}

#include "IdentifiedParticleAnalysis.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "IdentifiedParticleAnalysis.tcc"
#endif

#endif /* HERWIG_IdentifiedParticleAnalysis_H */
